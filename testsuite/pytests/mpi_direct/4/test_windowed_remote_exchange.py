# -*- coding: utf-8 -*-
#
# test_windowed_remote_exchange.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
EXPERIMENTAL, NOT STOCK NEST: tests of the windowed remote exchange (kernel attribute ``remote_exchange_slices``).

Four processes with two threads each. Connections inside a process are as short as one slice; every connection
between processes is at least as long as the largest window tested. The weights are whole numbers of pA, so every sum
of them is exact in any order, and the spike trains must be identical for every window. The simulation is split into
runs whose lengths are not multiples of any window, and, with two-step slices, not of the slice either.
"""

import nest
import numpy as np
import pytest

NUM_PROCESSES = 4
THREADS = 2
N = 64
INDEGREE = 10
WINDOWS = [2, 3, 7]
RUNS_MS = [23.3, 76.7, 100.0]


def node_process(node_id):
    """The process of a node, as NEST's MPIManager::get_process_id_of_node_id computes it."""
    return node_id % (NUM_PROCESSES * THREADS) % NUM_PROCESSES


def connectivity(min_delay_steps, cross_extra_steps=0, short_cross_steps=None):
    """Sources, targets, weights (pA) and delays (steps), the same on every process (fixed seed)."""
    rng = np.random.default_rng(20261028)
    targets = np.repeat(np.arange(1, N + 1), INDEGREE)
    sources = rng.integers(1, N + 1, size=targets.size)
    weights = rng.choice([60.0, 60.0, 60.0, -150.0], size=targets.size)
    same = np.array([node_process(s) == node_process(t) for s, t in zip(sources, targets)])
    delays = np.where(
        same,
        min_delay_steps + rng.integers(0, 3, size=targets.size),
        max(WINDOWS) * min_delay_steps + rng.integers(0, 5, size=targets.size) + cross_extra_steps,
    )
    if short_cross_steps is not None:
        delays[np.flatnonzero(~same)[0]] = short_cross_steps
    return sources, targets, weights, delays


def build(window, min_delay_steps=1, compressed=True, cross_extra_steps=0, short_cross_steps=None, plastic=False):
    """The test network; returns its spike recorder."""
    nest.ResetKernel()
    nest.set(
        resolution=0.1,
        local_num_threads=THREADS,
        rng_seed=12,
        use_compressed_spikes=compressed,
        remote_exchange_slices=window,
    )
    assert nest.num_processes == NUM_PROCESSES
    neurons = nest.Create("iaf_psc_exp", N)
    drive = nest.Create("poisson_generator", params={"rate": 6000.0})
    recorder = nest.Create("spike_recorder")
    nest.Connect(drive, neurons, syn_spec={"weight": 30.0, "delay": 0.1 * min_delay_steps})
    nest.Connect(neurons, recorder)
    sources, targets, weights, delays = connectivity(min_delay_steps, cross_extra_steps, short_cross_steps)
    ids = np.array(neurons.tolist())
    assert np.array_equal(ids, np.arange(1, N + 1))  # node_process() assumes the neurons are nodes 1 to N
    if plastic:
        syn_spec = {"synapse_model": "stdp_synapse", "weight": np.abs(weights), "Wmax": 1000.0, "delay": 0.1 * delays}
    else:
        syn_spec = {"synapse_model": "static_synapse", "weight": weights, "delay": 0.1 * delays}
    # node id arrays with repeated ids: NEST's array form of one-to-one Connect
    nest.Connect(ids[sources - 1], ids[targets - 1], "one_to_one", syn_spec)
    return recorder


def run(recorder):
    """Simulates in runs whose lengths are not multiples of any window; returns this process's spikes, sorted."""
    for t in RUNS_MS:
        nest.Simulate(t)
    events = recorder.get("events")
    order = np.lexsort((events["senders"], events["times"]))
    return events["times"][order], events["senders"][order]


def simulate(*args, **kwargs):
    return run(build(*args, **kwargs))


@pytest.mark.parametrize("min_delay_steps", [1, 2])
@pytest.mark.parametrize("compressed", [True, False])
def test_spike_trains_identical_for_every_window(min_delay_steps, compressed):
    reference_times, reference_senders = simulate(1, min_delay_steps, compressed)
    assert reference_times.size > 50  # this process's neurons spike enough for the comparison to mean something
    for window in WINDOWS:
        times, senders = simulate(window, min_delay_steps, compressed)
        np.testing.assert_array_equal(times, reference_times)
        np.testing.assert_array_equal(senders, reference_senders)


def test_one_step_later_across_processes_changes_the_spikes():
    """The comparison above would see a window delivered one step late: such a shift changes the spike trains."""
    reference_times, reference_senders = simulate(1)
    shifted_times, shifted_senders = simulate(1, cross_extra_steps=1)
    assert not (
        np.array_equal(reference_times, shifted_times) and np.array_equal(reference_senders, shifted_senders)
    )


def test_connection_between_processes_shorter_than_window_stops_the_run():
    build(3, short_cross_steps=2)
    with pytest.raises(nest.NESTError):
        nest.Simulate(10.0)
    build(2, short_cross_steps=2)
    nest.Simulate(10.0)


def test_window_fixed_after_first_simulation():
    simulate(2)
    with pytest.raises(nest.NESTError):
        nest.set(remote_exchange_slices=3)
    nest.set(remote_exchange_slices=2)


def test_plastic_synapses_refused():
    build(2, plastic=True)
    with pytest.raises(nest.NESTError):
        nest.Simulate(10.0)


def test_window_beyond_lag_field_refused():
    build(20000)
    with pytest.raises(nest.NESTError):
        nest.Simulate(10.0)
