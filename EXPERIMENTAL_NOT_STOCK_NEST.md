# EXPERIMENTAL NEST: NOT STOCK NEST 3.10, NOT THE REFERENCE

This tree is a separate, confined copy of NEST, made for one experiment. It is not NEST 3.10 as released, and it is
not our reference NEST. Never build it into a main NEST tree or install. Never present data from it as NEST data.

- **Experiment:** the lookahead experiment in
  `/home/magnus/Projects/sp4-direction-scoping/lookahead_percolation_2026-10-28/` (`NEST_DESIGN.md`).
- **Base:** our reference NEST, branch `v3.10-cycle-phase-timers`, commit 95e45a0ba. That is NEST 3.10 plus the
  per-phase cycle timers.
- **Branch:** `EXPERIMENTAL-NOT-STOCK-windowed-remote-exchange`.
- **Version string:** `VERSION` reads `3.10.0-EXPERIMENTAL-NOT-STOCK-windowed-remote-exchange`. So the build
  configuration, `nest.__version__` and every model log of a run name this copy.
- **Remotes:** none. The owner adds one only to push this branch.
- **Data:** every run of this copy is a run of a changed NEST. Its data go only to
  `nest/nest_cluster_data/CHANGED_NEST_NOT_STOCK/EXPERIMENTAL-NOT-STOCK-windowed-remote-exchange/`, with a marker in
  every run. They never go among the reference runs, and never into a table or figure without the label
  "EXPERIMENTAL, NOT STOCK NEST 3.10".

## The change

A new kernel setting, `remote_exchange_slices` = W, sets how often spikes are exchanged with other MPI processes. W is
a whole number of slices; one slice is `min_delay` steps. The default, 1, keeps NEST's own exchange and delivery
paths.

**With W > 1:**
- Spikes for targets on the own process are delivered at the start of the next slice, as in stock NEST, but without MPI.
- Spikes for targets on other processes are held for W slices. Their lag counts steps from the window's start.
  Then one all-to-all exchanges them, and they are delivered at the start of the slice after the window.
- **Exact when every connection between processes has a delay of at least W slices.** Such a spike reaches its
  target's ring buffer before the target is updated at its arrival step.
- **Checked before each run.** NEST computes the smallest delay of any connection whose source is on another
  process, over all processes, whenever the connections have changed. The run stops if that delay is below W slices.
- **The run also stops with W > 1 if any of these is used,** since none has been checked:
  - off-grid spiking;
  - secondary connections;
  - structural plasticity;
  - a synapse model other than `static_synapse` and its variants;
  - a source table that has been cleared.
- **Not bitwise identical to W = 1.** Spikes reach a ring-buffer slot in another order, so sums of floating-point
  weights can differ in the last bit. Spike trains are identical when weights are small whole numbers.
