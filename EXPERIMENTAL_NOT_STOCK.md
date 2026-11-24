# EXPERIMENTAL, NOT STOCK NEST: one flat array for the compressed spike data

Branch `EXPERIMENTAL-NOT-STOCK-flat-compressed-spike-data`, from the reference build's `v3.10-cycle-phase-timers`
(95e45a0ba). Never for NEST's main branch; data from it go to `nest_cluster_data/CHANGED_NEST_NOT_STOCK/`.

**The change.** Stock NEST keeps, for each synapse type, one separately allocated vector of entries per sender (one
entry per thread: 64 × 8 bytes), about 2 GB per process in the multi-area model, on 4 KB pages. This copy stores the
same entries in one flat array per synapse type, indexed `source index × threads + thread`, allocated once at its final
size. Every index and every entry is the same as in stock NEST. Only where the entries lie in memory differs, so the
simulation is unchanged bit for bit (to be shown before any timing is read).

**Why.** The mechanistic model's test of an intervention (`sp4-direction-scoping/mechanistic_model_2026-11-22/mechanism/
REGISTRATION_INTERVENTION.md`). On Hambach, NEST's lookup of the receiving thread's entry missed the L2 TLB on 90 % of
received spikes, with walks of about 185 cycles (IBS, `REGISTRATION_IBS_KINDS.md`). One allocation can sit on 2 MB
pages, and each lookup also loses a dependent load (the vector's header).
