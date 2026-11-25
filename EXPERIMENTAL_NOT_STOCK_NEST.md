# EXPERIMENTAL, NOT STOCK NEST 3.10

This tree is a changed copy of NEST, made for the order-certificate project
(`sp4-direction-scoping/new_realms_2026-11-25/learned_algorithm_search/work/`).

- **Parent:** stock `v3.10-cycle-phase-timers` at 95e45a0ba, from `/home/magnus/Projects/bennch-nest-simulator`.
  The stock tree was only read; it is never changed.
- **Branch:** `EXPERIMENTAL-NOT-STOCK-order-certified-step`. It has no git remotes; the owner adds one only to push.
- **VERSION:** `3.10.0-EXPERIMENTAL-NOT-STOCK-order-certified-step`, so every build and log names it.
- **Data from this tree** go only under `nest/nest_cluster_data/CHANGED_NEST_NOT_STOCK/` (branch name as the
  folder), with a marker in every run. They are never reference data, never pooled with stock data, and every table
  or figure that shows them says "EXPERIMENTAL, NOT STOCK NEST 3.10".
- **Purpose:** the step's delivery and drive paths are rewritten as plans that the order-type checker certifies as
  bit-identical to the parent at the same decomposition and seeds. Every change must pass a spike-for-spike
  comparison with the parent before any timing is reported.
