# EXPERIMENTAL, NOT STOCK NEST 3.10

**This tree is a changed NEST. Nothing built from it, and no data run with it, is stock NEST 3.10.**

- **What:** NEST's `fixed_total_number` rule with static connections, regenerated on demand instead of stored
  (Stage A of `../NEST_DESIGN.md`). The connection rule's law per virtual process is NEST's own; the realisation is
  not stock NEST's.
- **Base:** the stock ci-beNNch NEST, `Projects/bennch-nest-simulator`, branch `v3.10-cycle-phase-timers` at
  95e45a0ba (NEST 3.10 plus the per-phase cycle timers). That tree is never changed; this is a separate clone made
  without hardlinks (4 Nov 2026).
- **Branch:** `EXPERIMENTAL-NOT-STOCK-procedural-connections`. No git remotes; the owner adds one only to push.
- **Version string:** `3.10.0+experimental.not.stock.procedural.connections` (`VERSION`), so every build and log names
  it.
- **Data:** only under `nest/nest_cluster_data/CHANGED_NEST_NOT_STOCK/EXPERIMENTAL-NOT-STOCK-procedural-connections/`,
  with a marker in every run. Never among the reference runs, never pooled or tabled with stock data without the label
  "EXPERIMENTAL, NOT STOCK NEST 3.10".
- **Runs:**
  - laptop code checks only (tiny size, never data);
  - tests as the owner's 1-node Hambach exception of 29 Oct allows;
  - benchmark runs only through ci-beNNch, by the owner.
