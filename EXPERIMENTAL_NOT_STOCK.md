# EXPERIMENTAL, NOT STOCK NEST 3.10

A marked copy of NEST (branch `EXPERIMENTAL-NOT-STOCK-async-exchange`, on the reference commit 95e45a0ba of
`bennch-nest-simulator`) for exact asynchronous spike exchange: the fused half-step sweep (each iteration runs one
step's input half and the next step's decision half, so one step's spikes travel while the next step is computed),
then point-to-point exchange with counted completion. Switched on per run by the environment variable `NEST_ASYNC`
(`fused`: the sweep with NEST's own exchange; `async`: with the asynchronous exchange); without it the copy runs NEST's
own code path.

For laptop bit-for-bit checks against stock NEST (owner's approval of marked copies with laptop checks, 25 November
2026). Never pushed, never run on a cluster without the owner's word; nothing from it is NEST data. Method:
`sp4-direction-scoping/async_2026-11-24/METHOD.md`.
