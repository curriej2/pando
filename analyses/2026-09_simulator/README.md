# simulator — findings

Nothing has been built or run. `src/`, `figures/` and `results/` are empty by intent: the design
below was settled on paper and checked by Monte Carlo, but no analysis script has been written
because the approval gate has not been cleared for one.

See `CLAUDE.md` in this directory for the full design, and `notes/sciphy_notes.md` §S3 for the
theory record (SciPhy's own simulation design; the birth–death sampling model derived from the
generating function, with every closed form checked against simulation).

## What is established

- **SciPhy's simulation recipe**, read from Methods p. 13 — including that they already simulate
  both of our dropout axes (heritable silencing + sequencing dropout) but **filter** rather than
  model them, which is the row-A6 gap this project exists to close.
- **Language: Python**, because their optimisation is in the likelihood, not the simulator.
- **Tree shape and branching times factorise** — shape is parameter-free (Yule/Harding); $b$,
  $\delta$ and $\rho$ act only on the times. Verified across turnover 0 → 0.75.
- **Clone sizes cannot come from one homogeneous birth–death** — so take them from the data and
  condition the within-clone tree on observed size.
- **The branching-time CDF is closed-form and analytically invertible**, which is what makes an
  exact $O(n)$ sampler possible at 10,997 tips.

## What is owed before anything is built

1. Settle whether the branching-time density check passes, by self-calibration rather than a raw KS
   p-value (see `CLAUDE.md`, "NEXT" item 1).
2. A proposal, under the approval gate, for each script.
