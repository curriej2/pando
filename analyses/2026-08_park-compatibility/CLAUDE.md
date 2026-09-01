# Analysis: park-compatibility

**Question.** Does the cross-tape character-compatibility structure of the Park data support the
perfect-phylogeny route (§D.4), and how large is the gap between homoplasy-predicted and measured
incompatibility — i.e. how much dropout/error/model violation is really there?

**Protocol.** `notes/sciphy_notes.md` §D.4b, executed in order. §D.4 for the underlying argument.

**Inputs.** `/data1/choij10/justin/pando/data/cancer_metastasis/` (symlink target of `data/`;
gitignored, never copy out). Justin's own copy of the Park, Chang et al. 2026 tables from Jihye Park.
- `{Initial,Mouse1,Mouse2,Mouse3,Subclone}_EditTable_filtered.csv` — cell × (166 tapes × 6 sites),
  wide, 997 columns. Values are NNNNGGA 7-mers; missing is the literal string `None`.
- `clonalbc_percell_hamming1_corrected.csv` — 167,736 rows, `CellID,Sample,ClonalBC_raw,ClonalBC`.
  **This is the clone assignment**, needed for §D.4b Procedure step 1 (work within a clone).

**State (2026-09-01).** Diagnostics complete. Compatibility measured on Mouse3 and Initial;
Mouse1/Mouse2 never finished (5 OOMs, all on their single largest clone) and are now **moot** —
see the strategic reassessment in README.md. $C$ measured on Mouse3 (0.42–0.57 depending on
heuristic). Subclone ground-truth test passed.

**⚑ Direction has changed.** The skeleton is a step sideways: its output is arbitrary (45% spread
across heuristics), only 9 of 2,547 clones exceed 1,000 cells so the scale problem it solves barely
exists here, and Felsenstein pruning handles dropout natively. **Current work is a figure programme
motivating the likelihood route**, not further skeleton development.

**Analysis unit (settled §D.4c).** Per `ClonalBC` clone, cells = clonal-barcode table ∩ group edit
table, **with the paper's tape filter applied (≥100 recovered tapes for Initial/Subclone, ≥20 for
Mouse1–3) — it is NOT applied in the delivered CSVs.** Clone-size floor is ours to choose, for
statistical power, not inherited from the paper.

## Figure programme (for PI presentations)

| fig | status | script |
|---|---|---|
| 1 recorder / $q$ | ✅ done, 2 panels | `20_fig1_recorder.py` |
| 2 homoplasy | ✅ done, **two versions** — `simple` (present) and `mle` (reserve) | `21_fig2_homoplasy.py <simple\|mle>` |
| 3 dropout | ⚠ **needs redesign** — see below | `22_fig3_dropout.py` |
| 4 compatibility spread + homoplasy null | ⭐ next; **needs the simulator** | — |
| 5 method comparison under simulation | planned; needs simulator | — |
| 6 calibration / honest uncertainty | planned; needs simulator | — |
| 7 runtime–accuracy frontier | planned | — |

**Fig 3 critique (2026-09-01, from Justin).** Panel b (determined fraction vs depth, naive vs
corrected) is *uninterpretable* — ten overlapping lines with no arm distinction — and worse, it
argues against a **strawman**: nobody proposes discarding information from unedited sites. Panel c
(per-arm dropout) is largely evident from panel a. Panel a alone is worth keeping; the decomposition
is genuinely informative (Pre-TX is 37% edits / 37% unedited / 25% absent).
⇒ Redesign so the figure shows dropout's **magnitude and structure** — is it random, per-cell,
per-tape? does it correlate with edit depth? — rather than our handling of it. Its consequence
belongs in Fig 4.

**Figure conventions.** Light theme, palette from the `dataviz` skill (surface `#fcfcfb`, ink
`#0b0b0b`/`#52514e`, grid `#e1e0d9`, series blue `#2a78d6`, orange `#eb6834`, aqua `#1baf7a`).
Minimal annotating text — Justin's explicit preference. Panels should **teach the quantity**
(show the formula / the read-off), not just report it. Always render and inspect the PNG for label
collisions before calling a figure done.

## Conventions here
- Scripts in `src/`, figures in `figures/`, small tables in `results/`.
- Record findings in `README.md` as you go, not at the end.
- Scripts that only stream the CSVs use stdlib `csv` and run under any python.
- Numerics/figures need the project env: `/data1/choij10/justin/envs/pando/bin/python`
  (python 3.12 + numpy/scipy/matplotlib, built from the system `miniforge3` module).
