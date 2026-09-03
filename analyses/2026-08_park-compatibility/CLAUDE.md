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

**State (2026-09-02, end of session 2).** **Fig 3 done (a–e)**; **row A9 measured** — tape loss is
heritable, heritable *below* the clone (monotone subclade-depth gradient, 5/5 arms), and discrete
rather than a graded rate. Fig 4 panels a/b built. See README "Row A9" and the full $\rho$ derivation
recorded under "Methods, in full".

**⭐ NEXT, and it is the crux.** The effect is claimed to be a **tree** property, not a clonal one,
so the sharpest test is that it appears **WITHIN a single large clone** at subclade resolution —
where a clone-level or batch effect cannot reach. Start with **Mouse2's 3,387-cell clone**: it has
near-zero *clone-level* excess yet holds 98.9% of that arm's pooled weight, so if tree structure is
real it must be visible inside it. Then Subclone's 10,996-cell clone. Panels **c and d may carry the
figure alone**; a and b become supporting.

**⚠ Estimator rules, learned the hard way.**
0. **Every layer needs its own margin.** $\Lambda$ scored against $\alpha_c+\beta_z$ alone gave
   FDR 64–73% at *every* threshold, because clone-wide losses look spectacular in permuted data too.
   Adding a per-(clone,tape) $\gamma$ dropped it to ≈0%. Whatever layer you are testing, fit the
   layer above it.
1. **Equal-clone weighting, not pooled.** $\hat\rho$ is a ratio of sums, so clones enter weighted by
   $n_C(n_C-1)$ — one clone held 98.9% of Mouse2 and 84.0% of Mouse1, and pooling *hid* the signal
   (Mouse2 +0.009 pooled vs +0.246 equal-weighted).
2. **Significance saturates; report effect size.** With millions of within-clone pairs, 77–100% of
   tapes reach FDR 5% and median $z$ runs +8 to +422. $z$ is an ordering statistic only — the
   per-tape null is right-skewed. Dropout characterised on both margins — see
README "Fig 3 redesign". Compatibility measured on Mouse3 and Initial;
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
| 3 dropout | ✅ **done, 5 panels a–e**, each a standalone PNG | `27` a · `28` b · `30` c · `31` d+e |
| 4 dropout & lineage (row **A9**) | ⚑ a/b built; **event catalogue done — c/d next** | `32`–`43` |
| 5 compatibility spread + homoplasy null | needs the simulator | — |
| 6 method comparison under simulation | planned; needs simulator | — |
| 7 calibration / honest uncertainty | planned; needs simulator | — |
| 8 runtime–accuracy frontier | planned | — |

**No composed grid.** Justin does not want a–e assembled into one multi-panel figure; the panels are
presented individually. So `27/28/30/31` each own their PNG and there is no assembly script —
do not build one.

**Fig 3 as built (2026-09-01).** The old figure described *our handling* of `None`; the redesign
measures the assay. The original decomposition panel was dropped outright — nobody disputes that
trailing `None` is biology, and if challenged the answer is the flat per-cell depth profile (panel e),
not the internal-gap rate, which only excludes *random* site loss and not terminal truncation.

| panel | claim | headline |
|---|---|---|
| a | dropout is not a coin flip on entries | VIF 22.6, $\rho_{\rm cell}=0.131$; $R_c$ bimodal |
| b | which tape it is matters more than which cell | $\rho_{\rm tape}=0.250$; rates 0.006–0.962 |
| c | the shelf is a QC choice, not an arm difference | $\rho_{\rm cell}$ 0.13–0.18 → 0.012–0.024 at a common ≥100 cut |
| d | a tape we cannot see recorded less | $\rho=+0.34$; deciles 2.92 → 4.95 sites |
| e | …but a cell we see badly recorded as much | $\rho=+0.05$; 4.88 → 4.98 sites |

d and e **must keep their shared $y$-axis** — the contrast is the argument and it dies if they are
scaled independently.

**Figure conventions.** Light theme, palette from the `dataviz` skill (surface `#fcfcfb`, ink
`#0b0b0b`/`#52514e`, grid `#e1e0d9`). **Orange `#eb6834` means "the null" figure-wide** and is not
available as a series colour. Five-arm series use categorical slots 1,3,4,5,7 —
Mouse1 `#2a78d6` · Mouse2 `#1baf7a` · Mouse3 `#eda100` · Pre-TX `#e87ba4` · Subclone `#4a3aa7`
(validated: CVD $\Delta E$ 9.1, normal-vision 19.6 on the adjacent pairlist). Where a panel needs
both five arms and a null, the **null is drawn neutral** (`#cfcec6`/`#8f8e86`), because orange fails
the normal-vision floor against magenta and yellow.
⚠ The skill's `validate_palette.js` will not run on the cluster (node v10.24: no ESM, no `??=`).
Port it to Python — same thresholds, same Machado–Oliveira–Fernandes matrices — rather than eyeballing.
Minimal annotating text — Justin's explicit preference. Panels should **teach the quantity**
(show the formula / the read-off), not just report it. Always render and inspect the PNG for label
collisions before calling a figure done.

## Conventions here
- Scripts in `src/`, figures in `figures/`, small tables in `results/`.
- Record findings in `README.md` as you go, not at the end.
- Scripts that only stream the CSVs use stdlib `csv` and run under any python.
- Numerics/figures need the project env: `/data1/choij10/justin/envs/pando/bin/python`
  (python 3.12 + numpy/scipy/matplotlib, built from the system `miniforge3` module).
