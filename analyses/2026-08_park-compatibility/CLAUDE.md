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

**State.** Step 0 done. $\xi$, $q$, dropout decomposition and the $m$-identifiability check all done.
Clone structure resolved against the paper (§D.4c) — our ClonalBC handling reproduces theirs exactly.
**Character-set construction is the next step**; compatibility check not started.

**Analysis unit (settled §D.4c).** Per `ClonalBC` clone, cells = clonal-barcode table ∩ group edit
table, **with the paper's tape filter applied (≥100 recovered tapes for Initial/Subclone, ≥20 for
Mouse1–3) — it is NOT applied in the delivered CSVs.** Clone-size floor is ours to choose, for
statistical power, not inherited from the paper.

## Conventions here
- Scripts in `src/`, figures in `figures/`, small tables in `results/`.
- Record findings in `README.md` as you go, not at the end.
- Scripts that only stream the CSVs use stdlib `csv` and run under any python.
- Numerics/figures need the project env: `/data1/choij10/justin/envs/pando/bin/python`
  (python 3.12 + numpy/scipy/matplotlib, built from the system `miniforge3` module).
