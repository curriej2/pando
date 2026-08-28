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
Character-set construction is the next step; compatibility check not started.

## Conventions here
- Scripts in `src/`, figures in `figures/`, small tables in `results/`.
- Record findings in `README.md` as you go, not at the end.
- Scripts that only stream the CSVs use stdlib `csv` and run under any python.
- Numerics/figures need the project env: `/data1/choij10/justin/envs/pando/bin/python`
  (python 3.12 + numpy/scipy/matplotlib, built from the system `miniforge3` module).
