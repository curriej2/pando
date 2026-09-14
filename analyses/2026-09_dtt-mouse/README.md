# dtt-mouse — findings

Nothing has been run. This records what was pulled on 2026-09-14 and what it is good for.
See `CLAUDE.md` for the full assessment.

## Purpose

A **low-tape stress test**. Park has 166 tapes; this has 11. Every statistic in
`2026-08_park-compatibility` scales with tape count, not cell count, so this is where we find out
how far they generalise.

## What was pulled

`/data1/choij10/justin/pando/data/dtt_mouse/` — dated trees (1,281,141 and 655,701 tips),
cell metadata (1,753,895 rows), and the GEO manifest. Fetched by `src/01_fetch.sh`, which must
run on a compute node because NCBI does not resolve from the login nodes.

## ⚠⚠ The tape matrix is not public

GEO `GSE341627` holds transcriptomes only; the NextCell host and the GitHub repo likewise. The
`tape_pipeline/` code, the integration-barcode whitelist, the insertion vocabulary and the
epegRNA↔TAPE-BC pairing are all absent. **The trees are the only lineage product released** — they
encode the relationships *inferred from* the tape data, but not the per-cell × per-tape genotypes
or their missingness, which is exactly what our dropout methods consume.

⇒ To use this dataset for anything beyond tree-shape work we need to request
`e3v5v6.{B1,B2}_tape_consensus.tsv.gz` and the pre-threshold molecule counts. Junhong Choi is a
co-author and is at MSK.

## What the trees alone could still support

- tree-shape statistics and clade-size distributions at 1.3 M tips
- a sanity check of our relatedness metric against *their* inferred topology, if the tape matrix
  ever arrives
- nothing about dropout, since missingness is not in the released products
