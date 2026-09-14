# Analysis: dtt-mouse — a DNA Typewriter dataset with 11 tapes instead of 166

**Why this exists.** A stress test for the methods built in `2026-08_park-compatibility`.
Every statistic there — the disjoint A/B tape split, the variogram, the per-tape unanimity
scatter, the within-clone ladder — is limited by the **number of tapes**, not the number of
cells. Park gives 166 tapes × 6 sites. This gives **11 × 6**. If the methods degrade, this is
where we will see it, and it is the honest test of how far they generalise.

**⭐⭐ STATE (2026-09-14): DATA PULLED, NOTHING RUN. Justin's instruction — park it here for later.**

**The paper.** Yu, Kim, Seidel, …, **Junhong Choi**, Qiu, Shendure (2026). *In vivo reconstruction
of the cell lineage history of a developing mouse with DNA Typewriter, from zygote to late
organogenesis.* bioRxiv `10.64898/2026.07.29.741625`, CC-BY 4.0, posted 2026-07-30.
PDF at `refs/mouse_lineage_tracing.pdf` (gitignored).
⚑ **Junhong Choi is a co-author and is at MSK** — the obvious route to anything unreleased.

## The recorder, against Park

| | Park (cancer metastasis) | **this (embryo #3)** |
|---|---|---|
| cells | 99,451 (38,652 Subclone after filter) | **1,340,794 in the tree; 1,753,895 profiled** |
| integrations per cell | **166** | **11** |
| sites per tape | 6 | 6 |
| writable positions | **996** | **66** |
| symbol alphabet | ~100 carrying 99.95% | **8** (13 observed ≥0.5%) |
| clone structure | many `ClonalBC` clones | **one embryo**; 2 blastomeres (B1/B2) as replicates |
| missing rate | 21.8–40.0% | **~27% of locus calls** (but see the ⚠ below) |

⇒ the A/B split would be **5 vs 6 tapes**, against 83 vs 83 in Park. The per-tape scatter would
be **11 points**. There is **no clone layer**, so the ladder's $M_1$ and the whole "hold the clone
fixed" design need rethinking against tree depth instead.

## ⭐ The one thing this dataset has that Park does not: the insert ↔ tape link

The construct is a single cassette, **`epegRNA::TAPE-BC::circTAPE`** — the epegRNA that writes a
symbol, a 12-bp `TAPE-BC` (`NNNNNAANNNNN`), and the blank 6-unit TAPE are **physically linked**.
Park's `TargetBC` (10 nt) and `NNNN` (4 nt) were never linked, which is exactly why Direction 2
was parked there. The authors sequenced full cassettes and "recovered 10 of the 11 TAPE
integrations and showed them to collectively encode all 8 insertions present at ≥0.5%".
⚠ **That pairing table is NOT released** (not in the supplement S1–S7, not in the public repo).
⚑ But recovering it here is an **11 × 8** matching problem against Park's 166 × ~100, and the
silencing signal itself identifies it: silencing cassette $z$ kills both its TAPE **and** its
epegRNA, so symbol $s(z)$ should deplete across *all* tapes in exactly the cells missing tape $z$.

## ⚠⚠ WHAT IS PUBLIC, AND WHAT IS NOT — read before planning anything

**Pulled (2026-09-14), in `/data1/choij10/justin/pando/data/dtt_mouse/`:**
| file | size | content |
|---|---|---|
| `trees/nextcell_full_tree_dated.nwk.gz` | 19 MB | dated newick, **1,281,141 tips** |
| `trees/nextcell_backbone_tree_dated.nwk.gz` | 12 MB | dated newick, **655,701 tips** |
| `meta/cell_metadata.annotation.txt` | 150 MB | 1,753,895 rows: `cell_id, major_trajectory, celltype, UMAP_1, UMAP_2` |
| `geo/filelist.txt`, `series.html`, `df_gene.csv.gz` | small | GEO manifest |

**⚠⚠ THE TAPE CONSENSUS MATRIX IS NOT PUBLICLY AVAILABLE.** Checked all four plausible sources:
- **GEO `GSE341627`** (token `gbwpuyeqhrypxyt`) — **transcriptomes only**. `filelist.txt` shows 7
  samples, each `df_cell.csv.gz` + `gene_count.mtx.gz`. No circTAPE data.
- **NextCell** (`nextcell.pages.dev`, files on `shendure-web.gs.washington.edu/.../NextCell/`) —
  `sc_transcriptome/` only, plus the two trees.
- **GitHub `seidels/dtt-mouse-analysis`** — tree-building and downstream only. The
  `tape_pipeline/` directory is **absent**, as are `tape/merge.py`, `tape/sc.py`,
  `tape/editchain.py`, `blastomere_route.py`, the integration-barcode whitelist and the
  insertion vocabulary.
- **Supplement** — S1–S7 are blastomere edits, cell counts, sibling enrichment, couplings. No map.

⇒ **what we would need to request** (short email to Junhong Choi is likely faster than anything
else): `e3v5v6.{B1,B2}_tape_consensus.tsv.gz`, the pre-threshold per-(cell, integration) **molecule
counts**, the whitelist + vocabulary, and the epegRNA↔TAPE-BC pairing.

## ⭐ The shape we would get, documented from their own parser

`tree_building/1_build_nj_backbone/parse_tape_consensus.R` pins the contract of
`e3_tape_consensus.tsv`:
- col 1 `cell_id`; **cols 2–12 are the 11 integration barcodes, as column names**
- each value `"s1|s2|...|s6"`, or **`NA` = integration not recovered**
- site token `"U"` = unedited; anything else is an inserted symbol (`"AAG"`, `"GATG"`)
- later deliveries add `n_loci`, `n_doublet_loci`, `mean_dominance`, `pass_qc`

⇒ **this is our `dropout_matrix` / `prefix_codes6` structure**, one tape per column instead of
wide. Conversion is trivial. ⚑ **They ship `n_loci` per cell** — the per-cell capture statistic
($\alpha_c$ / $R_c$) we have to compute ourselves in Park — and `mean_dominance`, which Park has
no analogue for.

## ⚠⚠ THE METHODOLOGICAL CATCH: missingness here is partly a deliberate threshold

> "A locus is emitted only when the called lineage is supported by ≥3 molecules and its dominance
> is ≥0.9; otherwise it is **left missing, on the principle that a wrong call harms the tree while
> a missing one does not** … these thresholds retain ~98% of cells while removing **~27% of the
> lowest-confidence locus calls**."

In Park, missing = *the tape was not recovered*. Here, missing = *not recovered **or** recovered
below a confidence threshold*. Our entire claim is that missingness is partly biological rather
than technical, so a large deliberate analytic component inside it is a real complication.
**Separating them needs the pre-threshold molecule counts**, which live in the unreleased pipeline.

⚑ **They already treat shared missingness as a hazard without naming it.** Blastomere routing
bipartitions only cells with "≥8 recovered integrations, **so that shared missingness could not
drive the split**", and an opposite-side call is "treated as ambient bleed and set to missing".

## Notes for whoever picks this up

- ⚠ **NCBI does not resolve from the login nodes.** `github.com` and
  `shendure-web.gs.washington.edu` do; `ncbi.nlm.nih.gov` needs a **compute node**
  (`sbatch -A lesliec -p cpu`). `src/01_fetch.sh` is written to be run that way and is re-runnable.
- Data lives in `/data1/choij10/justin/pando/data/dtt_mouse/` (gitignored, as `data/` always is).
- Embryo #2 had ~35 integrations but was **not** profiled at single-cell resolution; all
  single-cell circTAPE data derive from **embryo #3** (11 integrations).
- The 13-vs-11 symbol excess is explained as transiently expressed epegRNAs acting shortly after
  pronuclear injection, restricted to the 5′-most sites — a characterised caveat, not a mystery.
