# SciPhy / ENGRAM signal-history project

**Reader:** Justin, MSK. Code and data on the MSK HPC; Claude Code via VS Code remote tunnel.

**The project:** extend the SciPhy Bayesian phylogenetic framework to infer **signalling history**
from ENGRAM + DNA Typewriter data. Secondary goal: make the inference scale.

## Working style — follow this

- **Derive step by step; do not summarise.** Unpack equations term by term.
- **Verify numerically before asserting.** Closed forms get checked against Monte Carlo.
- **The notes files are the running record**, not scratch. Append as work proceeds; keep the
  register consistent with what is already there.
- Flag corrections to earlier notes explicitly rather than silently editing conclusions.
- Commit at the end of each working session; the commit message summarises what changed.

## Repository map

```
CLAUDE.md                     this file — loaded every session
notes/
  sciphy_notes.md             literature + theory record (§0–§I). Largely complete.
  analysis_log.md             running record of analyses. Active.
analyses/
  YYYY-MM_short-name/         one self-contained directory per analysis
    CLAUDE.md                 what this analysis is, its state, its inputs
    README.md                 findings — write this as you go, not at the end
    src/  figures/  results/
src/                          shared code, graduated out of analyses/ on second use
refs/                         paper PDFs — GITIGNORED, never committed
data/                         GITIGNORED, symlinks to scratch
```

**Notation (house register, §H.6.0).** $N$ sites/tape · $k$ tapes/cell · $M$ alphabet size ·
$M_{\rm obs}$ observed distinct symbols · $\xi$ insert probabilities · $q=\sum_i\xi_i^2$ ·
$\lambda$ **editing rate** · $(b,\delta)$ birth–death. **$f$ is reserved for functions**, so the
ENGRAM extension can write $\xi_i(t)$, $\xi_i(\text{signal})$. Swept through §0–§I on 2026-08-28.

**Read `notes/sciphy_notes.md` by section — it is long.** Section map: §0 DNA Typewriter molecular
background · §1a–1c SciPhy model, transition probabilities, likelihood · §2/§D/§G inference
background · §D.4/§D.4b perfect-phylogeny route and the Park compatibility protocol · §E–§F design
notes · §H literature (§H.6 = full Mulberry & Stadler close read) · §I experimental design
(§I.6 rate heterogeneity, §I.7 ENGRAM close read).

## Analysis index

| directory | question | status |
|---|---|---|
| `2026-08_park-compatibility` | Does Park's cross-tape character compatibility support the perfect-phylogeny route (§D.4b)? | $\xi$/$q$ measured, dropout decomposed, homoplasy quantified per prefix node. **Character-set construction + compatibility check next.** |

## Already settled — do not re-derive unless asked

- DNA Typewriter mechanism at sequence level, incl. key-vs-PAM and *cis*/*trans* rate decomposition (§0).
- SciPhy's editing model, transition probabilities, likelihood, prefix/lcp structure, Felsenstein
  pruning (§1a–1c).
- Architecture decision: **point-estimate topology + full posterior on continuous parameters** (§F.3, §H.2).
- **The A1–A9 table** of SciPhy assumptions vs molecular reality — the master list of extension
  points, referenced throughout as "row A*n*".
- Mulberry & Stadler 2026, every equation (§H.6), incl. three findings not in the paper.
- Chen et al. 2024 ENGRAM (§I.7), incl. the ~7× rate deficit and the shared-tape $q$ constraint.

## Key quantitative facts to keep in hand

- $q=\sum_i\xi_i^2 \ge 1/M$ is the **sole channel** by which symbol composition affects topology.
- Feasibility inequality: reconstruction needs $\lambda\ell \gg q/(1-q)$.
- A single tape's dynamic range is $\approx N$-fold ($N=5$–6 for Typewriter — small).
- ENGRAM's published editing rate is **~7× below** Mulberry-optimal.
- On a shared tape, the induced signal symbol must stay **below ~30%** of insertions (ceiling ~40–45%).
- **Measured on Park (2026-08-28), correcting earlier estimates:**
  - $q\approx\mathbf{0.0170}$, *not* the 0.004 the notes long assumed — effective alphabet
    $1/q\approx\mathbf{57}$. Slightly **worse** than Mulberry's $q=1/64$ "high diversity" regime.
  - The real alphabet is **~100 symbols carrying 99.95% of edits**, not the design's $M=256$;
    the other ~97 observed symbols are artifacts.
  - Homoplasy needs **one** collision at the last prefix position, not $L$ — birthday threshold
    $m^{*}=\sqrt{2/q}\approx\mathbf{11}$ independent writes per prefix clade.
  - 96% of prefix nodes fall below $m^{*}$, and the **top 1% of nodes carry 65% of all homoplasy**
    — concentrated, hence removable (§D.4b step 5).

## Live threads (priority order)

1. **Park et al. data** — matrix loaded, $\xi$/$q$ measured, homoplasy quantified (§D.4b Step 0 and
   the $m$ estimation, both done). **Next: build the cross-tape character sets per clone and run the
   compatibility check + the flat-composition negative control.**
2. Simulate the joint design space $(p,\lambda,m,k,j,\ell)$ under ENGRAM-measured parameters,
   sweeping $\lambda$ upward 5–10× from published.
3. Derive the signal-share / dilution trade-off curve (§I.7.5).
4. The corrected reconstruction bound (§H.6.14).
5. Read: PATH close read (brief in notes), then Zwaans/GABI, LAML, ConvexML, Liao et al. 2024.

## HPC notes

**Cluster.** MSK **iris**, Slurm 25.11.5. `bsub`/`lsid` exist but are Slurm shims (`lsid` prints
"Slurm 25.11.5") — write Slurm, never LSF.

### Account vs. partition

- `-p/--partition` = **which pool of nodes** the job runs on.
- `-A/--account` = **which lab's bank the job is billed to**. No money moves; it feeds Slurm's
  fairshare priority (a lab that has used a lot recently queues lower) and the accounting DB.
  Private partitions gate on it — `lesliec` has `AllowAccounts=hpcadmins,lesliec`, so `-p lesliec`
  without `-A lesliec` is rejected.

**Use `-A lesliec` for everything.** It is the only account that clears every partition wanted here:
the private `lesliec` partition allows it, and `cpu`/`cpushort`/`gpu`/`gpushort`/`interactive`/
`cpu_highmem` are open to all accounts except a deny list `lesliec` is not on. (`normantm` — not
this lab, do not use. `preemptable` is locked to the `preemptable` partition, where jobs get killed;
only for cheap restartable sweeps.)

**⚠ Correction to the first draft of this section: `choij10` has no compute.** A `choij10`
*account* does exist in the Slurm DB, but this user is not associated with it — and there is
**no `choij10` partition at all**: 0 dedicated nodes, 0 CPUs, 0 GPUs. The Choi lab's allocation is
**storage only** (28 T on `/data1/choij10`). All compute comes from `lesliec` or the general
partitions.

### Partitions worth using

| partition | nodes | CPUs | RAM/node | GPUs | walltime |
|---|---|---|---|---|---|
| `lesliec` | 4 (`iscb017–020`) | 64/node, **256 total** | ~1 TB | **4× A100 per node, 16 total** | 7 d |
| `cpu` | 239 | 14,264 total | 1 TB+ | — | 7 d |
| `cpushort` | 234 | 13,944 | 1 TB+ | — | 2 h |
| `cpu_highmem` | 7 | 400 | ~4 TB | — | 7 d |
| `gpu` | 59 | 4,176 | 1–3 TB | A100 / A40 / L40S / H100 (4–8) / H200 / H200-NVL | 7 d |
| `interactive` | 15 | 1,088 | — | — | 1 d |

### Submitting — never on a login node

```bash
# interactive session
srun  -A lesliec -p interactive -c 8 --mem 64G -t 12:00:00 --pty bash

# CPU batch — lab nodes first, general cpu as overflow
sbatch -A lesliec -p lesliec,cpu -c 16 --mem 128G -t 2-00:00:00 job.sh

# GPU batch — list both; Slurm takes whichever frees first
sbatch -A lesliec -p lesliec,gpu --gres=gpu:1 -c 8 --mem 64G -t 1-00:00:00 job.sh
```

⚠ For GPU jobs keep the GRES **untyped** (`--gres=gpu:1`). `lesliec` is A100-only; the general
`gpu` partition is mixed. Typing it (`--gres=gpu:a100:1`) forfeits the general partition's
H100/L40S/H200 nodes and defeats the point of listing both.
Inspect with `sinfo -p lesliec,gpu -N -O NodeList,CPUs,Gres,StateLong` and `squeue -u $USER`.

### Environment

**Project env: `/data1/choij10/justin/envs/pando/bin/python`** — python 3.12 + numpy/scipy/
matplotlib, built 2026-08-28 from the **system** `miniforge3` module (`/admin/software/miniforge3`),
*not* normantm's install. Recreate with:
```bash
module load miniforge3/latest
mamba create -y --prefix /data1/choij10/justin/envs/pando python=3.12 numpy scipy matplotlib
```
Scripts that only stream CSVs use stdlib `csv` and run under any python.
`module avail`: R/4.3.0, python/3.8.0, gcc/12.2.0, cuda/12.0, gurobi/9.5.2, miniforge3/latest.

### Storage — ⚠ this corrects the repository map above

| path | size | free | use |
|---|---|---|---|
| `/data1/choij10/justin` | 28 T (lab) | 7.8 T | **everything: data, envs, outputs** |
| `/home/curriej2` | 100 G | **11 G — 90 % full** | dotfiles only; not envs, not data |
| `/scratch` | 2 T | **0 B — 100 % full** | unusable; no `/scratch/curriej2` exists |
| `/localscratch` | node-local | — | per-job temp only, not persistent |

⇒ `data/` symlinks point into **`/data1/choij10/justin/`**, *not* scratch. The repo map's
"symlinks to scratch" line at the top of this file is wrong for this cluster.
Purge/backup policy for `/data1`: **TBD**.

### Network

Outbound HTTPS works from compute nodes with no proxy variables set (checked: `ncbi.nlm.nih.gov`
200, `zenodo.org` 200) — data can be fetched on-cluster directly.

### Park data

**Uploaded 2026-08-28** — Justin's own copy from Jihye Park, at
`/data1/choij10/justin/pando/data/cancer_metastasis/` (627 MB, 6 CSVs; `data/` is gitignored).
Five `*_EditTable_filtered.csv` (cell × 166 tapes × 6 sites, wide) + `clonalbc_percell_hamming1_corrected.csv`
(clone assignment). 99,451 cells. Missing = the literal string `None`.
*(A lab copy of what appears to be the same experiment sits in `/data1/choij10/jihye/Cancer_Lineage/`
— `Mouse1–4`, `Subclones`. Not being used.)*
⚠ **Clone structure does not match §1629's "~75 clones × ~74 cells".** The delivered `ClonalBC`
column has **3,294 barcodes, median 7 cells, max 27,537**, five clones holding half the cells.
Clone-barcode dropout is uneven: 1.7% (Subclone), 5.9% (Initial), **44.9% (Mouse1)**. Resolve
before quoting per-clone results against the paper.

Step 0 (§D.4b) and the $m$/homoplasy quantification **done** — see `analyses/2026-08_park-compatibility`.

---

- Do not run long sessions on a login node.
- Data stays on the cluster. Never copy patient-derived or controlled-access data into notes,
  commit messages, figures, or anything that leaves the cluster.
- Paper PDFs live in `refs/` and are gitignored — do not commit publisher PDFs.
- **Derived per-node/per-cell tables are gitignored too** (`analyses/*/results/*.tsv.gz`): they
  carry clone barcodes from unpublished collaborator data. Only aggregate JSON summaries, figures,
  code and prose are committed. Reproduce the tables by rerunning `src/` on the cluster.
