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
| `2026-08_park-compatibility` | Does Park's cross-tape character compatibility support the perfect-phylogeny route (§D.4b)? | **Diagnostics complete.** Mouse3 94.66%/63.56% (spread +31 pts), Initial 91.91%/80.60% (+11 pts) ⇒ **dropout binds, row A6**. $C=0.57$ on Mouse3. Subclone ground truth **passed**. ⚑ **Strategic pivot: skeleton is a step sideways — only 9 of 2,547 clones exceed 1,000 cells, so likelihood inference is already in range. Figure programme under way: figs 1–3 done.** ⚑ **Dropout characterised (fig 3): two axes — the tape axis is larger ($\rho_{\rm tape}=0.25$ vs $\rho_{\rm cell}=0.13$), reproducible ($r=0.997$ between libraries) and informative about edit depth, but only through a removable 15% of tapes; the cell axis is a pure QC artefact and carries no editing information, so it can be marginalised.** |

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
- **⚠ Two conventions answer different questions (§D.4d).** Missing-excluded compatibility is
  *pair-specific* and cannot be assembled into a skeleton; only missing-as-absent can. So the
  spread is the gap between the ideal and the achievable, not an error bar.
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

0. **⭐ Figure programme for PI presentations** (`analyses/2026-08_park-compatibility`, **figs 1–3
   done**; fig 3 is five standalone panels a–e, deliberately *not* composed into a grid).
   **Next: plan a figure for the A9 lineage test** — do related cells lose the same tapes? Mulberry &
   Stadler name this as their reason for punting on dropout (§1c.2) and it is unmeasured in the
   literature. Attenuated here because shelf cells largely lack clone barcodes, and its honest null
   is the simulator. **The remaining figs are blocked on that simulator**: birth–death
   tree at Park clone sizes, sequential editing at measured $\lambda$/$\xi$, $N=6$, $k=166$, with
   **dropout as a switchable layer**. That simulator also supplies the **homoplasy null** — still the
   single biggest gap in the "dropout not homoplasy" argument.
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
✅ **Clone structure resolved against the paper 2026-08-31 (§D.4c).** The "~75 clones × ~74 cells"
in the old notes was a *misreading*: it is Metient's **migration** subset — Mouse 1 only, clones with
≥10 cells that appear in ≥2 organs (93 → 76 → 75 after the collision screen; 5,551 cells). We
reproduce their 93 and 76 exactly, so our ClonalBC handling **is** their pipeline. Tree
reconstruction has no such threshold: one clone at a time, cells = clonal-barcode table ∩ group edit
table. ⚠ **Their cell filter (≥100 recovered tapes for Initial/Subclone, ≥20 for the mice) is NOT
applied in the delivered tables** — apply it ourselves; it drops 2.1% of cells. Pooled `ClonalBC`
has 3,294 barcodes, median 7 cells, max 27,537, five clones holding half the cells.

Step 0 (§D.4b) and the $m$/homoplasy quantification **done** — see `analyses/2026-08_park-compatibility`.

---

### ⚠ Right-size the request — big asks queue behind small ones

`submit.sh`'s defaults (128 G, 12 h) are a **ceiling, not a recommendation**. Oversized requests sit
in `PD` behind jobs that would have run immediately, and on a busy day that is the difference between
results tonight and results tomorrow. Two rules:

**1. Estimate before submitting.** Ask what actually drives peak memory, not what the input weighs.
For the compatibility check it is the **single largest clone**, because clones are processed one at a
time — not the table size. Measured 2026-08-31 (`sstat -j <id>.batch -o MaxRSS`):

| table | characters | largest clone | **peak RSS** | requested |
|---|---|---|---|---|
| Mouse3 | 55,875 | 2,372 chars / 210 cells | **2.1 G** | 64 G |
| Mouse1 | 157,356 | 14,640 chars / 1,607 cells | **3.4 G** | 128 G |
| Mouse2 | 75,733 | 15,533 chars / 3,387 cells | **4.2 G** | 128 G |
| Initial | **1,006,226** | 3,963 chars / 127 cells | **1.9 G** | 256 G |

⇒ Initial has **18× more characters than Mouse2 and uses less than half the memory**, because its
biggest clone is small. Sizing on table size would have been exactly wrong. Every one of these
requests was 30–134× oversized.

⚠⚠ **A mid-run `MaxRSS` is worthless when the job processes work in size order.** The table above
was read ~4 minutes into runs that sort clones by *ascending* size, so it measured the cheapest
clones and nothing else. The rule of thumb derived from it was wrong: **Mouse2 then died
`OUT_OF_MEMORY` at 128 G**, and Mouse1 climbed from 3.4 G to >20 G once it reached its large clones.
Only a **completed** job's `MaxRSS` (`sacct -j <id> -o MaxRSS`) is a size estimate; a partial reading
tells you nothing except a lower bound.

Bound script 14 analytically instead: peak is driven by `nnz(MM^T)` for the single largest clone,
`≤ Σ_cells k(k−1)/2` with k = characters containing that cell, at ~16 bytes per non-zero.

**2. If a job will not start, resubmit smaller rather than waiting.** Check why first:

```bash
squeue -u $USER -t PD -o "%.10i %.16j %.40R"      # reason: Priority / Resources / QOSMaxJob…
sinfo -p lesliec,cpu -o "%.12P %.6a %.15F %.10m"  # nodes A/I/O/T -- is anything idle?
sacct -j <id> -o JobID,MaxRSS,Elapsed,State       # what a past run of this actually used
```

- `Resources` or a long `Priority` wait with nothing idle ⇒ **cancel and resubmit at a fraction of
  the memory and walltime.** A 16 G / 2 h job backfills into gaps a 256 G / 24 h job cannot.
- Drop `-p lesliec,cpu` to just `cpu` if the four `lesliec` nodes are full — `cpu` has 239 nodes.
- For a first run on unfamiliar input, submit the **smallest input** with a small request, read
  `MaxRSS` off it, then size the rest from measurement.
- Slurm kills a job that exceeds its `--mem` (`OUT_OF_MEMORY`), so under-asking is cheap to detect
  and costs one resubmit; over-asking costs queue time silently.

- **⚠ Run every analysis as a Slurm batch job — `scripts/submit.sh <script.py> [args]`.**
  Anything long-running in the VS Code tunnel takes the tunnel *and* the Claude Code session
  down with it when it dies, losing both the computation and the transcript. This happened on
  2026-08-31 during the §D.4b step-4 run. Defaults: `-A lesliec -p lesliec,cpu -c 4 --mem 128G
  -t 12:00:00`, logs to `logs/`; override with `--mem/--time/--cpus`.
- Do not run long sessions on a login node.
- Data stays on the cluster. Never copy patient-derived or controlled-access data into notes,
  commit messages, figures, or anything that leaves the cluster.
- Paper PDFs live in `refs/` and are gitignored — do not commit publisher PDFs.
- **Derived per-node/per-cell tables are gitignored too** (`analyses/*/results/*.tsv.gz`): they
  carry clone barcodes from unpublished collaborator data. Only aggregate JSON summaries, figures,
  code and prose are committed. Reproduce the tables by rerunning `src/` on the cluster.
