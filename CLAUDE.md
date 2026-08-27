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

**Read `notes/sciphy_notes.md` by section — it is long.** Section map: §0 DNA Typewriter molecular
background · §1a–1c SciPhy model, transition probabilities, likelihood · §2/§D/§G inference
background · §D.4/§D.4b perfect-phylogeny route and the Park compatibility protocol · §E–§F design
notes · §H literature (§H.6 = full Mulberry & Stadler close read) · §I experimental design
(§I.6 rate heterogeneity, §I.7 ENGRAM close read).

## Analysis index

| directory | question | status |
|---|---|---|
| *(add each analysis here as it starts — one line, keep current)* | | |

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

- $q=\sum_i\xi_i^2 \ge 1/j$ is the **sole channel** by which symbol composition affects topology.
- Feasibility inequality: reconstruction needs $\lambda\ell \gg q/(1-q)$.
- A single tape's dynamic range is $\approx k$-fold ($k=5$–6 for Typewriter — small).
- ENGRAM's published editing rate is **~7× below** Mulberry-optimal.
- On a shared tape, the induced signal symbol must stay **below ~30%** of insertions (ceiling ~40–45%).

## Live threads (priority order)

1. **Park et al. data** — load the per-cell × per-tape genotype matrix, reproduce one published
   summary number, then run the §D.4b compatibility check and the flat-composition negative control.
2. Simulate the joint design space $(p,\lambda,m,k,j,\ell)$ under ENGRAM-measured parameters,
   sweeping $\lambda$ upward 5–10× from published.
3. Derive the signal-share / dilution trade-off curve (§I.7.5).
4. The corrected reconstruction bound (§H.6.14).
5. Read: PATH close read (brief in notes), then Zwaans/GABI, LAML, ConvexML, Liao et al. 2024.

## HPC notes

*(fill in: scheduler and interactive-job command, module loads, conda/venv activation, scratch vs
home paths and their purge/backup policy, where the Park data lives, proxy env vars if needed)*

- Do not run long sessions on a login node.
- Data stays on the cluster. Never copy patient-derived or controlled-access data into notes,
  commit messages, figures, or anything that leaves the cluster.
- Paper PDFs live in `refs/` and are gitignored — do not commit publisher PDFs.
