# simulator — findings

**State 2026-09-20: the tree layer is BUILT and VALIDATED at both capture fractions. The cost
pilot (`03`) has not been run.** Editing and dropout are separate later layers, not in scope here.

## ✅ Validation COMPLETE (`02`, 2026-09-20/21) — PASS in every configuration

**84,000 trees, zero structure faults.** Three turnovers × three $n$ × two capture fractions, plus
four large-$n$ structural points. Worst statistic in any run sat at **0.96×** its own p99 critical
value; every other run at 0.64–0.71×.

| run | $\rho$ | $\theta$ swept | $n$ | large-$n$ (check F) | worst / crit |
|---|---|---|---|---|---|
| pass 1 | 0.5 · 8e-4 | 0.3 | 4, 5, 8 | — | 0.64 · 0.71 |
| pass 2 | 0.5 | 0, 0.3, 0.75 | 4, 5, 8 | 1,024 · 3,387 | 0.96 |
| pass 2 | 8e-4 | 0, 0.75 | 4, 5, 8 | 210 · 1,024 | 0.70 |

⭐ **Coverage is wide where it matters.** Nodes traversed per branch point — how hard degree-2
suppression works — ran from **3.4** ($\rho$=0.5, low turnover) to **65.2** ($\rho$=8e-4,
$\theta$=0.75), a **19× span**, monotone in turnover as it must be. Structure and topology held
across all of it, and at $n$ up to 3,387 — closing the "validated at $n\le8$, used at $n\le10{,}997$"
gap.

⚠ **Remaining calibration subtlety: per-CHECK multiplicity is not handled.** Each check is flagged at
its own p99 and a pass-2 job runs ~15 checks, so ~14% of correct runs will trip one. Per-*cell*
multiplicity is handled; per-check is not. The 0.96 near-miss (root split, $n$=8, $\theta$=0.75:
2.90 vs 3.03) is exactly what that rate predicts.

## ✅ Cost measured (`03`, 2026-09-20) — forward simulation is feasible everywhere

**⇒ THE BDS DENSITY CAN BE DROPPED.** §S3.3.5's route (c) and the owed branching-time verification
are not needed: nothing in this project exceeds a single Slurm task.

**The atomic unit** — wall-clock for one tree of each arm's largest clone, the thing that cannot be
parallelised away (core-hours can be, by clone and by replicate):

| arm | largest clone | $\rho=0.5$ | $\rho=8\times10^{-4}$ |
|---|---|---|---|
| Initial 157 · M3 373 | | 0.1–0.2 s | 9–48 s |
| M1 1,847 | | 2.5 s | 19.6 min |
| M2 4,443 | | 11.9 s | 1.9 h |
| Subclone 27,224 | | 6.3 min | **70.7 h** (7-day cap) |

Total core-hours at $R$=100: design regime **0.05–25**; Park regime **4.7–192** for four arms and
**16,292** for Subclone alone.

⭐ **The cost model validated OUT OF SAMPLE.** Fitted on $n$ = 4, 32, 210 only, it predicted check
F's independent $n$=1,024 point at $\rho$=8e-4 to **0.91×** on seconds and **1.00×** on peak live
lineages — a 4.9× extrapolation beyond the fitted range. Both halves of the model are now confirmed:
rejection (analytic $P(k{=}n)$ = 32.1 attempts vs predicted 30.3) and population ($n/\rho$).

**Memory is a non-issue**: peak RSS never exceeded **0.18 GB** anywhere. Ask for 4 G, not 16.

⚠ **Caveats on the Park-regime projection.** (i) $\gamma$=1.76 fitted on three points, with Subclone
a 27× further reach even after the $n$=1,024 confirmation. (ii) Clone sizes are **unfiltered** (the
paper's ≥100/≥20-tape cell filter needs the edit tables), inflating cost by 1.3× (M1) to **6.1×**
(Subclone) — so Subclone's 16,292 core-h is likely nearer 2,700. Direction is conservative.
(iii) $R$=100 is still a placeholder.

## ✅ Tree library BUILT (`04`/`05`, 2026-09-21/22)

**1,490/1,490 cells · 29,780 trees · ZERO structure faults · 97 MB** in
`results/tree_library/` (gitignored — reproduce by rerunning `src/`). Check B runs on every stored
tree, so any pruning or slot-assignment failure would surface.

**Design: stock the RANGE of clone sizes, not Park's clone list.** 20 log-spaced sizes 2–11,081
pinning the five observed arm maxima (133 · 210 · 1,607 · 3,387 · 11,081), **plus every integer
2–64**, which is where 97.1% of Park clones sit (median 8) — so nearest-size substitution is exact
there rather than up to 34% off. An arm is reassembled later by drawing at its **true size
distribution**; range is what we stock, distribution is what we draw. ⚠ Not optional: the statistics
this feeds are strongly clone-size dependent (best-$k$ inverts with clone size; the clone floor was
monotone 4/4).

**Per tree:** `parent` int32 + **branch lengths** float32 (not node times — float32 differencing
costs ~9% on the shortest measured branch, $1.4\times10^{-6}T$), seed + accepting attempt index
(replay verified **bit-identical**), a hash, and edit-free summaries (tree length, B1, LTT).
⚠ Prefix-clade-size-by-depth is deliberately absent — prefix clades are defined by shared *edit*
prefixes, so it needs the editing layer.

**⭐⭐ Finding: $\rho$ dominates $\theta$ by 3–5× on coalescent depth.** Median coalescent depth
(fraction of $T$ by which half the lineages exist) spans **+0.27 to +0.34** over $\rho$ 0.25→0.002
against **+0.04 to +0.11** over $\theta$ 0.3→0.7, monotone in 20/20 $\theta$ rows and 5/5 $\rho$
comparisons. ⇒ the unknown we could not get from the paper (mouse $\rho$) matters more than the one
derived from growth arithmetic ($\theta$). ⚑ The $\theta$ effect also shrinks with clone size, so
small clones are the informative ones about turnover.

**In flight (2026-09-22):** `14228126` adds the **$\theta=0$ Yule reference** (8.8 core-h);
`14228127` extends $\rho$ to **{0.05, 0.02, 0.0005}** (212 core-h), because SciPhy found a nominal
sequenced/population ratio over-stated effective $\rho$ by **4–16×** — so treat every nominal
$\rho$, including Park Initial's 0.24, as an upper bound.

## ⚠⚠ $\rho$ per arm — two established, three not (2026-09-21)

| arm | $\rho$ | basis |
|---|---|---|
| Initial / Pre-TX | **~0.24** (upper bound) | 37,810 edit-table cells from a stated ~160,000 pool |
| Subclone | **0.1–0.3** (bounded) | 8 colonies from single founders, 295–11,081 cells each, vs ~35,000 expected after 35 d |
| M1 · M2 · M3 | **unknown** | growth tracked by IVIS bioluminescence = relative flux, no cell count or tumour mass reported |

⚠⚠ **$\rho=8\times10^{-4}$, used throughout this analysis before 2026-09-21, is SciPhy's fixed value
for HEK293T** (`sciphy_notes.md` lines 899, 1742) — **not a Park measurement.** It is accidentally
the right order for the mice and 2–3 orders off for Initial and Subclone.

⚑ **$\theta$ from the only edit-independent handle:** the pool went ~8,000 → ~160,000 in 10 days
(paper p.4) = a realised net doubling of **55.5 h**, against NCI-H1299's intrinsic 22–30 h, giving
$\theta = 1-r/b =$ **0.46–0.60** — a *lower* bound, since the literature figure is itself net.
⚠ Confounded: the deficit could be slower division rather than more death. Hence a grid, not a point.

## Scope change, recorded (2026-09-20)

**Forward simulation is now the primary route, not the validator.** §S3.3.5 listed direct
conditioned sampling from the BDS density as "the one to build", with forward+reject as validation
only. That ordering assumed the tree generator serves adjudication against Park. On the stated
priority — **purposes I and IV, with IV leading** — it is the other way round: the design sweep has
$n$ and $\rho$ as *design variables*, not quantities to be matched, so it needs no conditioning, no
rejection and no density. The $n/\rho$ blow-up that makes forward simulation expensive is entirely
an artefact of Park's $\rho \approx 8\times10^{-4}$; a design study at $\rho \sim 0.1$–$1$ is cheap.

⇒ build forward first, **measure** what it costs, and reach for the density only where the
measurement says it is needed — which is the sequence `03_cost_curve.py` exists to decide.

## The three scripts

| file | what it does |
|---|---|
| `src/01_bdtree.py` | the simulator. Gillespie loop over the population size → genealogy → $\rho$-sampling → reconstruction. Two implementations: `simulate_tree_naive` (one event at a time, nothing pruned until the end — the oracle) and `simulate_tree` (production) |
| `src/02_validate_tree.py` | five correctness checks against closed forms already verified in §S3.3, reported in **se** not p-values: count law, structure, root split, labelled histories, naive-vs-fast |
| `src/03_cost_curve.py` | the cost measurement and the per-arm projection. Decides whether the BDS density is needed at all |

## Two design points worth keeping

**⭐ The two-pass split, and why it is exact.** At each event the lineage it lands on is chosen
uniformly and independently of everything else, so the law factorises as (count trajectory) ×
(uniform lineage assignments | count trajectory), and the acceptance test reads only $N(T)$ and a
$\mathrm{Binomial}(N(T),\rho)$ draw — never the genealogy. So pass 1 can get $N(T)$ with **no
genealogy at all**, vectorised in numpy, and the $\approx 2.7n$ rejected attempts never allocate a
node; pass 2 replays the identical trajectory from the same seed and builds the genealogy only on
acceptance. Given $\mathrm{Binomial}(N,\rho)=k$ the captured set is a uniform $k$-subset, which is
how pass 2 draws it.

**⭐ The cost study has exactly two free numbers.** $\alpha$ and $\beta$ depend on $(b,\delta,T)$
only through $u=e^{(b-\delta)T}$ and the turnover $\theta=\delta/b$, so at fixed $\rho$ the whole
rejection cost is a function of $(n,\theta)$ and $T$ is a pure scale choice, fixed at 1
(`rate_params`). ⚠ That is **not** true for adjudication against Park, where $T$ is pinned by
$\Lambda=\lambda T$ at the measured $\hat\Lambda$ — it is true for the cost study only.

## ⚠ Corrections to my own proposal of 2026-09-20, before anything ran

1. **Online extinct-pruning was proposed as the memory mitigation; it is not needed.** Working the
   arithmetic: total nodes $\approx 2N(T)/(1-\theta)$, so online pruning saves only ~2× at
   $\theta=0.3$ and ~4× at $\theta=0.75$ — the same order, not the difference between feasible and
   OOM. At Subclone scale the simple "keep everything, prune at the end" path is a few GB.
   ⇒ built the simple version, and `03` records peak RSS so the question is settled by measurement.
2. **Memory is not the binding constraint; wall-clock is.** The genealogy loop (the only sequential
   part) is paid once per acceptance and is a fraction of a percent of the total, so it is
   deliberately left unoptimised.

## What is owed

1. Run `02`. It decides whether the simulator is correct; `03`'s numbers mean nothing until it passes.
2. Run the `03` pilot at small $n$, read $\gamma$, $\kappa$, $\tau$ and `MaxRSS` off a **completed**
   job, then size the large-$n$ points from measurement (CLAUDE.md, "Right-size the request").
3. $R$ — replicates per clone — is deliberately unset. It multiplies the whole cost linearly and is
   set by what the rung 0/1/2 comparison needs. Justin's call, deferred until the timing is known.
4. The forward-vs-BDS equivalence check, **only if the cost measurement says the density is needed**.
   `01_bdtree.py` already carries $h$ and $G$ (fenced, unused by the simulator) and `02` persists
   branching times, so that check needs no re-simulation.
