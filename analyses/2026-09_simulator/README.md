# simulator — findings

**State 2026-09-23: the tree layer is BUILT, VALIDATED, COSTED, and the library is COMPLETE
(1,870 cells, 37,380 trees, zero structure faults). Figures S1 (method) and S3 (pull of the
present) are built, with a LaTeX write-up in `writeup/`.** Editing and dropout are separate later
layers, not yet started.

## ⭐ Session 16 (2026-09-22/23) — figures, write-up, and what the trees record

**Write-up:** `writeup/simulator_figures.tex` (pdflatex; `% !TEX program` magic comment, and a
pdflatex×2 recipe in the gitignored `.vscode/settings.json` because iris has no `latexmk`).
One section per figure, with every defining equation of the tree layer. Justin edits it directly.
Figure PDFs and the compiled document are **not committed** (the repo ignores `*.pdf`); the
PNGs, `.tex` and scripts are. Rebuild the figures with `src/10` and `src/11`.

**Fig S1 (`src/10_fig_method.py`, three standalone panels).** (a) One tree through three stages:
full process → Binomial capture → reconstructed tree, with the spliced divisions marked (21 of
them at $n$=8, $\rho$=0.25, $\theta$=0.5; 77 genealogy nodes → 15). (b) Conditioning by
rejection: 60 attempts, 4 accepted, 30 extinct; beside it the law of $K$. (c) Library inventory.
⭐ **New closed form, verified:** given $K\ge1$, the capture count is geometric,
$P(K=k)=s(1-\beta')\beta'^{k-1}$ with $\beta'=\rho\beta/c$, $c=1-\beta(1-\rho)$,
$s=P(K\ge1)=(1-\alpha)\rho/c$. Against 40,000 attempts: $P(K=8)$ 0.0190 exact vs 0.0183 simulated,
$P(K=0)$ 0.515 vs 0.518, worst of 59 cells 2.6 se. Reproduces the log's 32.1 attempts per accepted
tree (32.0). ⚠ The first draft of `10` (never committed) drew an unconnected panel 3 and
captured exactly $n$ of the survivors instead of Binomial($N,\rho$); both fixed.

**Fig S3 (`src/11_fig_pull_present.py`) — an ILLUSTRATION, not a test** (Justin's call; a
closed-form-vs-library test and a matched-pair confounding test were proposed and **declined** in
favour of this). Mean $\ln L(t)$ over 20 trees at $n$=63, from stored branch lengths on a
400-point grid; $\rho$=1 trees simulated fresh (not in the library). Tip rate = $d\ln L/dt$ over
the last 10% of $T$, per lineage per unit $T$; baseline $\theta$=0, $\rho$=1 is **4.0**.
- **(a) the pull of the present exists:** at $\rho$=1 the tip rate is 4.0 / **7.5** / **10.2** at
  $\theta$ = 0 / 0.5 / 0.7 — 1.9× and 2.6× the baseline.
- **(b) capture reverses it:** at $\theta$=0.5, 7.5 → 3.2 → 0.49 → **0.008** at $\rho$ = 1 / 0.25
  / 0.02 / 0.002. At mouse-like $\rho$ nothing branches in the last tenth of $T$.
- **(c) ⚠ MY PREDICTION WAS WRONG.** I predicted the $\theta$ lines would bunch together at low
  $\rho$; they do not. Turnover changes the curve's **shape** at $\rho$=1 and **shifts it later**
  (~0.05–0.1 $T$) at $\rho$=0.002. Raising $\rho$ also shifts it later, so at low capture higher
  turnover and higher capture move the tree the same way — a qualitative, untested reading.
- ⚠ **What is held fixed:** same $n$ and $T$, so lowering $\rho$ also raises the rates (the clone
  must grow $1/\rho$ larger): $b$ = 8.3 at $\rho$=1 vs 20.7 at $\rho$=0.002 ($\theta$=0.5). Panel
  (b) compares the same observed clone under different capture, not capture at fixed rates.

**⭐ Why it matters — the framing agreed with Justin: turnover inference is NOT a central aim;
the point is which timepoints are recorded well.** $L(t)$ is the number of **independent records
of time $t$** — captured cells descending from one lineage share its tape — so the LTT curve is the
recorder's sampling depth over time, and $n/L(t)$ the redundancy. Four consequences, the first two
about precision and the last two about bias:
1. Precision of any $\xi_i(t)$ estimate follows $L(t)$, not $n$. Early windows are always thin
   within a clone; at low $\rho$, late windows have lineages but almost no branch points, so late
   edits sit on long terminal branches and are poorly timed. ⚑ Sparse capture *decorrelates* cells:
   $L(0.5)\approx40$ at $\rho$=0.002 vs $\approx3$–4 at $\rho$=1, $\theta$=0.7 (read off the figure).
2. Pooling cells without the tree weights lineages by clade size — random and parameter-dependent.
   The pruning likelihood removes this if the tree is right.
3. **Survivorship bias (genuine bias):** deep time is seen only through lineages that survived and
   were captured. If the signal affects survival or division (metastatic fitness in Park), the
   recorded history over-represents the winners. **Not visible in the current simulator** (rates
   are signal-independent); needs signal-dependent $b$/$\delta$.
4. Hypothesis, needs the editing layer: tape saturation (≈4.5–5 of 6 sites) is the recorder's own
   bias towards early events, so the best-recorded window sits between the tree's thin early
   stretch and the tape's full late one — a design result for purpose IV.
⚠ Assumes clones start recording together; if so, early windows are replicated across clones
(149–2,946 per arm). Not checked against the Park protocol.

**⚠ Corrections to the record found this session (none moves a scientific number):**
1. **72,000 distinct validation trees, not 84,000** — pass 1 and pass 2 reran $\rho$=0.5,
   $\theta$=0.3 with the same seed (identical statistics), so 12,000 were counted twice.
2. **Check A (count law) never ran at the Park regime.** It is hard-coded to $n$=64, $\rho$=0.5
   ($E[N(T)]$=128); all four result files hold the identical check-A numbers. B–F at
   $\rho=8\times10^{-4}$ exercise the batched count path but test shape, not the count law, at
   $N(T)\sim10^6$. "Validated at both capture fractions" holds for B–F only.
3. **The attempts envelope in `04`/`03` is only asymptotic.** $2.7n/(1-\alpha)$ vs exact
   $1/P(K=n)\approx n e^{s}/s^2$, ratio $\approx e^{s-1}/s$ with $s\approx1-\theta$: 1.06× at
   $\theta$=0.3, **1.65× at 0.7** (1.24× at the S1 settings, 42.5 vs 52.7). The ratio across
   turnovers (1.56×) matches defect 8's measured actual/predicted trend (0.79→1.27, 1.61×), which
   the log attributed to the Python genealogy loop. **Plausibly most of defect 8; untested against
   per-task timings.** Affects job sizing only.
4. ⚠ In discussion I claimed the textbook upturn is never visible in our library ($\rho<1-\theta$
   everywhere). That criterion is for the unconditioned per-lineage rate; for trees conditioned
   on $n$ the branching-time density piles up at the present when $(1-\theta)/\rho\le2$, which
   includes $\rho$=0.25, $\theta\ge0.5$. Fig S3's $\rho$=0.25 panel is consistent with that.
5. `results/validate_rho8e4.json` still says FAIL — the pass-1 verdict under the old flat 3-se rule,
   superseded by `validate_pass2_rho8e4.json` (PASS, 0.70).
6. The recorded pruning-load range "3.4–65.2" is **2.9–65.2** in the pass-2 JSON.

⚑ **For the PI / Jihye Park:** the mouse capture fraction is the one number that decides how the
mouse trees' timing is read (Fig S3c). Ask: total viable cells per dissociation, or tumour mass,
per organ at day 45.

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

**✅ Both COMPLETE (2026-09-23), zero structure faults:** `14228126` added the **$\theta=0$ Yule
reference** (24/24 tasks); `14228127` extended $\rho$ to **{0.05, 0.02, 0.0005}** (140/140; its
tail tasks ran 1.4–2.6× over prediction — see correction 3 above). Motivated because SciPhy found a nominal
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
