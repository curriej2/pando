# Analysis: simulator — birth–death tree → sequential editing → dropout

**Question.** Build the forward simulator that supplies (i) the **homoplasy null** — the oldest
outstanding debt in the project — and (ii) the **design sweep** over $(p,\lambda,m,k,j,\ell)$.

**Opened 2026-09-15**, when the figure programme for heritable silencing was closed (see
`analyses/2026-08_park-compatibility`). Theory record: `notes/sciphy_notes.md` **§S3**
(SciPhy's own simulation design, the birth–death sampling model derived in full, both verifications).
Running record: `notes/analysis_log.md`, session 13.

**⭐⭐ STATE (2026-09-23): TREE LAYER BUILT, VALIDATED, COSTED; LIBRARY COMPLETE (1,870 cells,
37,380 trees, zero faults). Figures S1 (method) and S3 (pull of the present) built; LaTeX write-up
in `writeup/simulator_figures.tex`, which Justin edits directly. Editing and dropout layers NOT
started.** README is the findings record; this file's design sections below are the 2026-09-17
design and are partly superseded (the BDS-density route was dropped 2026-09-20 in favour of forward
simulation + rejection — see README).

⭐⭐ **2026-09-24: switch detectability MEASURED in closed form (`src/12`, Fig S4, README "Session 17
(cont.)").** At $k=30$ a lineage-specific switch is callable when it is **long ($\ge$1–2 edits/tape),
strong ($\ge$15×) and early (first ~60% of the experiment)**, in clades of only **2–5 cells**; short
($\le$0.5 edits/tape), weak (2×) and late (last ~20–30%) switches are out of reach. Branch-level
timing costs little ⇒ **the tree side is not the bottleneck; build the editing layer, then a
fitness-neutral state layer.** The information-budget pair it replaced ("PROPOSAL (2026-09-24)"
below) was never run — four defects found in review. Theory in `notes/sciphy_notes.md` **§S4**
(§S4.8 records the answer).

⚑ **Framing agreed 2026-09-23:** turnover inference is NOT a central aim. Fig S3 matters because
$L(t)$ = independent records of time $t$, so tree shape decides which timepoints of the signalling
history are recorded well (precision), and survivorship decides whose history is recorded (bias).
See README "Session 16".

---

## The four purposes this simulator serves (agreed with Justin, 2026-09-14)

Justin chose **I and IV** to lead.

| # | purpose | what it buys |
|---|---|---|
| **I** | **Adjudication** — quantify the *consequence* of things already established | converts "dropout is heritable and structured" into "and it costs you X"; decides whether the A9 absorbing state earns its place in the pruning |
| II | Validation — our own statistics against known ground truth | every statistic so far is tested only against permutation nulls, which probe $H_0$ and say nothing about the alternative |
| III | Method development — does the extended likelihood work | §F.3's architecture-licensing experiment; ML tree search (Program C); SBC |
| **IV** | **Design** — what recorder should be built next | the thesis question (§I.7.7); publishable without new data |
| V | Generalisation — the 11-tape regime | `2026-09_dtt-mouse` without its unreleased tape matrix |

⚠ **Standing warning (§H.6.12): a simulator is also a model.** Omit A5/A6 and it is exactly as wrong
as the closed-form bound, only more confidently so, because it ships error bars describing Monte
Carlo noise rather than model error.

## ⚠⚠ Scope corrections already paid for — do not reintroduce

1. **The compatibility statistic is NOT the deliverable.** Max-compatible-set died in the 2026-09-01
   pivot. The pivot preserved compatibility as a *diagnostic* ("the diagnostics are the deliverable"),
   but that is now stale too — the dropout diagnostic has been superseded by the ladder, $Q$, the
   variogram and unanimity, on five arms, declared sufficient 2026-09-14. ⇒ **An Initial rerun and a
   Subclone pricing were proposed and WITHDRAWN.** Compatibility survives only as a free end-stage
   check: Mouse3 is complete on disk, so if the finished simulator reproduces its 63.56% / 94.66% /
   +31.09-point spread through `../2026-08_park-compatibility/src/14`, that costs nothing.
2. **The homoplasy null is narrower than first framed.** For the mouse arms it is largely
   unnecessary — script 06 already *measured* level-0 recurrences at 0.03–0.11 per prefix node with
   0.4–1.1% of nodes carrying any. It is load-bearing only for **Subclone** (median clade 933 cells,
   16.0 mean recurrences per node, 40.8% of nodes with ≥1), where the strongest dropout signals also
   live so the two processes are genuinely entangled, and for the **design sweep**.
3. **The currency should change.** Incompatibility is a perfect-phylogeny quantity; the likelihood
   route cares about *reconstruction accuracy*. The question worth answering is "at what clade size
   does homoplasy start to cost accuracy" — the $m^{*}\approx11$ birthday threshold in RF/PI error
   against known truth. ⇒ **fig 5 to be renamed and rescoped** away from its skeleton-era title.

## The design, as settled

**Three rungs** — the ladder *is* the figure, mirroring `fig6a` which already works on a PI:

| rung | simulated | what the step measures |
|---|---|---|
| 0 | tree + editing, complete observation | the pure homoplasy floor |
| 1 | + dropout with **no** lineage structure (per-tape $\beta_z$ × per-cell $\alpha_c$, independent) | what ordinary technical dropout manufactures |
| 2 | + heritable Dollo loss (row A9) | what the silencing we spent five sessions establishing actually costs |
| — | observed − rung 2 | genuine model violation, error, doublets |

Rungs 0–1 are the first deliverable; rung 2 needs a generative loss model we have not fitted (we
measured the phenomenon, not a rate).

⚠ Both observed conventions are needed because they demand different nulls: **missing-as-absent**
(63.56% on Mouse3) scores every missing entry as "does not carry this prefix", so dropout
*manufactures* incompatibility and a dropout-free null is meaningless there; **missing-excluded**
(94.66%) conditions on $D_1\cap D_2$, which looks dropout-free but is not, since lineage-structured
missingness makes the conditioning non-random.

## Language: Python. Settled, with reasons

SciPhy's simulator (`src/sciphy/evolution/simulation/SimulatedSciPhyAlignment.java`) is **260 lines**
and takes the tree as an **input** — it only decorates a given tree; trees come from BEAST's
birth–death simulators or R's TreeSim. Forward simulation touches each edit once (Park scale
≈ 32 M draws per replicate — seconds in vectorised numpy). The paper's optimisations (subtree
likelihood caching 3.5×, tape-level threading, 8× overall at 1,000 sequences) are all **MCMC
likelihood evaluation**, a different workload.

⇒ simulator in Python; if SciPhy *inference* is wanted on simulated data, drive BEAST as an external
process on generated XML (templates in `examples/`). Java is forced only if the A9 absorbing state
has to live inside BEAST's MCMC. Repo clone lives in the session scratchpad, not committed.

## The tree model — settled, both parts verified

**⭐ Shape and time factorise.** For a constant-rate birth–death, conditional on tip count the
topology is independent of $b$ and $\delta$. Verified across turnover $\delta/b$ = 0, 0.14, 0.75 at
$n=8$: root-split distributions agree with Yule–Harding ($P=1/(n-1)$ per ordered split) to within
2.0 se on the worst of twelve cells. ⇒ **shape is parameter-free; all of $b,\delta,\rho$ act on the
branching times.**
⚠⚠ "Random binary tree" is ambiguous — uniform over *topologies* (PDA) ≠ uniform over *labelled
histories* (Yule). PDA is far more imbalanced, most libraries mean PDA, and clade sizes set $m$.

**⭐⭐ Clone sizes do NOT come from one birth–death.** A homogeneous BD gives a geometric
(exponentially bounded) tail. Fitting a geometric to each arm's median clone size and asking how many
clones it predicts at least as large as the largest observed — null expectation if the model were
right is ~1, observed is 1 in every arm:

| arm | clones | median (cells) | max (cells) | expected # ≥ max |
|---|---|---|---|---|
| Mouse3 | 149 | 4 | 210 | $2.8\times10^{-14}$ |
| Mouse1 | 295 | 4 | 1,607 | $4.0\times10^{-119}$ |
| Mouse2 | 216 | 2 | 3,387 | underflow |
| Initial | 2,946 | 6 | 127 | $1.4\times10^{-3}$ |
| Subclone | 15 | 997 | 10,997 | $7.2\times10^{-3}$ |

Mouse2's largest clone holds 62.9% of that arm's cells. ⇒ **DESIGN DECISION: do not simulate clone
sizes. Take them from the data; simulate the tree within each clone conditioned on its observed
size.** (What SciPhy's own benchmark does for comparator trees.) This absorbs across-lineage rate
heterogeneity, leaving only constancy *within* a clone — testable against prefix-clade-size-by-depth.

**⚠ Identifiability: fit two effective parameters, not three.** $(b,\delta,\rho)$ collapse — verified
that BD$(1.0,0.5)$ at $\rho=0.7$ and BD$(0.7,0.2)$ at $\rho=1$ give matching conditional clade-size
distributions, though different $P(0$ sampled$)$ (0.498 vs 0.283); the equivalence is about the tree
*given it exists*. The rescaling needs $\delta>b(1-\rho)$, which fails at Park-like $\rho$ (equivalent
$\delta=-0.499$ at $\rho=0.0008$), so it is a statement about the density, not a simulable process.

⇒ With shape free and $\rho$ collapsed, **the tree contributes essentially one effective number**:
how deep the coalescences sit. Constrained by prefix-clade sizes by depth (1,567,321 nodes, on disk).

## How to generate a tree — three routes, one choice

| route | verdict |
|---|---|
| forward + reject to $n$ tips | correct by construction, but **0.7% acceptance at $\delta/b=0.75$** (measured) and hopeless at $n=10{,}997$. **Validation only** |
| forward, stop at $n$ lineages | conditions on $n$ but **not** on $T$ — and $T$ is pinned, because $\Lambda=\lambda T$ is the editing depth we measured. **Unusable** |
| ⭐ **direct conditioned sampling** | (a) shape by uniform random joining; (b) $n-1$ branching times from the BDS density. Exact, $O(n)$, scales. **Build this** |

### The BDS branching-time density — derived and largely verified

Both building blocks fall out of the pgf (§S3.3.2) with no new machinery. A lineage at time $t$
before the present leaves no sampled descendant iff all its extant descendants are missed, so
$p_0(t)=\mathbb{E}[(1-\rho)^{N(t)}]=F(1-\rho,t)$. Writing $h=1-p_0$ and $K=b(1-\rho)-\delta$:

$$h(t)=\frac{\rho(b-\delta)}{b\rho+Ke^{-(b-\delta)t}},\qquad
p_1(t)=\frac{\rho r^2e^{-rt}}{[b\rho+Ke^{-rt}]^2},\qquad h'=K\,p_1$$

Because $p_1$ is (up to $K$) the derivative of $h$, the branching-time CDF is closed-form:

$$G(x)=\frac{h(x)-h(0)}{h(t_{or})-h(0)},\qquad h(0)=\rho$$

⇒ **inverse-CDF samplable analytically** ($h$ is Möbius in $e^{-rx}$): no rejection, no root-finding,
$O(1)$ per branching time.

**Verified so far** (150,000 forward replicates, $b=1.0,\delta=0.3,\rho=0.5,T=3.0$, $n=5$):
$p_0(0)=1-\rho$ and $h(0)=\rho$ exactly; $p_1=h'/K$ true by finite differences; exactly 4.00
branching times per tree (confirms the reconstructed-tree extraction); empirical-vs-$G$ quantile
deviations −0.0026 to +0.0077 across the 10th–90th percentiles; $\sup|{\rm emp}-G|=0.00800$ on
31,468 times from 7,867 accepted trees.

⚠ **NOT yet settled, and this is the next task.** The naive KS critical value at $n=31{,}468$ is
$1.36/\sqrt n=0.0077$, so 0.00800 sits fractionally above it — but the four times from one tree are
dependent, so that $n$ is not the effective sample size; against 7,867 *trees* the critical value is
0.0153 and 0.008 is comfortably inside. The five deviations also run $-,-,+,+,+$, a mild sign pattern
that is either noise across correlated quantiles or a small systematic bias.

## PROPOSAL (2026-09-24b) — how far can a lineage-specific switch be detected? ✅ APPROVED, RUN

Approved by Justin 2026-09-24 with three additions (a $p_1$ scan, $n=8$, $\rho$ up to 0.8). Ran as job
14592174. Full design, verification and results: README "Session 17 (cont.)". In one paragraph: a
switch turns on at a point on the tree, is inherited, and lasts $\Lambda_e$ edits/tape; expected
evidence is $k(1-d)\sum_bW_b\,\mathrm{KL}(\bar p_b\|p_0)$ nats with edits placed only to their branch
(or timing known, as the upper bound); callable ⇔ $\ge3+\ln H$. Scanned over $n$ (6) × $\rho$ (5) ×
$\Lambda_e$ (4) × $F$ (4) × $p_1$ (5) on 20 trees per cell. ⚠ The agreed headline (fraction of all
switch locations) turned out to be dominated by late switches; the by-start-time reading is recorded
alongside it, flagged as found after the run.

## PROPOSAL (2026-09-24) — the information budget, in closed form, before the editing layer

**⚠⚠ SUPERSEDED 2026-09-24, NEVER RUN — replaced by (2026-09-24b) above.** Kept as the record. Four
defects found in review (README "Session 17 (cont.)"): Table 2 never stated $p_0,p_1$ and used
$p_1=0.20$ in two columns but 0.24 in the third, which moved the $k=30$, 24× call across the 3-nat
line (3.36 vs 2.80); (b) scored a branch in isolation, which is not a bound; (a) was called both a
lower and an upper bound; and (a)'s row-normalised kernel makes $\sigma_i\propto1/\sqrt G$, so $R$
moves with $G$ by construction. `src/12` and `src/13` are now the switch-evidence scripts, not these.
Theory behind it: `notes/sciphy_notes.md` §S4 (derived 2026-09-23/24). Placed *ahead* of the editing
layer at Justin's request, because it decides that layer's specification for a fraction of the cost.

### 1. The question

Can the regime Justin expects to work in — **$k\approx30$ tapes, $N=6$ sites** — carry *per-lineage*
signal history at all? If yes, the editing layer must represent and be validated against
lineage-resolved $\xi(t,\sigma)$. If no, the modelling target narrows to a **population-level**
$\xi(t)$ plus high-contrast on/off calls, and the design sweep's headline becomes "what $k$ and
$\lambda$ would be needed" — a recommendation to the experimentalist rather than an inference method.
Either answer changes what we build, so it is worth answering first. ⚑ A negative answer is itself
the §I.7.7 thesis question answered with numbers attached, i.e. purpose IV, not a failure.

### 2. The math — every symbol defined

Two independent closed forms. **Neither uses the tree.**

**(a) Temporal rank $R$ — how many independent temporal features a tape can carry.**
Work in $\Lambda$-time $\tau=\Lambda(t)/\Lambda(T)\in[0,1]$, the fraction of the experiment's total
integrated editing rate spent by time $t$ ($\tau=t/T$ if the rate is constant); $\Lambda$ is in edits
per tape, $\tau$ dimensionless. Let $g(\tau)\in[0,1]$ be the **signal symbol's share of insertions**
at $\tau$ — the quantity ENGRAM modulates; no-signal value is $g$ constant. Let
$D=\min(N(T),N)$ be a tape's observed depth (sites filled, 0–6). The $\Lambda$-time of that tape's
$j$-th edit has a closed-form density $K_j(\tau)$:

- $D=d<N$ (unsaturated): $K_j=\mathrm{Beta}(j,\,d{+}1{-}j)$, from the order-statistics result (§S4.2);
- $D=N$ (saturated): $K_j(\tau)\propto f_{\Gamma(j,1)}(\Lambda\tau)\,F_{\Gamma(N-j,1)}\big(\Lambda(1-\tau)\big)$,
  the $j$-th jump of a unit-rate process conditioned on the $N$-th landing inside the experiment.

The mean share observed at site $j$ is then a **linear smoothing of the whole history**,
$\bar g_j=\int_0^1 K_j(\tau)\,g(\tau)\,d\tau$. Discretise $\tau$ on $G$ grid points (rows normalised
to sum to 1) to get $\bar g=Kg$ with $K$ of size $N\times G$, and take its singular values
$\sigma_1\ge\cdots\ge\sigma_N$ (dimensionless). With $n_B$ edits observed at a site, measurement
noise is $\sigma_\varepsilon=\sqrt{\bar g(1-\bar g)/n_B}$, and mode $i$ is **resolvable** iff
$\sigma_i>\sigma_\varepsilon/A$, where $A$ is the signal amplitude (peak-to-trough swing in $g$,
dimensionless). $R$ = the number of modes clearing that bar. $R=1$ means only the experiment-average
signal is recoverable — no history at all.
*Assumptions:* site index observed without error; edits at a site treated as independent draws; only
the $N$ **linear** (marginal) functionals used. Bigrams are *quadratic* in $g$ and are excluded, so
$R$ is a **lower bound** on what the full likelihood could extract.

**(b) Per-branch evidence $\mathcal I$ — whether one lineage's state is callable.**
For one branch: $k$ = tapes per cell (count); $\Delta$ = the branch's integrated editing rate (edits
per tape); $d$ = per-tape dropout probability; $U$ = fraction of tapes not yet saturated. Expected
edits recorded on that branch is $m=k(1-d)\Delta U$ (edits). Signal symbols number
$\mathrm{Poisson}(mp)$ with share $p=p_0$ (off) or $p_1$ (on). The expected log-likelihood ratio in
favour of the truth — the Kullback–Leibler divergence between the two Poissons — is

$$\mathcal I_{\rm on}=m\big[p_1\ln(p_1/p_0)-(p_1-p_0)\big],\qquad
  \mathcal I_{\rm off}=m\big[p_0\ln(p_0/p_1)-(p_0-p_1)\big]\quad\text{nats}$$

**No-effect value 0 nats** (at $p_1=p_0$). One nat = one $e$-fold of odds. Both directions are
reported because they are asymmetric: $p_0$ small means few symbols, so "off" is the weaker call.

### 3. What could make it wrong

**The confound both defeat:** neither uses relatedness, so neither can be contaminated by the
project's recurring circularity — relatedness read from the same edits whose readability dropout
controls. Every input is a design parameter or a directly measured marginal, not an estimate.

**The confound neither defeats, and it is the important one:** both are **upper bounds on achievable
resolution**. They credit the site index with locating edits in $\Lambda$-time and never ask whether
$\Lambda$-time can be converted to *calendar* time — which needs the tree, and which Fig S3(b) already
shows fails over the last tenth of $T$ at the mouse capture fraction (§S4.7). They also treat the $N$
site estimates as independent when they come from the same tapes, which is optimistic again.
⇒ **Only a negative result is decisive. A positive one licenses building, not claiming.**

**Abandon the approach if:** $R$ moves materially with the grid size $G$ or with the depth mixture.
That would mean the rank is a discretisation artefact rather than a property of the recorder, and the
right response is to go to simulation, not to reinterpret the number.

### 4. Output — the mock, and the reading rule

**Table 1 (`src/12_temporal_rank.py` → `results/temporal_rank.json`).** Six rows, one per singular
mode, at $N=6$ and a stated $\Lambda$. **Placeholder numbers — the spectrum cannot be done by hand.**

> *Reading rule:* mode $i$ is resolvable when your per-site edit count $n_B$ exceeds column 3. The
> **temporal rank $R$** is the number of rows your experiment meets. $R=1$ = only the time-average
> of the signal; $R=2$ = average plus a monotone trend; $R\ge3$ = a genuine shape. **Decisive:**
> $R\le2$ at $k=30$ kills per-lineage *graded* history. **Marginal:** $R=3$.

```
  Lambda = 5.0,  N = 6,  tau-grid G = 200,  depth mixture from Poisson(Lambda) truncated at N
  i    sigma_i / sigma_1     n_B needed to resolve mode i at amplitude A = 0.10
  1          1.000                        9
  2          0.31                        94
  3          0.08                      1400
  4          0.02                     23000
  5          0.004                   580000
  6          0.0009                 1.1e+07
```

**Table 2 (`src/13_branch_evidence.py` → `results/branch_evidence.json`).** Four rows, indexed by
tapes per cell $k$ (11 = mouse-embryo DTT, 30 = Justin's expected regime, 100, 166 = Park).

> *Reading rule:* no-effect value 0 nats. **$\ge3$ nats ($\approx$ 20:1) = a branch call worth
> making; 1–3 nats = usable only pooled over several branches; $<1$ nat = not callable.** Decisive
> question: does the row for *your* $k$ clear 3 nats at the fold-change your recorder delivers?
> (ENGRAM measures 15–24× swings: 23.8× for 3' FT ±TNF, 15.1× for a high/low CRE pair.)

```
  Delta = 0.5 edits/tape/branch,  d = 0.44 (Park per-tape dropout),  U = 0.75
                            2-fold          10-fold         24-fold
    k       m = k(1-d)DU   on    off      on    off      on    off
   11           2.3       0.09  0.07     0.65  0.31     1.23  0.46
   30           6.3       0.24  0.19     1.77  0.84     3.36  1.25
  100          21.0       0.81  0.64     5.89  2.81    11.2   4.16
  166          34.9       1.35  1.07     9.79  4.68    18.6   6.91
```
⚠ These Table-2 entries are the formula evaluated by hand at the stated assumptions, shown so the
reading rule can be judged — **not results.** The run replaces the assumed $\Delta$ with the branch
-length distribution measured from the tree library, $U$ with the saturation profile at the arm's
$\hat\Lambda$, and $d$ with each arm's measured per-tape dropout, and sweeps all three.
⚑ Note what the hand arithmetic already hints: **dropout and saturation roughly halve the effective
tape count**, dropping the 10-fold call at $k=30$ from 4.2 nats (raw $k\Delta$) to 1.8 — callable to
marginal. That is the swing the run has to pin down.

### Cost, and where it runs

Numpy SVD on a $6\times200$ matrix and closed-form arithmetic: **seconds, one core**. Per the
standing rule ("run every analysis as a Slurm batch job"), submit anyway via
`scripts/submit.sh src/12_temporal_rank.py --mem 4G --time 00:10:00 --cpus 1`, and likewise `13`.
Right-sized: 4 G, 10 min, partition `cpu`.

---

## ⇒ NEXT, in order (updated 2026-09-24; the 2026-09-23 list is superseded)

0. ✅ **Switch detectability** (PROPOSAL 2026-09-24b) — DONE. ⇒ editing layer first, then a
   fitness-neutral state layer validated on long, strong, early switches; tree-side timing
   refinement deprioritised (item 3 is about rate inference and stands on its own).
1. **The editing layer** — decorate library trees with sequential edits at measured per-tape
   rates, site 6 handled separately, depth 0 censored; editing and dropout must be fitted
   **jointly** (park-compatibility CLAUDE.md). ⭐ Build **route (ii)**, the partial-sum sampler of
   §S4.3 — bounded at $N$ draws per tape-branch, sorted by construction, no waste — and keep route
   (i) as the cross-check. Carry **edit times** and $\Lambda$-time from the start (§S4.1): a
   time-varying *composition* is then free, and a time-varying *rate* is one extra knot.
2. Then, with edits: tape saturation vs $L(t)$ — the "best-recorded window" hypothesis (README
   Session 16, point 4); a proposed figure is $L(t)/n$ with the tape-saturation profile overlaid.
   ⭐ §S4.7 sharpens it: $L(t)$ sets precision, $L'(t)$ sets temporal resolution.
3. **$\lambda(t)$ vs $(\theta,\rho)$ confounding** — §S4.7's named confound, testable on the existing
   library by putting a known $\lambda_0(t)$ on trees across the $\theta$/$\rho$ grid.
3b. ⭐ **The state layer (fitness-neutral)** — paint a two-state CTMC path onto library trees
   (§S4.5, two-pass), then check with the §S4.4 likelihood whether the switches Fig S4 calls
   recoverable survive tree error and homoplasy. Target regime: $\Lambda_e\ge1$–2, $F\ge15$, start
   in the first ~60%, clades $\ge$2–5 cells.
4. Survivorship bias needs signal-dependent $b$/$\delta$. ⚠⚠ §S4.5: **painting states onto library
   trees is exact only if the state is fitness-neutral.** Proliferative states require a multi-type
   birth–death in the *tree* layer, and the library becomes a stepping stone.
5. Open, lower priority: check A at a Park-scale point; refit the attempts envelope with the exact
   $1/P(K=n)$ (README correction 3); S2 validation figures were proposed and Justin declined them.
6. ⚠ Owed to Jihye Park: total viable cells or tumour mass per organ at day 45 (mouse $\rho$).

### Decisions taken with Justin, 2026-09-23/24

- ✅ **Time-varying composition $\xi(t)$ in from the start** — even though it will not be varied yet.
- ✅ **Time-varying rate $\lambda(t)$ too** — Justin's experiments will span states of differing
  proliferation, and inferring time-varying rates from trees is a stated goal.
- ✅ **Route (ii)** (partial sums) as the production sampler. ⚠ This corrects my own earlier
  recommendation of route (i) as "the fast one" — in partial-sum form route (ii) vectorises equally
  well and has a *bounded* draw count, so it dominates.
- 🕐 Still open: piecewise-constant vs spline $\lambda_0(t)$; whether the first build includes
  $\lambda_\sigma$; whether edit times are persisted or replayed from a seed; Park vs design-sweep
  first; whether the two-channel symbol partition is in the first build.

## Verification policy for this analysis

⚠ **Raw p-values are useless here** — with $10^5$ pooled quantities a KS test rejects almost
anything, including two samples from the same distribution. The project already recorded this as
"significance saturates; report effect size". **Calibrate the yardstick by self-comparison** (split a
correct-by-construction sample in half) and report the statistic against that spread.

Statistics to compare are the ones $m$ depends on, not generic tree distances: the branching-time
CDF, the **clade-size distribution by depth**, the root-split distribution, and SciPhy's own
summaries (tree height, tree length, B1 balance) for free comparability with their validation.
