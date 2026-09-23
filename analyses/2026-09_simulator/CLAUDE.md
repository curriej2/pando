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

## ⇒ NEXT, in order (updated 2026-09-23; the 2026-09-17 list is superseded)

1. **The editing layer** — decorate library trees with sequential edits at measured per-tape
   rates, site 6 handled separately, depth 0 censored; editing and dropout must be fitted
   **jointly** (park-compatibility CLAUDE.md).
2. Then, with edits: tape saturation vs $L(t)$ — the "best-recorded window" hypothesis (README
   Session 16, point 4); a proposed figure is $L(t)/n$ with the tape-saturation profile overlaid.
3. Survivorship bias needs signal-dependent $b$/$\delta$ — a model extension, not yet designed.
4. Open, lower priority: check A at a Park-scale point; refit the attempts envelope with the exact
   $1/P(K=n)$ (README correction 3); S2 validation figures were proposed and Justin declined them.
5. ⚠ Owed to Jihye Park: total viable cells or tumour mass per organ at day 45 (mouse $\rho$).

## Verification policy for this analysis

⚠ **Raw p-values are useless here** — with $10^5$ pooled quantities a KS test rejects almost
anything, including two samples from the same distribution. The project already recorded this as
"significance saturates; report effect size". **Calibrate the yardstick by self-comparison** (split a
correct-by-construction sample in half) and report the statistic against that spread.

Statistics to compare are the ones $m$ depends on, not generic tree distances: the branching-time
CDF, the **clade-size distribution by depth**, the root-split distribution, and SciPhy's own
summaries (tree height, tree length, B1 balance) for free comparability with their validation.
