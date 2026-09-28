# simulator — findings

**State 2026-09-28: the tree layer and the EDITING LAYER are BUILT and VALIDATED; the editing layer is illustrated in Fig S5 (a–d) and written up in `writeup/editing_layer.tex`.** Tree library
complete (1,870 cells, 37,380 trees, zero structure faults); editing layer `src/14_editing.py`
passes all 219 closed-form cells and every hard assertion (`src/15`). Figures S1, S3, S4 built,
LaTeX write-up of S1/S3 in `writeup/`. **Dropout is the next layer; not started.**

## ⭐ Session 19 (cont., 2026-09-28) — the editing layer written up in LaTeX

`writeup/editing_layer.tex` (new, `\input` from `simulator_figures.tex`, which now covers both
layers; 24 pages, builds clean with `pdflatex` ×2). Three sections: **§4 the model and its equations**
— the recorder, the three-factor rate and $\Lambda$-time, the exact inverse of a piecewise-constant
$\lambda_0$, the route-(ii) sampler step by step with why each step is right, the channel
construction and its budget, the collision formula, what is stored, and the four closed forms;
**§5 validation** — the three configurations, the pre-fixed reading rule, the check table, and both
disclosed changes (E$_u$ and A$_r$); **§6 Fig S5** — the four panels with full captions. Intuition
first, then the equation, throughout; every symbol is in the section's own notation table.
⚠ Written for Justin to edit directly, like the tree-layer sections.
⚑ **Added after Justin's questions:** the rate subsection now states that "$\lambda_0$ averages 1"
and "$\lambda_0$ integrates to 1" are the same condition *only because time runs over $[0,1]$*; that
**$\lambda_0$ is a relative rate with no upper bound** (the constraint fixes the area, not the
height — the Fig S5d burst runs at 1.57), against **$W(t)$, which is the quantity running 0 → 1**;
and that **$\Lambda_T$ is an input** fixed before a run (5.5, from the Park mouse arms' 5.39–5.69
edits per tape), never inferred.

## ⭐ Session 19 (2026-09-28) — Fig S5: the editing layer, illustrated

`src/16_fig_editing.py` → `figures/figS5{a,b,c,d}_*.png` (+pdf), numbers `results/figS5_numbers.json`.
Justin chose panels a, b, c, then asked for d (a programmed editing rate). **The homoplasy cost of a
signal channel — collision probability vs signal share, first offered as "d" — is DEFERRED to a
fuller homoplasy test** and will live in that figure. Illustrations, not tests: where a closed form
exists it is drawn as a line over the simulated points; the formal record stays Session 18. 200
library trees of 210 cells (the validation's $\rho\times\theta$ spread), 30 tapes, $\Lambda_T=5.5$,
tapes empty at the clone founder; standard errors across trees.

- **S5a — cells inherit their ancestors' edits as a shared prefix.** One 8-cell library tree
  ($\rho$=0.25, $\theta$=0.5), 3 tapes per cell drawn as rows of 6 slots; every edit shaded by how
  many cells share it (a magnitude, so one sequential ramp — no per-branch colours). The stem's edits
  shade slots 1–2 dark in every cell; later branches add lighter, private suffixes.
- **S5b — the later two lineages split, the more slots their cells share.** One comparison per
  (internal node, tape), one cell from each side — never random cell pairs, which mostly re-compare
  the same edits (Session 18). The line is the **full** expectation: slots filled before the split,
  $\mathbb E[\min(\mathrm{Poisson}(\Lambda_Tt),6)]$, plus chance matches after it,
  $\sum_m q^m P(\text{both lineages write}\ge m)^2$ with $q=q_L=0.01664$ (lineage symbols only). All
  20 time bins sit within **$|z|\le1.04$** of it, mixed signs; chance matches add at most **0.017**
  slots.
  ⚠ **Corrected after the first render, disclosed:** drawing only the inherited depth, at bin
  centres, put the points above the line in **19 of 20 bins** (~0.03 slots early, ~0.01 late). Two
  causes, both in the drawing, not the generator: splits bunch late within a bin, where the line
  climbs ~5.5 slots per unit time; and the line omitted the chance matches the points include. Now
  plotted at the mean split time in each bin, against the full expectation.
  ⚠⚠ **A prediction of mine was WRONG.** Recommending this panel, I said the line "hits the 6-slot
  ceiling past roughly 0.7 of the experiment". It does not: cells that split at harvest share **4.8**
  slots on average. Saturation *bends* the curve: its slope (new shared slots per unit time) is
  $\Lambda_T\,P(\mathrm{Poisson}(\Lambda_Tt)\le5)$, falling from **5.5 to 2.9 by harvest**, when 47% of
  tapes are full. So a late split leaves about half as many distinguishing edits per unit time as an
  early one — fewer, not none: saturation halves the recording rate near harvest, it does not stop
  it. How much of Fig S4b's late collapse it explains, against truncation at harvest and
  single-lineage terminal branches, is still **unapportioned** (Session 17 cont.).
- **S5c — a signal pulse is written into slot order.** Share 0.02 → 0.30 → 0.02 (on from $t$=0.4 to
  0.7) at a **constant** editing rate — unlike the validation's TV set-up with its zero-rate window,
  chosen so slot order reads as time order. Read from tapes that filled all 6 slots (47% of tapes).
  Top: when each slot was written ($K_{j,6}$ of §S4.9.5 in calendar time) under the shaded
  on-window — mean write times 0.12 / 0.24 / 0.36 / 0.48 / 0.61 / 0.73, with wide, overlapping
  windows. Bottom: signal share by slot, **2.8 / 6.3 / 12.2 / 16.8 / 16.2 / 11.8%** simulated against
  2.8 / 6.1 / 12.0 / 16.7 / 16.4 / 11.9% predicted, **all within $|z|\le1.21$**. Every slot averages
  the signal over its write window, so a 30% pulse returns as at most ~17%: the recorder blurs
  time, it does not lose it. Per tip, not deduplicated — the closed form is stated per tip.
- **S5d — the editing rate can be programmed over time.** The validation's pause-and-burst schedule:
  **6.47 / 0 / 8.63 / 4.31 edits per tape per experiment** on 0–0.25 / 0.25–0.35 / 0.35–0.6 / 0.6–1.
  *Top:* the realised rate — edits in a 0.025-wide bin divided by the tape-time spent **open** (not
  yet full) in it, over every lineage alive — sits on the programmed step: **exactly 0 edits in the
  pause**, and over the 36 bins with editing $z$ has mean 0.17 and SD 1.00 (22 of 36 positive), worst
  $|z|=2.83$ (a maximum that large among 36 bins happens by chance ~16% of the time; the three
  beyond $|z|=2$ are two adjacent early bins, +2.2/+2.8, where trees have few lineages, and one late
  bin at −2.2). The same edits divided by **all** tape-time, open or full, sag as tapes fill though
  the rate holds: **8.48 → 7.17** across the burst, **3.48 → 2.31** across the last phase — the
  distinction between *the rate changed* and *the tapes filled* that panel b's corrected claim turned
  on. *Bottom:* filled slots per tape along lineages, within $|z|\le1.1$ of
  $\mathbb E[\min(\mathrm{Poisson}(\Lambda(t)),6)]$ at all 41 time points: **flat at 1.62 slots
  through the pause**, steeper through the burst, bending as tapes fill, and ending at 4.80 — the
  same as the faint constant-rate curve, since both spend the same total.
  ⚠ **Uses the simulator's true edit times**: it shows the rate can be *imposed*, not that it can be
  *recovered* from sequenced tapes. Recovery goes through branching times and meets §S4.7's confound
  with $(\theta,\rho)$ — NEXT item 3, its own proposal.

## ⭐⭐ Session 18 (2026-09-27) — the editing layer, built and validated

Approved by Justin ("PROPOSAL (2026-09-24c)", CLAUDE.md) with **$\Lambda_{\rm pre}=0$: tapes empty
at the clone founder** (Park revisited later). Theory `notes/sciphy_notes.md` §S4.9.
`src/14_editing.py` is the generator; `src/15_validate_editing.py` its validation (job 15236481,
1 min 56 s, 0.11 GB); results `results/validate_editing.json`.

**What was built.** Route (ii) of §S4.3 on library trees, parents before children, vectorised over
tapes; route (i) as the independent second implementation. $\Lambda$-time throughout, with the
exact two-step inverse of a piecewise-constant $\lambda_0$ (§S4.9.2); the channel partition
(§S4.9.3) with $A$ signal channels on their own share schedule, drawn channel-then-symbol;
per-tape speeds $r_z$; edits **stored per branch** (count, hand-off slot, symbols, times), and a
cell's tape assembled from that record alone. $\xi^L$ = Park's three mouse arms pooled by edit count
— **113 symbols, $q_L=0.01664$** (per-arm 0.0163–0.0169; the recorded "0.0170" is a rounding of these
and nothing moves). Only frequency values are read, never sequences.

**Validation — reading rule, fixed before the run.** Every cell is a deviation from a closed form
in standard errors, the se taken across **trees** (200 at $n=210$, spread over $\rho\in\{0.0005,
0.002,0.02,0.1,0.25\}\times\theta\in\{0.3,0.7\}$). No-effect value $|z|=0$. **PASS = worst
$|z|\le3.0$ over every tested cell (added checks included) AND zero hard-assertion violations.**
Three configurations: **CONST** (constant rate, signal share $p=0.05$), **TV** (a rate with a
**zero-rate window** at $t$ 0.25–0.35, and $p$ = 0.02/0.30/0.02 switching at 0.4 and 0.7 on knots
deliberately not aligned with the rate's), **RHET** (CONST with $r_z\sim$ lognormal, sd 0.5).

⇒ **PASS: 219 cells tested, 0 excluded, worst $|z|=2.67$, none above 3 (0.59 expected by chance);
0 violations in $1.9\times10^8$ hard-assertion comparisons; storage round trip identical.**

| check (what it proves) | cells | worst \|z\| |
|---|---|---|
| A — tip depth $=\min(\mathrm{Poisson}(r_z\Lambda_T),N)$ *(tree-blind)* | 21 | 1.48 |
| A′ — internal-node depth at $\Lambda_z(t_{\rm node})$, 3 time bins | 63 | 1.81 |
| A_r — per-tape mean depth under spread $r_z$ ⚠ *added* | 5 | 1.61 |
| B — signal share by slot given final depth, vs $\int K_{j,d}\,p$ *(the only order check)* | 42 | 2.03 |
| B_L — within-lineage composition, 10 equal-mass groups | 20 | 2.26 |
| C — route (i) vs route (ii), paired per tree | 38 | 2.67 |
| E — collision past the LCA, random tip pairs *(as proposed)* | 15 | 1.58 |
| E_u — the same over unique comparisons ⚠ *added* | 15 | 1.65 |
| D — hand-off slot = recomputed depth; depth ≤ N; times inside branch; queue order; symbol ids; prefix identical to the LCA | — | **0 of $1.9\times10^8$** |

⭐ **B has teeth, and passes.** Under TV, depth-6 tapes carry the signal symbol in **3.2 / 7.2 / 14.0 /
19.9 / 20.6 / 15.6%** of slots 1–6 — the on-window rising and falling along the queue — and the
closed form tracks it to ~0.001 (slot 3: 0.1404 observed vs 0.1411). B passing under both TV and
CONST means **the abandon criterion was not triggered**: §S4.2's queue-vs-bag ordering holds.

⭐ **§S4.9.4 is VERIFIED.** When both writes fall in the on-window ($p=0.30$, one signal symbol), two
lineages collide at **0.0976 against the predicted 0.0982** ($z=-0.66$); one write in the window,
0.0176 vs 0.0174; neither, 0.0165 vs 0.0164. So with one signal symbol at $p=0.3$ the collision
probability is **~5.9× the lineage channel's $q_L$**, and the effective alphabet $1/q$ falls from
**60 to 10**. This is now a quotable property of the model — ⚠ of the *model*: it says nothing
about whether ENGRAM data behave this way.

⚠ **Two changes to the validation, both disclosed.**
1. **E_u added after the 20-tree smoke run**, because E as proposed tripped there ($|z|=4.30$).
   Cause: most random tip pairs compare the *same two edits* — a slot written just below the LCA
   is shared by every tip beneath it — so per tree E rests on a handful of independent rare events
   (a match has probability ~0.018) and its t-statistic over 20 trees is not calibrated. The
   minimum-count rule overstated E's information **~11×** (115,828 expected matches at LCA depth 0
   in CONST, against 10,111 unique comparisons). E_u counts each (branch of edit 1, branch of
   edit 2, tape) once. Both are unbiased (the weighting never depends on symbols); **E as proposed
   passes at 200 trees too**, so the trip was a small-sample artefact, not a generator defect. The
   verdict's scope was widened to include E_u **before** the 200-tree run.
2. **A_r added** so the per-tape speed $r_z$, a code path the $r\equiv1$ runs never exercise, is
   tested at all.

⚠ **The yardstick is ~10% generous.** The half-split null (the same statistic between two random
halves of the trees) has SD **0.90** over 219 cells rather than 1, worst 2.25. So the across-tree se
is slightly over-estimated and every $|z|$ above is, if anything, ~10% small — the checks are
marginally less sensitive than nominal, not more permissive of a real defect at the scale that
matters.

**Cost.** 0.05 s per tree at $n=210$, $k=30$; **0.45 s** (route ii) / 0.53 s (route i) at
$n=1{,}976$, plus 0.09 s to assemble tapes; peak 0.11 GB. ⇒ decorating a whole library cell is
seconds; the tree layer, not editing, remains the expensive one. Not measured beyond $n=1{,}976$.

**What passing does NOT establish.** ⚠ A and A′ are tree-blind; B/B_L are marginals; only D and E
see the joint. ⚠ A simulator is also a model (§H.6.12): passing means it generates SciPhy's
recorder, not that the recorder is right for Park — dropout, heritable silencing (A9) and the site-6
mechanism are absent by design, and **no Park-calibrated run was made**.

## ⭐⭐ Session 17 (cont., 2026-09-24) — how far can a lineage-specific switch be detected?

**Question.** At $k\approx30$ tapes per cell, which switches in a lineage's signalling state are
detectable, and is the limit set by tape count or by timing? Scripts `src/12_switch_evidence.py`
(job 14592174, 6 min 21 s, 0.37 GB) and `src/13_fig_switch_evidence.py`; results
`results/switch_evidence.json`, plotted numbers `results/figS4_numbers.json`.

**⚠⚠ Supersedes the information-budget pair (CLAUDE.md, "PROPOSAL (2026-09-24)"), which was never
run.** Reviewed before running, it had four defects:
1. **(b)'s mock never stated $p_0,p_1$.** Back-solving its own entries: $p_1=0.20$ in the 2× and 10×
   columns but $0.24$ in the 24× column. At $k=30$, 24× that is **3.36 vs 2.80 nats** — either side
   of the 3-nat decision line, so the decisive row was set by an unstated, inconsistent parameter.
2. **(b) scored one branch in isolation**, which is not a bound: a state persists across branches and
   every sublineage records it independently, so evidence accumulates. "Only a negative is decisive"
   did not hold for it.
3. **(a) was called a lower bound (bigrams excluded) in its §2 and an upper bound in its §3** — both
   true, so neither outcome was decisive; and it addressed population $\xi(t)$, which §S4.7 had
   already called over-determined.
4. **(a)'s kernel rows were normalised to sum 1**, so singular values scale as $1/\sqrt G$ while the
   noise bar does not — $R$ would move with the grid $G$ by construction, tripping its own abandon
   criterion for a definitional reason.
⚠ Two corrections to my own reasoning in the review, recorded: (i) "the lower of $\kappa_{\rm on}$,
$\kappa_{\rm off}$ controls the call" is **wrong for detection** — power comes from the on-direction;
the off-direction governs affirming *no* switch. (ii) "a single-lineage short switch at $k=30$ sits
right on the line" was measured against 3 nats, not the searched bar (7–13 nats); it is well below.

### The design (approved by Justin; he added $n=8$ and $\rho$ up to 0.8, and a $p_1$ scan)

A **switch**: at a point (branch $x_0$, time $s_0$) a cell turns a pathway on; every descendant
inherits it; all turn off at $s_1=\min(s_0+D,1)$, $D=\Lambda_e/\Lambda_T$. Switch locations are
**uniform along the tree's total length, stem included** (a constant per-lineage switch rate).
Expected evidence in favour of the switch, nats:

$$\mathcal I=k(1-d)\sum_b W_b\,\mathrm{KL}(\bar p_b\|p_0),\quad W_b=E(e_b)-E(a_b),\quad
E(s)=\mathbb E[\min(\mathrm{Poisson}(\Lambda_Ts),N)]=\sum_{j=1}^{N}P(j,\Lambda_Ts)$$

- $s=t/T\in[0,1]$; branch $b$ spans $[a_b,e_b]$; $W_b$ = expected edits per tape written on it.
- $\bar p_b=p_0+(W^{\rm on}_b/W_b)(p_1-p_0)$ — the signal share of an edit known only to lie on $b$
  (**branch-only**, the realistic mode). **Timing-known** (upper bound):
  $\mathcal I=k(1-d)\,\mathrm{KL}(p_1\|p_0)\sum_bW^{\rm on}_b$.
- KL is Bernoulli per edit (signal symbol or not): ENGRAM moves $\xi$, not $\lambda$, so the edit
  count on a branch carries no state information. $p_0=p_1/F$, $F$ the fold change.
- **Callable** ⇔ $\mathcal I\ge3+\ln H$, $H=(2n-1)\times4$ candidate switches: **7.09 nats at $n=8$,
  12.67 at $n=1{,}976$**. The pseudo-LR has $\mathbb E_{H_0}[\mathrm{LR}]=1$ exactly, so the tail
  bound holds; "callable" means expected evidence clears the bar (~50% power).
- Fixed: $\Lambda_T=5.5$ edits/tape, $N=6$, $d=0.44$, $\theta=0.5$. Scanned: $n\in\{8,25,64,210,626,
  1976\}$ × $\rho\in\{0.002,0.1,0.25$ (library)$,0.5,0.8$ (simulated fresh by `04.run_cell`, into
  `results/tree_library_hi_rho/`, gitignored)$\}$ × $\Lambda_e\in\{0.25,0.5,1,2\}$ edits/tape (4.5,
  9.1, 18.2, 36.4% of the experiment) × $F\in\{2,5,15,24\}$ × $p_1\in\{0.05,0.1,0.2,0.3,0.4\}$; 20
  trees × 1,000 systematic switch locations per cell; $k$ evaluated afterwards, since $\mathcal I$ is
  linear in it.
- ⚠ Library trees store the root's own branch as 0, so absolute times are rebuilt from the tips
  (all at $s=1$); the stem $0\to t_{\rm root}$ is recovered and included.

**✅ Verification, run before any result.** $E(s)$ vs $10^6$-draw Monte Carlo at four times:
$|z|\le1.24$. Closed-form $\mathcal I$ vs the mean realised LLR over 40,000 simulated tapes on a real
$n=25$ tree (stem, internal and terminal switches × two durations × three $(F,p_1)$ × both modes):
**worst $|z|=2.43$ over 36 comparisons**, what 36 draws give by chance. Branch-only ≤ timing-known
asserted on every switch (KL is convex in its first argument). Tree-to-tree SD of the callable
fraction is at most 0.10 in any cell, so the abandon criterion (conclusions changing between
replicate trees) was not tripped.

### Result 1 — the agreed headline: fraction of ALL switch locations callable (Fig S4a)

*Reading rule:* each entry is the fraction of switch locations callable at $k=30$, $p_1=0.3$,
branch-only timing, as a **range over the 30 $(n,\rho)$ settings**. No effect = 0; the pre-set
positive bar was 50%; the pre-set decisive negative was "<10% at every $n$ even with timing known, at
24×, $\Lambda_e=2$".

| duration $\Lambda_e$ (edits/tape · % of experiment) | 2× | 5× | 15× | 24× |
|---|---|---|---|---|
| 0.25 · 4.5% | 0 | 0 | 0 | 0 |
| 0.5 · 9.1% | 0 | 0 | 0.1–0.4% | 0.5–1.4% |
| 1 · 18.2% | ≤6% | 0.7–2.3% | 5–19% | 9–36% |
| 2 · 36.4% | ≤6% | 8–27% | 20–58% | 23–63% |

⇒ **Decisive negative NOT triggered** (timing-known at 24×, $\Lambda_e=2$: 23–75%). **Positive bar met
only partly** (2-edit switches in some settings; never at 1 edit). **Marginal, and duration-driven**:
switches of ≤0.5 edits/tape, and 2-fold switches of any duration, are effectively invisible.

### Result 2 — ⚠ the same result read by start time (Fig S4b). Found AFTER the run.

⚠ **The headline metric is dominated by late switches.** Uniform-along-length puts **46–82% of switch
locations on terminal branches** and **43–92% in the last 30% of the experiment** — where almost
nothing is callable. The by-time breakdown was in the pre-specified output but **not** in the mock or
the reading rule; the statistic was not changed, and both readings are recorded. Which question is
the right one is Justin's call.

At 2-edit switches, $F\ge15$, branch-only: switches starting in the **first 60%** of the experiment
are callable in **≥89% of locations in every one of the 30 settings** (100% at 24×); those starting
in the **last 20%** in **at most 22% (15×) / 31% (24×)**. At 1-edit, 24×, callable falls steadily
with start time: median over settings **0.82** in the first decile, 0.50 at 0.5–0.6, 0.37 at
0.7–0.8, 0.00 in the last. ⚠ The saved output cannot
apportion the late failure between truncation at harvest, single-lineage terminal branches, and
tape saturation.

### The other findings

- **Clade size barely matters:** a 2-edit, ≥15× switch needs only **2–5 sampled descendants** ($m^*$ at
  $k=30$) — ~2 at $\rho=0.002$, ~5 at $\rho=0.8$ (at high capture small clades are young). Terminal
  (single-cell) switches rarely clear the bar.
- **Timing precision matters little:** branch-only keeps **83–97%** of the timing-known fraction at
  2-edit switches; only at 1 edit does it cost much (**24–88%** kept), mostly at low $\rho$.
  ⇒ **the tree side is not the bottleneck.**
- **Tapes:** $k$ at which half of switch locations are callable (24×, $p_1=0.3$, branch-only):
  **88–176** at $\Lambda_e=0.5$, **34–120** at 1, **18–120** at 2 (251–503 at 0.25). Park has 166.
- **$p_1$ is a ~10× lever:** 2-edit, 24× switches go from **2–5%** callable at $p_1=0.05$ to
  **23–63%** at $p_1=0.3$. The ~30% shared-tape ceiling (§I.7) binds.
- **The search penalty is large:** against a flat 3 nats (location known in advance), 1-edit, 24×
  switches are callable at **51–77%**, not 9–36%. Knowing which clade to test — e.g. from a
  transcriptomic readout on the tips — would recover much of it.
- **$\rho$ and $n$ are second-order.** Larger clones score lower mainly because $\ln H$ rises and more
  of their length is late terminal branch; in absolute terms they hold more callable switches.

### What it means for building (recommendation, 2026-09-24)

1. **Build the state layer after the editing layer**, with fitness-neutral states painted onto library
   trees; validate against the regime the closed form says is recoverable — long ($\ge$1–2
   edits/tape), strong ($\ge$15×), early switches in clades of $\ge$2–5 cells — to see whether tree
   error and homoplasy erase it. The multi-type birth–death can wait.
2. **Deprioritise tree-side timing refinement**; branch-level timing costs little.
3. **Short and late switches are design limits (purpose IV)**, not modelling targets: short switches
   need ~90–180 tapes and high $p_1$; the last ~20–30% of the experiment is out of reach at this
   editing depth at any capture fraction.

**Caveats, by direction.** Optimistic: the tree is known and homoplasy-free (the one serious one).
Conservative: Bonferroni over correlated switch positions; $(1-d)$ for internal branches (an edit is
seen if *any* descendant keeps the tape — at most $1/(1-d)\approx1.8\times$); site order unused.
Simplified: one $\Lambda_T$ with no per-tape rate spread; a fixed-length switch ending everywhere at
once, truncated at harvest.

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
