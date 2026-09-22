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

### ⚠⚠ Propose before you compute — the approval gate

**Do not write an analysis script, and do not launch a job, until Justin has said go.** Brainstorming
is collaborative and welcome; execution is gated. Never put the proposal and the run in the same turn,
and never treat "that sounds interesting" as approval.

A proposal is short — aim for under ~25 lines — and has four parts:

1. **The question, in one sentence.** Specifically: *what would change about the project's
   conclusions* depending on how it comes out. If nothing changes either way, say so and drop it.
2. **The math, written out.** The estimator or test statistic in symbols, its null hypothesis, and
   the assumptions it needs to be valid. Define **every** symbol against the house register (§H.6.0)
   — including ones defined in an earlier session — and give units. If a quantity is new, name it
   and say what one unit of it means.
3. **What could make it wrong.** Name the confound the design defeats *and* the one it does not.
   In this project the recurring one is circularity: relatedness is read from the same edit data
   whose readability dropout controls.
4. **The output, structurally — a mock, not a description.** Sketch the actual table or panel with
   placeholder numbers: exact column headers with units, how many rows and what indexes them, the
   grid a curve runs over, the axes of a figure. Then state **the reading rule**: what the
   no-effect value is, which direction means the effect is present, and what magnitude would count
   as decisive versus marginal.

**Why this is a hard rule.** A result whose test Justin did not understand *before* it ran costs more
to read than the run saved, and it tends to be overturned later — the ⚠⚠ corrections littered
through this file (self left in the null, stratifying on an estimate, decile instead of rank
matching, $\Lambda_{\rm soft}$ used as a detector) are all cases where the reasoning was finished
after the numbers arrived instead of before. Approval is cheaper than retraction.

- If the design turns out to be wrong **mid-run**: stop, say what broke, and re-propose. Do not
  patch the statistic and keep going.
- **Exempt:** cheap reversible probes whose purpose is to *check* rather than to conclude — a row
  count, whether a column exists, an assertion on a table already in memory, reading a log. The gate
  applies to anything whose numbers would be quoted, plotted, or written into the notes.

### Reporting numbers — definitions travel with them

- **No large tables of bare numbers.** Every number that appears carries: what it is, its units, its
  null / no-effect value, and what a larger value would mean. A reader who has forgotten the variable
  must still be able to read the line. Loosely-defined symbols are the main way output becomes
  unusable here.
- Prefer **one number carrying the point inside a sentence** over a grid of numbers. Where a table
  really is the right shape, keep it small, label every column with units, and put the reading rule
  immediately above it — not in a paragraph further down.
- Quote a threshold-dependent number **with its threshold**, an effect **with its comparison**
  (observed − null), and a per-$k$ or per-clone number **relative to the denominator that matters**
  ($k/n_C$, not $k$). These three have each already caused a retraction in this project.

### Caveats are wanted — raise them fully

- A side analysis that surfaces a real caveat, bug, or confound is **not** noise: report it
  unprompted and at length. Say what was found, how, **which recorded numbers it affects**, and
  whether it moves a number or overturns a claim we have been quoting. Mark it ⚠ in the notes, as
  the existing record does.
- Keep those two categories distinct in the write-up. "This shifts 16.3 to 15.5 nats" and "this
  means the gradient we reported does not exist" need different headlines.

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
| `2026-08_park-compatibility` | Does Park's cross-tape character compatibility support the perfect-phylogeny route (§D.4b)? | **Diagnostics complete.** Mouse3 94.66%/63.56% (spread +31 pts), Initial 91.91%/80.60% (+11 pts) ⇒ **dropout binds, row A6**. $C=0.57$ on Mouse3. Subclone ground truth **passed**. ⚑ **Strategic pivot: skeleton is a step sideways — only 9 of 2,547 clones exceed 1,000 cells, so likelihood inference is already in range. Figure programme under way: figs 1–3 done.** ⚑ **Dropout characterised (fig 3): two axes — the tape axis is larger ($\rho_{\rm tape}=0.25$ vs $\rho_{\rm cell}=0.13$), reproducible ($r=0.997$ between libraries) and informative about edit depth, but only through a removable 15% of tapes; the cell axis is a pure QC artefact and carries no editing information, so it can be marginalised.**
⚑ **Row A9 measured (2026-09-02): tape loss IS heritable, and heritable below the clone — so it is a
Dollo character needing one absorbing state in the pruning, not a reason to reach for SBI.**
⚑ **B = 1,000 permutations launched 2026-09-03** over all three layers × five arms; clone-wide layer
already back at $p=1/1001$ in 5/5 arms. ⚠ Two orthogonal follow-ups are specified but not run: a
`--lam 4` re-run with a size-stratified null (the floor, not $B$, is what hides weaker events), and
the **symbol-depletion test** of row A9's untested *trans* limb. |

| `2026-09_simulator` | Build the forward simulator (birth–death tree → sequential editing → dropout) for the homoplasy null and the design sweep. | **Design settled and verified; nothing built.** Tree model: shape/time factorise, clone sizes taken from data not simulated, two effective parameters not three. BDS branching-time CDF closed-form and invertible. Editing layer measured (per-tape rates real, site 6 anomalous, depth 0 censored). ⚠ Editing and dropout must be fitted jointly. |

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
   **A9 lineage test MEASURED (2026-09-02); Fig 4 panels a/b built.** Related cells lose the same
   tapes — **~1,980 individual clones, essentially all positive**, at 31–126% of a marginal-matched
   heritable control; a **monotone gradient over subclade depth in 5/5 arms**; and a **discrete**
   loss rather than a graded rate. ⚠ Its honest null turned out to be a permutation, *not* the
   simulator — that correction unblocked the work. ⚠ Pooled $\rho$ weights clones by $n_C(n_C-1)$ and
   **hid** the signal; use equal-clone weighting.
   **Event catalogue built (2026-09-03), FDR ≈ 0 in 5/5 arms.** ✅ Mechanism verified from the
   paper: the tape is read from **cDNA**, and pegRNA+tape are **co-integrated**, so silencing an
   integration removes both the readout and that symbol from the writing pool — re-explaining Fig 3b
   ($\beta_z$ = the integration's expression level) and Fig 3d. **Zero events in clones <20 cells**
   (power, not biology); where detectable, **6.5–19%** of missing entries sit in a called loss, plus
   a **clone-wide layer of 7–33%** — the same mechanism before vs after the day-3 bottleneck.
   $\hat\pi$ median 0.99–1.00 in 4/5 arms, so most losses are complete; 8–32% partial (an upper
   bound — Subclone's 57% is clade-resolution artefact, since prefix clades stop at depth 4 while its
   clones run 200–11,000 cells).
   **`MAX_D=6` test RUN (2026-09-03): partly artefact, partly real.** At full recorder resolution
   median $\hat\pi$ reaches 1.000 in all three mice and the partial fraction falls (Mouse3
   12.7%→5.3%, Pre-TX 15.1%→6.8%, Mouse1 22.9%→19.8%, Mouse2 37.5%→26.0%, Subclone 85.2%→37.7%).
   ⚑ **The residual partial fraction is monotone in max clone size** (210→5.3% … 10,996→37.7%) — the
   recorder runs out of resolution before the tree does, and where six levels suffice it converges to
   **5–7%**. ⇒ **the graded appearance is very largely clade coarseness**, so a per-tape absorbing
   state is the right and largely sufficient extension; the lineage-varying-rate case is weaker than
   the depth-4 numbers implied.
   **✅ B = 1,000 COMPLETE 2026-09-04, all 15 configurations.** $p=1/1001$ in all twenty runs.
   ⭐ **Quote the null MAX, not the p:** the most extreme of 1,000 random labellings gave **111**
   candidates in Subclone against **279,973** real (118–2,522× across arms); the null is right-skewed
   6–100×, which $B=3$ could not have seen. Every hard-catalogue event clears $q\le1.9\times10^{-4}$.
   **Cell × tape entries inside called blocks** (lead with this, not event counts): Subclone 266,688
   / Mouse2 21,277 / Pre-TX 22,256 / Mouse1 14,770 / Mouse3 2,916. ⚠ One marginal result to disclose:
   clone-wide Mouse2, FDR 4.36%, max $q$ 0.044.
   ⚠⚠ **CORRECTION: I predicted the thresholds would fall and they did not** — Mouse3 soft 16.3→16.2
   (d4) / 15.5 (d6), everything else stayed on the 10-nat floor. The $B=3$ null was unbiased in the
   *mean*; what it missed was the *tail*. $B=1{,}000$ bought resolution on $p$ and the null's shape,
   **not one new event and not one moved threshold** — and it sharpens the `--lam 4` case, since the
   threshold now sits *on* the floor in 14 of 15 configurations.

   *(launch record)* **B = 1,000 launched 2026-09-03** across all 15 configurations (5 arms × {`40` hard, `43` soft
   d4, `43` soft d6}) at 20 array tasks each, plus 5 direct `42` runs; `47_submit_B1000.sh` →
   `46_perm_merge.py` → `48_collect_B1000.sh`. Parts store **count vectors on a fixed grid**, never
   candidate lists, and permutation $b$ is seeded from $(\text{SEED},b)$ so slices are addable in any
   order. **The five `42` clone-wide runs are in: $p=1/1001$ in every arm** — not one of 1,000
   permutations came close — with loss counts unchanged from $B=200$.
   ⚠⚠ **$B$ is NOT the detection floor.** Three knobs: $B$ (resolution of $p$), the FDR *threshold*
   (falls where a 3-permutation null pushed it up — Mouse3 soft 16.3 nats should drop toward 10), and
   the scan **floor** `--lam`, hard-fixed at 10 nats, beneath which nothing is ever collected and
   which this run's grid cannot see. ⇒ **follow-up: `--lam 4`, which subsumes the floor-10 run**,
   plus a **null stratified by clade size** and per-combo exceedance counters (both only worth it at
   the lower floor — see the analysis README; per-combo p at floor 10 censors every real event at
   $1/1001$, and $p\ge1/\binom{n_C}{m}$ regardless of $B$).
   **⭐ NEXT: (i) read the collector's FDR curves, then launch the `--lam 4` run.
   (ii) panels c (a worked example inside one large clone) and d (the catalogue).**
   **The remaining figs are still blocked on the simulator**: birth–death
   tree at Park clone sizes, sequential editing at measured $\lambda$/$\xi$, $N=6$, $k=166$, with
   **dropout as a switchable layer**. That simulator also supplies the **homoplasy null** — still the
   single biggest gap in the "dropout not homoplasy" argument.
0a. **✅ SETTLED 2026-09-07 — the defensible catalogue, and the criterion is an ABSOLUTE budget
   not a rate.** Low-floor ($\Lambda\ge2$) + clade-size-stratified + expected false $\le2$ per
   stratum: **Subclone 6,597 events / 20.07% of all missing · Pre-TX 2,019 / 2.46% · Mouse2 452 /
   7.02% · Mouse1 424 / 4.34% · Mouse3 132 / 3.15%**, $\le11$ expected false candidates per arm,
   every event $q\le2\times10^{-3}$. **Quote these.** Against floor 10 that is only 1.1–1.8×
   (0.8× for Subclone) — the 3–14× the 5% target advertised was false positives.
   ⚠⚠ **Why a rate fails at a low floor:** FDR is computed on *candidates*, the reported quantity is
   *events* after the overlap collapse. Harmless at floor 10 (2.15 expected false vs 320 events),
   fatal at 5% (5,515 vs 2,083), and count vectors cannot dedup permuted sets to fix it. Use
   `--budget`.
   ⚠⚠ **Two corrections to my own reasoning, both recorded:** (i) "the floor is conservatism not
   rigour / loose lower bounds" was right in direction, badly wrong in magnitude — I measured the
   distance to the cliff in FDR units when it lives in nats; defensible thresholds are **5.9–12.4
   nats**, essentially where the arbitrary floor sat, so floor-10 was accidentally near-right.
   (ii) I declared the 4–5 cell stratum "DEAD in all arms" from the arm-*median* $\tilde p$; it
   calls 18–635 events per arm, because called small-clade events sit at $\tilde p$ 0.05–0.13 vs
   medians 0.12–0.42. ⇒ **"a small clade is detectable only on a tape that should have been there"**
   — capture-independence, not a dead zone.
   ⚑ **`MAX_D=6` survives**: partial fraction at the budget is 8.9–46.3%, unchanged or lower than
   floor 10. The rise at 5% was mechanical — a partial loss scores lower $\Lambda$, so any threshold
   drop inflates it. **Always quote the partial fraction with its threshold.**
   ⚠ Still open: `42` (clone-wide) is unstratified, and clone sizes vary far more than clade sizes —
   its Mouse2 FDR 4.36% is the last marginal number. Full audit + results in README.

0a-old. *(superseded)* **Event calling was pinned to an ARBITRARY floor.** Two thresholds: the scan floor `LAM0 = 10` nats (hard-coded, nothing
   below it is ever collected) and an adaptive `LAM` (smallest grid point with FDR $\le5\%$).
   **The adaptive step is degenerate in 13 of 15 configurations** — FDR at the floor is 0.002–0.019%
   for the hard catalogue, **262–2,511× below target**, so `LAM` snaps to the floor. Binds only for
   Mouse3 soft (16.2 / 15.5). ⇒ event counts and the dropout share are set by the constant, and the
   **arm ranking is not stable**: Pre-TX is 2nd at $\Lambda\ge10$ (1,783) and 4th at $\ge20$ (110).
   **Never quote a single count or fraction — quote the curve.** ⚑ The floor is *conservatism*, not
   rigour (we sit hundreds of times inside the false-positive cliff), so the numbers are loose lower
   bounds and `--lam 4` is simply applying the stated criterion. ⚠ Do not inherit 5%: with $10^5$
   candidates it admits ~14,000 false events in Subclone — prefer an absolute target (expected false
   events $\le10$). Threshold-independent claims that survive: the 118–2,522× null comparison,
   $q\le1.9\times10^{-4}$ per event, $\hat\pi$ median 1.000 vs expected 0.11–0.41. Full audit in
   README "How events are actually called".
0-CURRENT. **⭐⭐ 2026-09-20/22 — THE SIMULATOR IS BUILT, VALIDATED, AND THE TREE LIBRARY EXISTS.
   ⭐⭐ THE BDS DENSITY IS DROPPED.** `analyses/2026-09_simulator/src/` holds `01_bdtree.py`
   (Gillespie forward + reject → genealogy → $\rho$-sample → reconstruct, plus a naive oracle),
   `02_validate_tree.py`, `03_cost_curve.py`, `04_tree_library.py`, `05_submit_library.sh`.
   **✅ VALIDATION PASSES EVERYWHERE**: 84,000 trees, **zero structure faults**, 3 turnovers × 3 $n$
   × 2 capture fractions plus structural checks at $n$ = 210/1,024/3,387. Pruning load spans
   **3.4–65.2** nodes per branch point (19×), so degree-2 suppression is tested across the regime.
   **✅ COST MEASURED ⇒ the decision quantity is the ATOMIC UNIT** (one tree of the arm's largest
   clone), **not** core-hours, which parallelise over clones and replicates. Worst case in the whole
   project is Subclone at Park-like $\rho$, **70.7 h**, inside the 7-day cap. ⇒ §S3.3.5 route (c) and
   the owed branching-time density verification are **no longer needed**. ⭐ The cost model validated
   **out of sample** (0.91× on seconds, 1.00× on peak live lineages at an $n$ not used to fit it).
   **✅ TREE LIBRARY: 29,780 trees, zero faults, 97 MB** (gitignored). Stocks the **range** of clone
   sizes, not Park's clone list — 20 log-spaced sizes pinning the five arm maxima plus every integer
   2–64 (97.1% of Park clones, median 8). Stores `parent` + **branch lengths** (not node times) and
   seed provenance; replay is **bit-identical**.
   **⭐⭐ $\rho$ DOMINATES $\theta$ BY 3–5×** on coalescent depth (+0.27…0.34 over $\rho$ vs
   +0.04…0.11 over $\theta$, monotone 20/20 and 5/5). ⇒ pinning the mouse capture fraction matters
   more than pinning turnover.
   ⚠⚠ **$\rho=8\times10^{-4}$ was SciPhy's HEK293T constant, never a Park measurement.** Corrected:
   Initial ~0.24 (upper bound), Subclone 0.1–0.3, mice **unknown** — IVIS gives relative flux, not a
   cell count. ⚠ SciPhy found a nominal sequenced/population ratio over-stated effective $\rho$ by
   **4–16×**, so treat every nominal value as an upper bound. **⇒ ASK JIHYE PARK: total viable cells
   per dissociation, or tumour mass, per organ at day 45.**
   ⚑ $\theta$ from the pool's realised doubling (8,000→160,000 in 10 d = 55.5 h) against H1299's
   intrinsic 22–30 h gives **0.46–0.60**, the only edit-independent handle; hence a grid, not a point.
   ⚠⚠ **$\theta$ cannot be fitted before the editing layer exists** — every observation of
   within-clone structure in Park is edit-mediated (ClonalBC gives membership, no internal structure).
   ⚠ **EIGHT defects this thread, every one in a test, reading rule or projection, none in the
   simulator** — see `analysis_log.md` sessions 14–15. The recurring shape: a statistic fitted or
   calibrated in one regime and applied in another.
   **⇒ NEXT:** read arrays `14228126` ($\theta=0$ Yule reference) and `14228127` ($\rho$ extended to
   0.05/0.02/0.0005), compute the **Yule-differenced LTT** to decide whether turnover survives at the
   mouse capture fraction, then **the editing layer**.

0-PREV5. **⭐⭐ 2026-09-15/17 — THE SIMULATOR THREAD IS OPEN. Tree model SETTLED, editing
   layer MEASURED, nothing built yet.** New directory `analyses/2026-09_simulator` (its CLAUDE.md is
   the full design); theory in `notes/sciphy_notes.md` **§S3**; running record `analysis_log.md`
   session 13. **`src/` is deliberately empty — no script has cleared the approval gate.**
   **Purposes I and IV lead** (Justin's call): *adjudication* — quantify the consequence of what is
   already established — and *design*, the $(p,\lambda,m,k,j,\ell)$ sweep.
   ⚠⚠ **TWO SCOPE CORRECTIONS, both narrowing.** (i) I proposed compute for a **dead question** — an
   Initial compatibility rerun and a Subclone pricing — serving the max-compatible-set route that
   died in the 2026-09-01 pivot. **Both withdrawn.** Compatibility survives only as a free end-stage
   check, since Mouse3 is complete on disk. (ii) **The homoplasy null is narrower than framed**: the
   mouse arms are already measured as near-homoplasy-free (0.03–0.11 recurrences per level-0 prefix
   node, 0.4–1.1% of nodes with any), so it is load-bearing only for **Subclone** — where the
   strongest dropout signals also live, so the two processes are entangled — and for the design
   sweep. ⇒ **fig 5 to be renamed and rescoped**; incompatibility is a perfect-phylogeny quantity and
   the likelihood route wants *reconstruction accuracy*.
   **✅ SciPhy Session 3 read** (the validation methods, queued since session 2a). They simulate 100
   birth–death trees, draw clock rate $\sim\mathrm{LogNormal}(-2,0.5)$ and $\xi\sim\mathrm{Dir}(1.5)$,
   decorate 10 tapes, and infer with the *same* distributions as priors (that is SBC — it tests the
   implementation, not realism, so **do not inherit it**: our null needs parameters *fitted* to Park).
   ⭐ **They already simulate both our dropout axes by name** — heritable silencing + sequencing
   dropout — but **filter to a complete submatrix rather than modelling it**, which is exactly the
   row-A6 gap. ⚑ Reusable: **PI distance is "particularly sensitive" to tape loss; wRF rises only
   weakly** — tells us which metric will show the effect.
   **⇒ Language: PYTHON.** Their simulator is 260 lines of Java and takes the tree as an *input*;
   all the optimisation (caching 3.5×, threading, 8× overall) is in the **likelihood**, a different
   workload. Drive BEAST externally only if SciPhy *inference* is wanted.
   **⭐⭐ TREE MODEL, both parts verified.** (a) **Shape and time factorise** — conditional on tip
   count the topology is independent of $b,\delta$ (checked at turnover 0 / 0.14 / 0.75; root splits
   match Yule–Harding to ≤2.0 se). Shape is parameter-free; $b,\delta,\rho$ act only on the branching
   *times*. ⚠⚠ "Random binary tree" is ambiguous: PDA ≠ Yule, PDA is far more imbalanced, and clade
   sizes set $m$. (b) **Clone sizes CANNOT come from one homogeneous birth–death** — its geometric
   tail predicts $2.8\times10^{-14}$ (Mouse3) to underflow (Mouse2) clones as large as the largest
   observed; one is observed in every arm. ⇒ **do not simulate clone sizes; take them from the data
   and condition the within-clone tree on observed size.** ⚠ $(b,\delta,\rho)$ are not separately
   identifiable — **fit two effective parameters, not three**; with shape free the tree contributes
   essentially **one number**, how deep the coalescences sit.
   **⭐ The BDS branching-time density is closed-form and analytically invertible**:
   $G(x)=[h(x)-\rho]/[h(t_{or})-\rho]$ with $h=1-p_0$, because $p_1=h'/K$. Both $p_0$ and $p_1$ fall
   out of the generating function at $s=1-\rho$ — no new machinery. ⇒ exact $O(n)$ sampling, which is
   what makes a 10,997-tip clone feasible. ⚠ Verified to $\sup|{\rm emp}-G|=0.0080$ over 7,867 trees;
   **not yet settled** because within-tree dependence makes the naive KS critical value inapplicable.
   **⭐ EDITING LAYER MEASURED (all five arms, cached matrices).** ⚠ **Depth 0 is censored** —
   "recovered but unedited" is indistinguishable from dropout, so fits must be zero-truncated; at
   Initial's $\hat\Lambda\approx2.82$ that is up to **a quarter of its apparent dropout**. **A single
   rate is rejected**, $\chi^2/\mathrm{df}$ = 15,836–86,154. **Per-tape rates are real** (mean depth
   1.33–5.54 sites on Mouse3; $\mathrm{corr}$ with dropout −0.33 to −0.63, **unchanged** in the
   best-captured decile, so not a capture artefact — and SciPhy fits per-tape clock rates too).
   ⭐ **Site 6 is mechanistically different and survives the per-tape fit** (obs/exp at depth 5/6 =
   2.09/0.82 on Mouse3), matching the independently-measured elevated $q$ at site 6 in 5/5 arms.
   ⚑ **Initial vs the mice is a natural two-timepoint experiment** ($\Lambda$ 2.87 vs 5.39–5.69).
   ⚠⚠ **Editing and dropout are NOT separately fittable** — missingness is informative for depth, so
   fitting $\Lambda$ on clean reads conditions on a depth-dependent selection and is biased upward.
   This **revises** the "calibrate editing first, then dropout" ordering proposed earlier the same day.
   **⇒ NEXT:** (1) self-calibrate the density check by splitting the existing 7,867 trees — no new
   forward replicates; (2) build the Option C sampler, validating density and sampler *separately*;
   (3) the joint editing/dropout fit. ⚠ Each needs a proposal under the approval gate.

0-PREV4. **⭐⭐ 2026-09-14 — ⭐ THE SILENCING EVIDENCE IS SUFFICIENT (Justin's call). FIGURE
   PROGRAMME CLOSED; THE SIMULATOR IS NEXT.** `2026-08_park-compatibility`, README and its CLAUDE.md.
   **Four figures carry it**, five arms each: `fig6a` the ladder (held-out gain vs parameters spent),
   `fig6b` $Q$ (the lineage share of what the technical model left), `fig6c` per-tape unanimity
   (missing vs present against the diagonal), `fig6d` the companion (capture matched, relatedness not).
   **Ladder:** Subclone **+8.31 nats/cell, $Q$=15.9%** · Pre-TX **+4.59, 7.1%** · Mouse2 +2.47 ·
   Mouse3 +2.19 · Mouse1 +1.48 — positive 5/5 at 16–90 se, one added parameter against 3,232–64,975.
   **Per-tape unanimity** ($k$=5, common tapes): Pre-TX **9.77× missing vs 1.22× present**, Subclone
   **3.84× vs 1.15×**, mice ~1.2× vs ~1.03×.
   ⚠⚠ **Four corrections this thread paid for, all recorded**: the null must exclude self; a cell
   could be its own nearest relative; stratifying on an ESTIMATE does not hold the truth fixed (a
   no-structure simulation still rose 14.9→28.2%); and capture matching must be by rank, not decile
   (deciles inflated folds up to 20%).
   ⇒ **NEXT: the simulator.** Target the ladder's 4-vector including the SIGN of $M_2$, the
   variogram's slope and convexity, the nats-by-$k$ inversion, and the unanimity fold. ⚠ These are
   JOINT constraints on (tree, editing, dropout) — the variogram's $x$-axis is the editing process —
   so calibrate tree and editing first, then use these to falsify rather than to fit.

0-PREV3. **⭐⭐ 2026-09-11 — THE NESTED LADDER IS BUILT (`77`), AND IT FOUND A BUG IN THE NULL
   THAT ALSO AFFECTS `74`.** `2026-08_park-compatibility`, README "The nested ladder".
   **The design, agreed with Justin over three rounds:** fit **within clone** (so clone, batch,
   harvest and mouse are held fixed *by construction*, not by trusting a $\gamma$ term), hold out
   **entries** not cells (so a rung can LOSE — the only way the comparison means anything), and fit
   $M_4$ **two-stage** on a frozen offset (so $\Delta_4$ is a conservative lower bound and "one
   parameter" is literal). Justin's simplification — fit $\alpha$ and $\beta$ inside the clone
   rather than adding $\gamma_{Cz}$ — is **exactly equivalent**, since the likelihood factorises over
   clones, and it fixed an identifiability error of mine (one constraint is needed **per clone**).
   ⚠⚠ **THE NULL MUST EXCLUDE SELF.** The $M_3$ score equation forces $\sum_{c\in C}r_{cz}=0$, so a
   self-excluding neighbour set carries $-1/(n_C-1)$ of the cell's own residual. `74`'s null permutes
   residual ROWS across fixed neighbour SLOTS — a uniform subset of ALL $n_C$ cells INCLUDING $c$,
   expectation 0 — so it carries no leak and obs−null cannot cancel it. Negligible at $n_C=3387$,
   **fatal at $n_C=31$**: Pre-TX returned a negative lineage rung until fixed. The null is now
   **$k$ random clone-mates, self excluded**. ⚠ `74`'s small-clone numbers are in question.
   **⭐ RESULTS (floor 100, $k=20$, nats/cell obs−null):** Subclone **+8.31±0.09**, $Q$=**15.9%** ·
   Pre-TX **+4.59±0.12**, 7.1% · Mouse2 +2.47, 4.4% · Mouse3 +2.19, 3.8% · Mouse1 +1.48, 2.4%.
   **Positive 5/5 at 16–90 se**, $w$ = +0.78…+1.55. $Q$ = *of everything still unexplained after
   every technical correction, the share that knowing a cell's relatives removes* — the talk number.
   **⭐ M2 (per-cell capture) is NEGATIVE on Subclone and Pre-TX** and strongly positive on the mice:
   they carry a ≥100-tape cell filter against ≥20, so where QC is strict a per-cell parameter is pure
   noise out of sample. Fig 3c's "the shelf is a QC choice" as a held-out *loss*.
   ⭐⭐ **Subclone is FLAT across the floor sweep** (+8.29/+8.31/+8.31 — its clones are all large,
   so the floor drops 21 of 38,636 cells): the control proving the rise elsewhere is about *which
   clones are retained*, not about the sweep itself. ⭐ **$k$ optimum is near 20 even at $n_C$ in the
   thousands** (Subclone +6.82/+8.31/+7.57 at $k$=5/20/50), so the fall at large $k$ has two distinct
   causes: $k/n_C$ on small clones, and reaching past the genuinely related cells on large ones.
   **⭐ Clone floor monotone 4/4** (Pre-TX +0.75→+4.59 from floor 20→100), confirming Justin's
   intuition — but **clone size makes the measurement cleaner, not the effect bigger**, and the rule
   is $k\lesssim0.2\,n_C$.
   **⇒ NEXT: the $\gamma$/$\beta_z$ skew test** (is the clone-tape coefficient right-skewed, as
   one-directional Dollo loss predicts, against a permutation null?), then the talk figures.

0-PREV2. **⭐⭐ 2026-09-11 — BOTH STATISTICS MEASURED ON ALL FIVE ARMS AT 100% COVERAGE; THE
   OWED LIST IS CLOSED AND THE SIMULATOR IS THE ONLY THING LEFT.**
   `2026-08_park-compatibility`, README "Session 10". **⚠⚠ Two defects in `73`/`74` found by reading
   the code before sizing the jobs.** (i) *A cell could be its own nearest relative*: self and
   cross-clone pairs were both $-\infty$, so they tied at $+\infty$ under `argsort` and were ordered
   by index — exposure at $k=50$ was **71.5% of Pre-TX cells**, and it does not cancel in obs − null.
   ⚑ The published Mouse2 c76 result was unaffected (0.0% exposure) and reproduces to $10^{-15}$.
   (ii) *443× of the arithmetic was discarded*: a pooled $n\times n$ matrix where only within-clone
   blocks are used. Both fixed by working **one clone at a time** — Pre-TX **68 min 53 s → 33 s**,
   and `--nsub` is gone, so coverage is 100% (Subclone 5%→100%, Mouse2 13%→100%).
   **⭐ Variogram, obs − null top bin:** Pre-TX **+0.2618** ($t=95.9$) · Subclone **+0.1406**
   ($t=31.8$, 127.0 M pairs) · Mouse2 +0.0419 · Mouse1 +0.0369 · Mouse3 +0.0190. Null flat in every
   bin of every arm; interior monotone 5/5. **Pre-TX is the strongest arm in the thread.**
   ⚠⚠ **CORRECTION: "Subclone is a cliff, not a gradient" was a 5%-subsample artefact** — at full
   coverage it is strictly monotone 13/13. `--nsub` *reshapes* the curve, it does not merely add
   noise. ⚑ Only the SLOPE is a finding: the pair-weighted mean is pinned at the all-pairs reference
   by $\sum_c r_{cz}=0$, so the negative low bins are the complement of the positive high ones.
   **⭐ Prediction task, best $k$, nats/cell (obs − null).** Largest clone: Subclone c2 **+11.29** ·
   Pre-TX c7 **+10.64** · Mouse2 c76 +4.04 · Mouse3 c110 +3.33 · Mouse1 c36 +1.73. Pooled: Subclone
   **+9.26** · Pre-TX **+2.86** · Mouse2 +2.65 · Mouse1 +0.88 · Mouse3 +0.62.
   **⭐⭐ The best $k$ INVERTS with clone size** — large clones rise with $k$, small ones collapse
   (Pre-TX pooled +2.86 → +0.18 from $k=5$ to $k=20$), same $\gamma$ constraint. **Quote $k$ relative
   to $n_C$, never absolutely.**
   ⭐ *Subclone pooled, 38,636 cells: where the model says 7%, the tape is missing in **3%** of cells
   whose relatives all have it (n=268,916) and **90%** of those whose relatives all lack it
   (n=5,715)* — the same read-off as Mouse2 c76 but on fat counts, retiring the thin-cell caveat.
   **✅ SD floor swept: it does not matter** (a 20× change moves the curve <0.5%), because the floor
   binds where the model is *right*. **✅ Fig 4e built** (`75`, three standalone PNGs).
   **⇒ NEXT: THE SIMULATOR.** Calibration target, per arm, through the identical scripts: the
   variogram's slope and convexity (not its level), the top-bin value (+0.019 to +0.262, a **14×
   spread across arms**), and the nats-by-$k$ curve including its inversion.

0-PREV. **⭐⭐ 2026-09-10 — OFF EVENT COUNTING. Two statistics that need no events.**
   `2026-08_park-compatibility`, README "Structured dropout without events". Justin's call, recorded:
   events are hard to justify and probably useless for simulator calibration, prevalence (1.3–25%) is
   a weak headline, and the collapse is a greedy heuristic that fragments 4.2× on Subclone.
   ⚠⚠ **The confound both scripts defeat:** relatedness is read from the edit data and dropout decides
   which edits are readable ⇒ **DISJOINT TAPE SPLIT** (relatedness from half A, dropout from half B).
   **`73` variogram** — similarity of Pearson-residual dropout profiles against lineage relatedness.
   **Monotone in 4/4 arms run**, null flat in every bin: Subclone −0.014→**+0.152** ($|t|$ to 39),
   Mouse1/2/3 −0.03→+0.034. Pre-TX pending.
   **`74` prediction task — the digestible one.** Predict a cell's dropout from its $k$ nearest
   relatives (found on half A, self excluded), one scalar $w$ on top of the fixed predictor
   $\eta=\alpha_c+\beta_z+\gamma_{Cz}$, evaluated on held-out cells. **Mouse2 c76: +3.04 nats per
   cell** (null −0.001). ⭐ *Where the model says 14%, the tape is missing in 9% of cells whose
   relatives all have it and **99%** of those whose relatives all lack it.*
   ⚠ Quote observed − null only ($\gamma$ forces $\sum_c r_{cz}=0$, biasing clone-mates by
   $-1/(n_C-1)$). ⚠ Quote the table, not the odds ratio. ⚠ Print entry counts — extreme cells are thin.
   ⚠⚠ **The effect is CONCENTRATED, not diffuse** — that reconciles every seemingly-inconsistent
   number in this thread. **Stratify, never average**, to convey magnitude.
   ⚠⚠ **Completeness is partly the attribution rule** ($\Lambda_{\rm hard}$ prefers complete clades):
   power-controlled, strong in 3 arms, equivocal in Mouse1/Mouse3.
   **⇒ NEXT: Pre-TX variogram, `74` on all arms, sweep the Pearson SD floor, then the SIMULATOR —
   calibrated by pushing simulated data through `73`/`74` and matching the curve.**
   **Direction 2 PARKED** (evidential, not modelling; effect on $q$ is a fraction of a percent).
0-NEWEST. **⭐⭐ THE TEST IS ANALYTIC AS OF 2026-09-10 — `67_exact_merge.py`, and the machinery is gone.**
   No calibration draws, no clade-size strata, no permutations required (kept only to validate).
   The $\gamma_{C,z}$ fit forces $\sum_{c\in C}(X_{cz}-\tilde p_{cz})=0$, so under the within-clone
   permutation a clade is a simple random sample **without replacement** from its clone's residuals,
   giving closed-form moments $\mathbb{E}_\pi[T]=m\bar r$,
   $\mathrm{Var}_\pi[T]=m\frac{n-m}{n-1}\sigma^2$ (population = the clone's cells in that BLOCK).
   Then $p=\max(p_{\rm normal},p_{\rm hypergeometric})$ and one Bonferroni cut
   $p\le\mathcal{E}/N_{\rm test}$, which bounds expected false **events** because events are unions
   of combos and the collapse only merges or drops.
   ⚠⚠ **The model variance $V$ is wrong by ~2×** (permutation SD of $z$ is 0.51–0.67, not 1); my
   Monte Carlo check simulated the model null with $\tilde p$ FIXED and so assumed its conclusion.
   Calling was unaffected; the "standard deviations" interpretation was not.
   ⚠⚠ **Edgeworth dropped before building** — at $z=6$ the correction is 37× the leading term.
   ⚑ Strata unnecessary ($\mathrm{Var}_\pi$ conditions on $m$ and $n$); testability structural
   ($\sigma^2>0$, $m<n$; 53.3% on Mouse3 — the rest have $T\equiv0$ by construction).
   **Mouse 3: 92 events / 2.97% at $\mathcal{E}=1$**; at matched stringency the three routes converge
   (128 / 127 / 132 events, 85 clone×tape pairs shared by all three). **9 s / 1.76 GB**, zero
   permutations ⇒ laptop-scale for future data. Launched on all five arms 2026-09-10.
   ⚠ Owed: normal tail validated only to $p\approx10^{-3}$ vs events at $10^{-9}$; Bonferroni over
   dependent tests over-corrects; the hypergeometric is 13× conservative.
0-NEW. **⭐⭐ THE DETECTOR WAS REBUILT (2026-09-09) — an event is defined on two axes.**
   `2026-08_park-compatibility`, README "The detector, rebuilt". **Detect** with the one-sided score
   test $z=(k-E)/\sqrt V$ — the score of a NESTED family $X_c\sim\mathrm{Bern}(\sigma(\eta_c+\delta))$
   in which $\delta=0$ *is* $H_0$; **attribute** with $\Lambda_{\rm hard}$ (maximised at the largest
   clade that is still complete — the Dollo rule, via
   $\Lambda_{\rm hard}=\Lambda_{\rm soft}-m\,\mathrm{KL}(\hat\pi\|1-\varepsilon)$); **report**
   $\hat\pi=k/m$ instead of assuming it. Mouse 3: **127 events / 3.19% of missing, $\hat\pi$ median
   1.000 against 0.376 predicted** — the first non-circular completeness claim — converging on the
   committed routes (128/3.09%, 132/3.15%).
   ⚠⚠ **$\Lambda_{\rm soft}$ is not a completeness-agnostic $\Lambda_{\rm hard}$**: non-nested against
   $H_0$, so it mixes elevation with the dispersion of $\tilde p$ inside the clade, and **72% of its
   Mouse3 calls had no rate elevation at all** ($z<0.5$). Not a calibration failure — the alternative was
   a mixture. Do not resurrect it as a detector.
   ⚠⚠ **Supersedes two recorded numbers**: partial fraction **29.9%** (not 8.9% — the hard route cannot
   see partials, so that was circular), and the **$\hat\pi$-by-depth gradient is gone**, so "graded
   losses are largely clade coarseness" (09-03, re-confirmed 09-07) was a detector artefact.
   ⚠ Calibration draws are **appended**, not held out — `62 --calib N` alone, one scan each.
   **Fig 5 a/b/c built on Mouse 3**; panel a is **faceted by clade size**, one threshold line per facet —
   that line is the rule and every event clears it. ⚠ Never pool the ceiling across clade sizes: the first
   version did and put 58 of 127 events below a line they were never tested against, because the decision
   variable is on neither axis. **⭐ NEXT: the four other arms** (`64_submit_percombo_soft.sh <arm>`).
0b. **⭐ Does the co-integrated symbol vanish when a tape is silenced?** (`2026-08_park-compatibility`,
   README "Two directions from the B = 1,000 run"). The first **orthogonal** test of the silencing
   mechanism — every result so far infers it from missingness, which is also what technical dropout
   looks like. pegRNA and tape share one cassette and pegRNAs act in **trans**, so silencing
   integration $z$ should remove symbol $s(z)$ from **every other tape** in those cells. ⚠ The
   $z\mapsto s(z)$ map is unknown (`TargetBC` 10-nt vs `NNNN` 4-nt, never linked), so **recover it
   rather than assume it** and let its structure be the evidence: five-way cross-arm concordance
   (166 TargetBCs verified identical across arms), near-injectivity at the collision rate predicted
   for 166 draws from 256 ($E=122$ distinct, 100–106 observed), $\mathrm{corr}(\beta_z,\xi_{s(z)})>0$.
   **⚠ Step 0 is a power calculation**: what fraction of (tape, site) slots is polymorphic within a
   clone / within a clade — tapes are ~4.5–5 of 6 saturated, so only post-loss insertions can show
   the depletion.
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

**⚠⚠ Do not list `lesliec` for work that needs neither a GPU nor a ~1 TB node.** Measured
2026-09-03: 300 one-core 8 G array tasks submitted `-p lesliec,cpu` landed **194 CPUs on `lesliec`,
76% of the lab's four private nodes**, held by one user — while `cpu` had ~9,700 CPUs idle. The four
`lesliec` nodes are the lab's *only* GPU nodes, so pure-CPU work parked there can block a labmate's
A100 job on CPUs with the GPUs free. Resubmitting the same array `-p cpu` returned `lesliec` to
233/256 idle and cost nothing in queue time. **The `-p lesliec,cpu` pairing in the examples above is
for jobs that are few, fat, or GPU-bound; a wide array of small tasks goes to `cpu` alone.**

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
