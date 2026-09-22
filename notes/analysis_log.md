# Analysis log

Running record of analyses. One section per analysis, newest last.
Cross-reference the protocol section in `sciphy_notes.md` and the directory
under `analyses/`.

---


## 2026-08-28 — `2026-08_park-compatibility`, Step 0

Park data arrived (Justin's copy via Jihye Park), flattened out of a doubled path into
`data/cancer_metastasis/`. Six CSVs, 627 MB.

Ran §D.4b "Step 0, which precedes both" on all five edit tables
(`analyses/2026-08_park-compatibility/src/01_symbol_composition.py`).

**Structure confirms the design recorded in §1629:** 166 tapes × 6 sites = 996 columns, byte-identical
tape-barcode set across all five tables, 99,451 cells (Initial 37,810 · Mouse1 12,232 ·
Mouse2 6,899 · Mouse3 2,904 · Subclone 39,606).

**⚠⚠ Main result — $q=\sum_i \xi_i^2\approx0.0170$, not the $\approx0.004$ assumed throughout.**
Stable at 0.0169–0.0174 across all five independently sequenced tables. The old figure assumed a
flat distribution over $M=256$; only 167–244 symbols are observed and the distribution is skewed,
so the **effective alphabet is ~58**. This inverts the standing comparison to Mulberry — Park is
marginally worse than their $q=1/64$ "high diversity" regime, not better, and not "near
homoplasy-free". Correction inserted at §D.4b and flagged at §1681.

⇒ **The §D.4b homoplasy null must be recomputed at $q=0.0170$ before the predicted-vs-measured gap
is interpreted.** The residual framing survives; only the subtracted number changes.

Two handling decisions carried forward: missingness is 43–62% (encoded `None`), and ~2.3% of
observed entries are malformed non-NNNNGGA strings (1,939 distinct, largest class 6-mers) that are
neither edits nor missing — must not be silently coerced to `None`.

**Next:** recompute the null at the measured $q$; join `clonalbc_percell_hamming1_corrected.csv`
to work per clone (§D.4b step 1); then the compatibility check proper.

## 2026-08-28 (cont.) — $m$ estimation, first pass

Built the project env at `/data1/choij10/justin/envs/pando` (system `miniforge3` module, not
normantm). Extracted full $\xi$ vectors; $q$ reproduces Step 0 exactly.

**Identifiability of $m$.** $\mathbb{E}[s\mid m]=\sum_i(1-(1-\xi_i)^m)$ and its variance derived and
**verified against Monte Carlo** (agree to 3–4 s.f.). Inverting the ±1 sd band: $m$ recoverable to
better than 2× up to $m\approx1{,}199$, no upper bound above $m\approx5{,}997$.

**⚑⚑ The real alphabet is ~100 symbols, not the design's 256.** 100 symbols carry 99.95% of edits;
the other ~97 carry <0.05% and are artifacts. They nearly broke the estimator: being rare they keep
$s$ climbing, making the inversion *look* conditioned out to $m\sim10^6$ when it is not. $q$ is
untouched by them — again the robust statistic.

**⚑⚑ Main result — homoplasy is rare and concentrated.** Over 1,567,321 prefix nodes, 95.8% sit
below $m^\ast=10.8$; mean recurrences rise from 0.01 (clade 3) to 50.5 (clade >500); and the
**top 1% of nodes carry 73% of all homoplasy**. By §D.4b's step-5 rule that is the favourable
world — concentrated conflict is deletable. The metastasis mice are near homoplasy-free at their
actual clade sizes; the subclone arm (median clade 933 at level 0) carries almost all of it.

⇒ This substantially *improves* the outlook for the perfect-phylogeny route relative to the two
corrections earlier today. The $q$ correction and the one-collision mechanism both raised the
homoplasy prediction, but Park's clades are small enough that the birthday threshold is rarely
reached.

**⚠ Two data facts to resolve.** `ClonalBC` has 3,294 barcodes, median 7 cells, max 27,537 — not
the "~75 clones × ~74 cells" of §1629. And clone-barcode dropout is 45% in Mouse1 vs 1.7% in
Subclone.

**Next:** the set-dependent (Poissonised MLE) estimator — $s$ is not sufficient for $m$; then
character-set construction and the compatibility check proper.

## 2026-08-28 (cont.) — second pass: the Poissonised MLE

$s$ is not sufficient for $m$: the exact likelihood $P(A\mid m)=\sum_{B\subseteq A}(-1)^{s-|B|}W_B^m$
depends on the *mass* of the observed set, not just its size, but has $2^s$ terms. Poissonising the
draw count makes the per-symbol counts independent and the likelihood exact and factorised:
$\log L(m)=\sum_{i\in A}\log(1-e^{-m\xi_i})-m(1-W_A)$, strictly concave, unique root, and it returns
$\hat m=\infty$ exactly when $A$ is the whole alphabet. Derivation documented in full in the header
of `src/08_poisson_mle.py`.

**Validated against the true fixed-$m$ model.** Both estimators near-unbiased (Poissonisation costs
≤1.6%); the MLE's gain is variance — IQR 26.4% vs 35.0% at $m=500$, 17.6% vs 21.0% at $m=200$, and
no gain below $m\approx10$. Exactly the predicted regime split.

**On the data:** the two agree to <1% for $s\le25$ (99.5% of nodes) and diverge to a ratio of 0.799
at $s\ge81$ (411 nodes), where $\hat m\sim1/(1-W_A)$ makes the observed mass decisive.

**⚠ Correction to a check I nearly reported.** The first model check compared observed $W_A$ to
$\mathbb{E}[W_A\mid m]$ and appeared to show misspecification at large $s$. It was confounded: the
data are selected on $s$, and at fixed $m$ a large $s$ arises by hitting unusually many rare symbols,
depressing $W_A$. With the null conditioned on the realised $s$ by simulation, observed $W_A$ sits
inside the 90% band at every $s$ — **no misspecification detected.**

**Conclusions unchanged:** 96.0% of nodes below $m^\ast$ (was 95.8%), top 1% of nodes hold 65% of
homoplasy (was 73%), total recurrences 574,959 (was 629,428). MLE numbers are the ones to use.

## 2026-08-31 — clone structure resolved; junk values promoted to symbols

**§D.4c — the "~75 clones × ~74 cells" was a misreading.** Read the paper. It is Metient's
*migration* subset: Mouse 1 only, clones with ≥10 cells present in ≥2 organs (93 → 76 → 75 after
their collision screen, 5,551 cells). Tree reconstruction has no such threshold. We reproduce their
93 and 76 exactly, so our ClonalBC handling *is* their pipeline. Also recovered from the methods:
their cell filter (≥100 recovered tapes for Initial/Subclone, ≥20 for the mice) is **not** applied in
the delivered tables; applying it drops 2.1% of cells.

**Character sets built** (`src/10_build_characters.py`) with an explicit determination mask $D$
alongside each clade $S$ — $D$ is what stops dropout on one tape being read as absence on that tape.
Trailing `None` inside a *recovered* tape is a determined absence, not missing data, which is why the
determined fraction only falls 0.758 → 0.725 from $L=1$ to $L=6$ despite 24% tape absence.

**⚠⚠ CORRECTION — junk values are heritable characters.** The earlier decision to truncate prefixes
at any non-NNNNGGA value was wrong, and the reasoning ("`CACGGA` occurs 47,287 times, it would
manufacture false clades") assumed its own conclusion. Two tests, no tree needed: preceding-prefix
concordance and clone restriction against a per-table calibrated null. Junk is **more**
clone-restricted than real symbols on both Mouse1 (+0.444 vs +0.222) and Subclone (+0.436 vs +0.295),
and the library-artifact confound is excluded by the clone/sample ratio. Mechanistically expected —
pegRNA scaffold read-through physically writes scaffold sequence into the tape, irreversibly.

⇒ **Alphabet rule is now frequency alone: ≥1,000 occurrences in the table, no sequence-pattern
requirement.** Pipeline rerun. $q$ 0.0172 → 0.0166, alphabet 94–106 → 97–129, characters +4.1%,
$m$ recoverable to 2× out to 2,757 (was 1,199). **Concentration unchanged: top 1% of nodes hold
64.9% of homoplasy (was 65.0%).**

**Next:** step 4, the compatibility test itself. 19.9 B cross-tape pairs total, but 8 clones hold 95%
of them — brute-force the other 2,518 clones, tape-pair contingency for the rest.

## 2026-08-31 (cont.) — step 4 begun; session crashed mid-run

**Ran and kept:** compatibility check on **Mouse3, all 91 clones** (`results/compatibility_Mouse3.json`).
Missing-as-absent 63.56%, missing-excluded **94.66%**, **spread +31.09 points** — the last row of
§D.4b's decision rule, i.e. **dropout is the binding constraint (row A6)**. Spread widens with clone
size (+12 at 3–4 cells → +42 at 11–20).

**Two engines, both verified.** `13_compatibility.py` cross-tabs label vectors, checked against
explicit brute-force set operations on real clones under both conventions. `14_compat_sparse.py`
recasts it as two sparse products ($MM^\top$, $MN^\top$), agrees exactly, runs 5× faster and reaches
clones 13 could not.

**Lost to the crash:** the conflict-graph degree distribution — computed and printed, never
persisted. It indicated conflict is *not* concentrated (~58% of characters conflicting, top 10%
holding ~46% of edges), which is the opposite of the homoplasy result and would put us in §D.4b's
"spread thin" world. Script 14 now writes `conflict_degrees_{table}.npz`; **Mouse3 needs a rerun to
recover this.**

**Never ran:** Mouse2, Mouse1, Initial, Subclone.

⚠ **Process change:** the run was executing in the VS Code tunnel and took the tunnel and the session
down with it. All analyses now go through `scripts/submit.sh` as Slurm batch jobs. Recorded in
CLAUDE.md.

## 2026-08-31 (cont.) — step 4 run; the C route turns out to be invalid

**Engine.** `13_compatibility.py` (cross-tab, verified against brute force on real clones) and
`14_compat_sparse.py` (two sparse products, verified against 13). Script 14's first version built the
determination matrix as characters × cells, making $MN^\top$ near-dense ($C\times C$); Mouse1 and
Mouse2 both died `OUT_OF_MEMORY` at 128 G. Fixed by exploiting the fact — already used for storage —
that **$D$ depends only on (tape, level)**, so at most 166×6 = 996 distinct masks exist per clone.
Mouse3 retest: identical results at 16 G instead of 64 G.

**Mouse3, all 91 clones:** missing-as-absent **63.56%**, missing-excluded **94.66%**, spread
**+31.09 points**. Conflict graph **not** concentrated (58% of characters conflict; top 10% hold
45.9% of edges) — unlike homoplasy, whose top 1% held 65%. ⇒ most observed conflict is *not*
homoplasy; it is dropout.

**⚠⚠ Main methodological finding — see §D.4d.** Computing $C$ from the missing-excluded conflict-free
characters gave $C/(n-1)=2.107$, impossible. Missing-excluded compatibility is **pair-specific**
(laminar on $D_1\cap D_2$, which differs per pair) and does not compose into the single laminar family
a tree needs. Verified on real pairs. ⇒ $C$ must be computed against the missing-as-absent graph;
**94.66% is latent agreement, 63.56% is what a skeleton can be built from.**

**Process.** Everything now runs through `scripts/submit.sh` as Slurm jobs after the previous session
crashed the tunnel mid-run. Two sizing lessons recorded in CLAUDE.md: a mid-run `MaxRSS` is worthless
when a job processes work in ascending size order (it measures the cheapest units), and the analytic
$MM^\top$ bound missed the real driver entirely.

**Next session:** save both conventions' degrees, recompute $C$ on missing-as-absent, finish Mouse1 /
Mouse2 / Initial, and decide whether the well-determined-subset variant (§D.4d) is worth trying.

## 2026-09-01 — compatibility finished on two tables; strategic pivot; figure programme begun

**Compatibility measured.** Mouse3 (91 clones): as-absent **63.56%**, excluded **94.66%**, spread
**+31.09 pts**. Initial (1,780 clones): **80.60%** / **91.91%**, spread **+11.31 pts**. Initial's
clones are tiny (median 6), which is why its spread is smaller — the spread-scales-with-clone-size
relationship holds across tables, not just within Mouse3. Mouse1/Mouse2 never completed (5 OOMs at
128/96/48/32 G, each on its single largest clone).

**$C$ measured** on Mouse3 against the missing-as-absent graph (the only constructible one, §D.4d):
$C_{\rm greedy}=1{,}003/1{,}765 = 0.568$ via conflict graph + min-degree greedy, $783$ via near-linear
insertion with 500 restarts, $689$ without restarts.

**Near-linear skeleton built** (`17_skeleton_linear.py`): deduplicate clades, then insert
largest-first maintaining `owner[cell]`, where the whole test is "all cells of $S$ share one owner".
**0.1 s and 153 MB against script 14's 9 min and 24 GB** — the cost we fought all session was an
artifact of computing the full pairwise matrix, which a skeleton never needs.

**Ground-truth test passed, and inverted a conclusion.** Pooled subclone colonies, skeleton built
blind to barcodes. After implementing Park's collision screen (3.77% of cells pruned), median clade
purity **1.000**. The large apparent false clades were *correct*: `ClonalBC` over-splits colonies,
and collapsing the four high-similarity group pairs gives **exactly 8 colonies — the paper's number**.
Group-consensus similarity is bimodal (off-diagonal 0.048 vs 2.18–2.54).

**⚑⚑ STRATEGIC PIVOT.** The skeleton is a step sideways. (i) Its output is arbitrary — 45% spread
across three reasonable heuristics on identical data. (ii) Only **9 of 2,547 clones exceed 1,000
cells**, 7 of them in Subclone, so 382 of 384 mouse clones are already within likelihood range;
median clone size is 4. (iii) It is a hard combinatorial device on soft data, and Felsenstein pruning
handles dropout natively — a cell's ~120 observed tapes place it while the missing tape marginalises.
⇒ **The diagnostics are the deliverable**; the skeleton's residual role is soft decomposition for the
~9 large clones.

**Figure programme begun** (PI presentations, motivating the likelihood route). Fig 1 (recorder/$q$)
and Fig 2 (homoplasy, `simple` + `mle` versions) done. Fig 3 needs redesign — its panel b argues
against a strawman and is unreadable. **Figs 4–6 all depend on one simulator that does not exist
yet**, and the homoplasy null it would provide is still the single biggest gap in the argument.

**Next:** redesign Fig 3 around dropout's magnitude and structure; then build the simulator
(birth–death tree, measured $\lambda$/$\xi$, $N=6$, $k=166$, **dropout as a switchable layer** so
homoplasy-only / dropout-only / both can be separated) for Fig 4.

## 2026-09-01 (cont.) — Fig 3 rebuilt around dropout structure

Fig 3 redesigned from "how we handle `None`" to a measurement of the assay. Panel a (the old
decomposition) dropped at Justin's call — nobody disputes that trailing `None` is biology.
`src/23` caches the (cell × tape) recovery/depth matrices once; 24–31 read the cache.

**a — dropout is not a coin flip.** $R_c\sim\mathrm{Bin}(166,p)$ gives sd 6.30 tapes; observed 29.9.
**VIF 22.6**, $\rho_{\rm cell}=0.131$. $R_c$ is **bimodal**: a mode near 120 plus a shelf from the
QC cut at 20 up to ~100 holding 38% of Mouse1 cells.

**b — the tape axis is bigger, and knowable.** Per-tape rates span 0.006–0.962, sd 0.244 against a
null sd of 0.0044 (55×). **$\rho_{\rm tape}=0.250 > \rho_{\rm cell}=0.131$.** ⚠ The VIFs (3,052 vs
22.6) are *not* comparable — $\mathrm{VIF}=1+(m-1)\rho$ and $m$ is 166 vs 12,232. Counting noise is
0.033% of the between-tape variance, so the spread is real analytically; replicate libraries agree
at $r=0.9968$ (confirmatory, so not plotted).

**⚑ $\rho$ verified three ways** (`src/29`): variance identity 0.13084/0.24948, brute-force Pearson
over 3 M sampled entry pairs 0.13064/0.24727, ANOVA between/total 0.130/0.249. Spearman ≡ Pearson
exactly (ranking a 0/1 variable is affine). $P(\text{miss}\mid\text{miss})$ measured 0.4729/0.5445
against $(1-p)+\rho p$ = 0.4733/0.5452.

**c — the shelf is a QC choice.** Cut all five arms at ≥100 and $\rho_{\rm cell}$ collapses to
0.012–0.024 everywhere, mice most homogeneous. Batch ruled out first: $\rho_{\rm sample}=0.022$
pooled over 12 harvest samples, a seventh of $\rho_{\rm cell}$.

**d/e — informative on the tape axis only.** Per tape $\rho=+0.34$ (+0.28 controlled), decile means
2.92 → 4.95 sites. Per cell $\rho=+0.05$, 4.88 → 4.98 — flat. **⚑ The tape coupling is a threshold,
not a gradient:** 25 tapes below $\hat\beta_t=0.3$ average 3.58 sites, the other 141 average 4.86
with $\rho=+0.12$ among them ⇒ **the informative part of dropout sits in a removable 15%.**

⚑ Also measured: $P(\text{no ClonalBC}\mid R_c)$ falls 64%→17% in Mouse1 and is flat in Pre-TX, so
the shelf and barcode loss are one phenomenon — the analysis population is selected on capture twice.

**Palette.** Five arms now coloured from the `dataviz` reference palette, slots 1,3,4,5,7; orange
reserved for the null. The shipped JS validator cannot run (cluster node v10.24), so it was ported
to Python with identical thresholds and CVD matrices.

**Next:** assemble a–e into the 2×3 grid; decide whether the A9 lineage test (do related cells lose
the same tapes?) is a sixth panel or its own figure with the simulator's null.

**Session close.** Fig 3 delivered as **five standalone panels (a–e)**, not a composed grid —
Justin's call, recorded in the analysis `CLAUDE.md` so no assembly script gets built later.
`22_fig3_dropout.py` is kept as the record of the superseded design. **Next session: plan the A9
lineage figure before building it.**

## 2026-09-02 — row A9 measured: tape loss is heritable, and heritable below the clone

**The question** (Mulberry & Stadler's stated reason for punting on dropout, §1c.2): are DNA tapes
"simultaneously lost for groups of related cells"? Unmeasured in the literature.

**Statistic.** $r_{cf}=Y_{cf}-\hat p_{cf}$ with $\hat p=\sigma(\alpha_c+\beta_z)$ — the *surprise*,
so a clone of uniformly poor cells leaves no residual and only **tape-specific** structure survives.
Intraclass correlation of those residuals within a group, computed without enumerating pairs via
$\sum_{c\neq c'}r_cr_{c'}=(\sum r)^2-\sum r^2$: $O(nk)$ for all 166 tapes, at clone and sample level
so between-clone comes free. **T0 is a calibration curve**, not a number — the same statistic on
"tape reached site $L$", whose marginals overlap the missingness marginals.

**T1, clone level** (excess over a within-sample permutation null): Subclone **+0.263**, Mouse1
**+0.166**, Mouse3 **+0.163**, Initial **+0.161**, Mouse2 +0.009 — **heritable in every arm**, at
30–130% of a marginal-matched heritable character. ⚠ Mouse2's small value is a *power artefact*: one
clone holds 3,387 of its 5,382 cells, so the permuted group is nearly the real clone.

**T2, below the clone** — subclades defined by the depth-$d$ prefix of an anchor tape, missingness
measured on all *other* tapes. The within-clone null is **analytic** ($m(m-1)/n(n-1)$ times the
clone's own pair sum) and was **verified against explicit permutations: mean ratio 1.0000**.
⚑⚑ **Monotone rising in 5/5 arms** (Initial +0.063→+0.155, Mouse3 +0.019→+0.059, Mouse1
+0.007→+0.020, Subclone-screened +0.007→+0.012, Mouse2 +0.0004→+0.0038). A flat batch effect predicts
zero at every depth. ⇒ closer relatives agree more about which tapes they lost.

⚑ **The depth test bounded the Subclone colony≡well confound** rather than arguing about it: Subclone
has the largest clone-level signal but the flattest gradient (sub-clone structure ~5% of its
clone-level number), while Initial's sub-clone excess nearly equals its clone-level one.

⚑ **Unexpected, and recorded:** the collision screen barely moves the clone-level test but **halves**
the subclade one. A misassigned cell carries a foreign prefix *and* foreign dropout, so it
*manufactures* subclade agreement instead of diluting it — the screen matters more for the fine test,
and in the opposite direction.

**T3, capture quality itself** ($\rho_{\rm clone}$ on $R_c$, the birth–death sampling question):
Initial **+0.153**, Mouse1 **+0.090**, Subclone +0.038, Mouse3 +0.032, Mouse2 +0.004, against
between-clone values ≤+0.021. Related cells share overall capture. ⚠ Not yet callable as biological
— needs a control for co-encapsulation / sub-lane structure.

**⚠ Two corrections to earlier claims in this project.**
1. Last session I said this figure's honest null is the simulator. **Wrong** — permutation nulls are
   sufficient for detection; the simulator is only needed to quantify the *consequence*. That
   unblocked the whole analysis.
2. $\rho_{\rm between}$ is contaminated as a comparator: with no per-sample term in the fit, the
   sample-level mean residual is not zero and gets squared over group sums of thousands of cells.
   The permutation null carries the same offset and is the sound contrast.

**⇒ Modelling consequence.** Heritability *below* the clone is what makes this a Dollo character on
the tree rather than a per-clone nuisance, and the remedy is a per-tape irreversible loss process —
**one absorbing state in the pruning recursion, a constant factor**. On this evidence **A9 does not
force SBI**, contrary to the §1c.3 table; only the uniform-sampling half (clustered *cell* loss) does.

**Housekeeping.** A `--verify` run submitted with `--depths 2` overwrote the full Mouse2 depth JSON;
re-run, numbers never lost. **Next: dig further, then plan the Fig 4 panels.**

## 2026-09-03 — Dollo vs propensity, and a null that had to be replaced

The depth gradient established that dropout concordance follows the topology, but **two mechanisms
predict a gradient**: a single irreversible loss (⇒ one absorbing state per tape, cheap) or a
heritable per-lineage *rate* at that locus (⇒ a latent field over the tree, a real SBI case).
`src/36_dollo_test.py` discriminates them from the distribution of $a$ = fraction of a subclade
missing a tape, restricted to subclades whose **own clone still carries it**.

**⚠⚠ Correction made mid-analysis.** The first pass compared $a$ to Bernoulli draws from
$\sigma(\alpha_c+\beta_z)$ and produced enormous enrichments (thousands-fold). That null assumes the
additive-logit fit is right; **any** cell×tape interaction makes $a$ U-shaped against it, and the
enrichment at $a\approx0$ was itself evidence the fit is too smooth. Replaced with a **within-clone
permutation of whole cell rows** — preserves each cell's full profile and each tape's marginal,
destroys only cell↔subclade alignment. The permutation null runs ~1.7× the Bernoulli one, so much of
the first-pass signal was cell×tape structure, not lineage.

**Against the honest null.** Complete absence below the clone: **1.4–12×** in Mouse1/2/3 and
Pre-TX (4–25× on tapes recovered in ≥50% of cells), but **0.9–2.2× in Subclone**, at or below the
null at shallow depths — consistent with its clone-level signal being largely colony-≡-well batch.

**⚑⚑ The shape is the evidence, 5/5 arms.** Mass moves *out of* "nearly all missing" and *into*
"exactly all missing": $[0.85,0.90)$ 0.70–0.91×, $[0.90,0.95)$ 0.76–1.00×, $[0.95,1]$ **1.34–1.50×**,
middle 0.92–1.01×. A graded propensity fills the near-complete bins; a discrete irreversible loss
empties them into the complete bin. And the sharper the criterion the larger the excess ($a=1$
exactly gives 2–4× where $a\ge0.95$ gives 1.3–1.5×) — a discrete state, not a rate.

⇒ **Dropout is a mixture**: a technical component described by $\alpha_c\beta_z$, plus a heritable
discrete loss. Cheap remedy confirmed — a per-tape irreversible loss process, one absorbing state.
Effect sizes are modest over the honest null, so it should be modelled, not treated as dominant.

⚠ **Open, and it bears on Fig 3b:** if part of per-tape missingness is heritable locus loss rather
than per-observation failure, $\beta_z$ is not a purely technical constant. Cross-library $r=0.997$
does not settle it — all arms descend from one engineered line, so an ancestral loss is shared by
construction.

**Also this session:** the three-centring decomposition for Mouse1/Mouse3/Pre-TX, generated last
session and never read, shows $\rho_{\rm within}$ *rising* when per-cell capture is removed
(0.134→0.168, 0.140→0.185, 0.150→0.161) — the heritable signal is partly masked by capture, not
produced by it, even in Pre-TX where capture is strongly clone-clustered.

## 2026-09-02 (cont.) — per-tape nulls, and the pooled estimator turns out to be the problem

**Per-tape permutation nulls** (`src/37`, B = 5,000 / 2,000). Each tape gets its own null
distribution, z-score and permutation p-value, with BH-FDR and `depth>=6` as a known-heritable
control. ✅ **Null calibration passes cleanly**: zero significantly-*negative* tapes across five arms
and two feature families — ten chances for a bad null to manufacture anti-concordance, none taken.

⚠⚠ **But significance saturates and must not be the headline.** 77–100% of tapes significant, median
$z$ from +8 to +422, because $\rho$ pools over millions of within-clone pairs. The p-value answers
"is there *any* detectable excess", not "is it large". The per-tape null is also right-skewed
(median skew +0.14 to +1.58), so $z$ is an **ordering statistic only**. **Effect size is the
quantity.**

**⚑⚑ What the per-tape nulls did buy.** Dropout heritability is **2–2.5× more concentrated across
tapes than edit-depth heritability measured identically** (top decile 43–58% vs 16–24%, in 3 of 5
arms). Edit depth accumulates along lineages on every tape alike; dropout heritability piles into a
minority of loci — what discrete losses at particular integration sites predict and a diffuse
technical effect does not. Also links to Fig 3: Spearman(recovery, per-tape excess) is negative in
every arm (−0.11 to −0.70) — poorly recovered tapes are the heritably lost ones.

**⚑⚑ The pooled $\hat\rho$ was hiding the signal** (`src/38`, prompted by Justin asking what summing
over clones does). $\hat\rho$ is a ratio of sums, so clones enter weighted by $n_C(n_C-1)$: the
largest clone holds **98.9% of Mouse2's weight**, 84.0% of Mouse1's, 1.6% of Pre-TX's.

| arm | pooled | drop largest | **equal weight per clone** | clones ≥5 cells, % positive |
|---|---|---|---|---|
| Mouse2 | +0.009 | +0.326 | **+0.246** | 78, 99% |
| Mouse1 | +0.166 | +0.226 | **+0.222** | 137, 100% |
| Mouse3 | +0.163 | +0.148 | **+0.237** | 67, 100% |
| Pre-TX | +0.161 | +0.163 | **+0.214** | 1,685, 100% |
| Subclone | +0.263 | +0.342 | **+0.353** | 13, 100% |

⇒ **Mouse2 was never a weak arm** — its dominant clone has near-zero excess and swamped the pooled
estimate. Equal-weighted the five arms collapse to **+0.21 to +0.35**. The replication claim becomes
**~1,980 individual clones, essentially all positive**, not "5/5 arms". ⇒ **Use equal-clone weighting
from here on.**

**Panels a and b built** (`src/39`): `fig4a_heritable.png` (excess vs marginal frequency against the
`depth≥L` control curve the $\sigma^2\le p(1-p)$ bound makes necessary; missingness at 31–126% of the
matched control) and `fig4b_perclone.png` (per-clone distributions, each arm's largest clone ringed —
Mouse2's 3,387-cell clone sits *on the null line* while its 77 smaller clones sit above it).

**⭐ Direction for next session (Justin's call).** Panels **c and d may carry the figure alone**. The
effect as hypothesised is a *tree* property, not a clonal one, so the most visible and defensible
demonstration is that **it appears WITHIN large clones** — at subclade resolution inside a single
clone, where a clone-level effect cannot reach. The decisive case is already identified: Mouse2's
3,387-cell clone shows near-zero *clone-level* excess yet holds 98.9% of the arm's weight; if the
structure is real it should be there internally. Same for Subclone's 10,996-cell clone. Dissect
these next.

## 2026-09-03 — the event catalogue, and the mechanism verified from the paper

**✅ Construct and readout checked** (`refs/metastasis_lineage_recording.pdf`, Methods p39–42).
Park's recorder is one cassette, `PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC` — **pegRNA and tape
co-integrated**, as in the other lineage configurations — and the tape is read from **cDNA** (10x 3′
v4, *"TAPE cDNA co-amplified"*). ⇒ **a tape is recovered only if its integration is transcribed**, so
silencing removes the tape from the readout *and* the co-integrated pegRNA's symbol from the writing
pool, both heritably. The mouse paper already blames low per-tape recovery on *"epigenetic silencing
of a subset of circTAPE-encoding integrations"*; what is new is that it is lineage-resolved.
⇒ **Re-explains Fig 3b** ($\beta_z$ is the integration's expression level, not primer efficiency) and
**Fig 3d** (one chromatin state suppressing transcription and editing at the same locus).

**The statistic** (`src/40`): a log-likelihood ratio for a Dollo loss on a clade's stem,
$\Lambda=\sum_{\rm miss}\log\frac{1-\varepsilon}{\tilde p}+\sum_{\rm pres}\log\frac{\varepsilon}{1-\tilde p}$.
Capture-independence is built in (a good cell missing a reliable tape earns +2.98 nats, a bad cell
+0.09) and all-or-none is enforced (a present cell costs up to −4.55, and *more* for a good cell).

**⚠⚠ The first version was useless and the fix is instructive.** Without a per-(clone,tape) margin, a
**clone-wide** loss made every subset of that clone look spectacular in real and permuted data alike;
the FDR sat at **64–73% at every threshold**. Adding $\gamma_{C,z}$ fitted to the clone's own margin
— which fixes *how many* cells lack the tape but not *which* — dropped the FDR to **≈0%**.

**Results, 5/5 arms at FDR ≈ 0:** 73–8,220 events per arm, **inside missing rate 1.00 everywhere**
against 0.11–0.41 expected, and **event cells never worse captured** ($R_c$ 113–131 vs 114–131).

**⚑⚑ Detection power is the whole story.** Stratified by clone size (`src/41`): **zero events in any
clone under 20 cells, in any arm** — $\gamma$ is fitted from that clone's cells and absorbs
everything. Share of a band's missing entries inside a called event climbs monotonically to
**6.5–19% in clones ≥200 cells**. Pooled percentages were dilution artefacts.

**Clone-wide layer** (`src/42`): **7.1–32.9% of all missing entries**, 155–4,763 (clone,tape) losses
per arm, inside rate 0.94–1.00 vs 0.17–0.34 expected. ⚑ Founder-39 was monoclonal and mRFP-sorted, so
all 166 integrations were active at cloning; ClonalBC came day −11 and the bottleneck day 3 ⇒
**clone-wide = silencing before the bottleneck, sub-clone = after. One mechanism, two epochs.**

**Soft variant** (`src/43`), fitting the clade's own rate $\hat\pi=k/m$ instead of pinning it at
$1-\varepsilon$: $\hat\pi$ median **0.989–1.000** against 0.26–0.51 expected, **~half exceed 0.99**,
and genuinely partial events ($\hat\pi<0.90$) are **8–32%**. ⇒ the hard test is strict but not badly
so; an absorbing state covers most of it, with a real minority wanting a lineage-varying rate.

**⚠ Three corrections recorded.** (1) "Fraction of missingness explained" is retired — $\gamma$
matches clone×tape totals, so excess inside an event is balanced by deficit elsewhere; the quantity
measures *concentration*, not an additive share, and the layers never partition a total. (2) Wilks
does not apply to $\Lambda_{\rm soft}$: $H_0$ is non-nested in $H_1^{\rm soft}$, so permutation is
the only calibration. (3) "Fitting $\varepsilon$" was wrong — any estimate from called events is
selected on having few present cells; a **sensitivity sweep** is honest, and shows $\varepsilon$ is
not load-bearing (10× range moves the hard/soft overlap by <7 points).

**Next:** Pre-TX and Subclone soft runs still in flight; then panels c (one worked example inside a
large clone) and d (the catalogue).

**Soft variant, final two arms.** Pre-TX 67,628 events, $\hat\pi$ median 1.000 (expected 0.254),
67.3% ≥0.99, 13.8% partial. **Subclone 516,549 events, $\hat\pi$ median 0.857, only 26.4% ≥0.99 and
57.2% partial** — inverting the other four arms.

⚠ **Subclone's inversion is most likely clade resolution, not biology.** $\hat\pi<0.90$ pools genuine
graded silencing with **a complete loss on a smaller clade than the one tested**; prefix clades stop
at depth 4, and Subclone's clones run 200–11,000 cells, so they are coarse relative to the true tree
(and its expected rate, 0.115, is the lowest of any arm). ⇒ **8–32% is an upper bound on genuine
partial silencing in the mice.**

⭐ **Cheap decisive test, not run:** rebuild the prefix cache at `MAX_D = 6` (tapes have six sites)
and re-run. If partial events resolve into complete events at depth 5–6, the graded component is an
artefact; if they persist at maximum resolution it is real — absorbing state vs lineage-varying rate,
decided event by event.

## 2026-09-03 (cont.) — the MAX_D = 6 test: partly artefact, partly real

Rebuilt prefix codes to the full six sites (`src/44`, new filename, alignment with the depth-4 cache
asserted) and reran the soft catalogue at depths 1–6 (`src/43 --maxd 6`), stratifying $\hat\pi$ by
clade depth. If "partial" events are complete losses scored on too-coarse clades, $\hat\pi$ must rise
with depth.

**It does, and not all the way.** Median $\hat\pi$ reaches **1.000 at depth 6 in all three mice**, and
the partial fraction falls: Mouse3 **12.7% → 5.3%**, Mouse2 37.5% → 26.0%, Mouse1 22.9% → 19.8%.

⇒ **Clade resolution explains much of the partial fraction but not all.** Mouse3 is essentially all
artefact; Mouse1/Mouse2 keep ~20–26% partial at **maximum recorder resolution**, so a genuine graded
component survives. ⇒ absorbing state for the majority, lineage-varying rate for a bounded minority.

⚠ Determination cost: (cell,tape) pairs determined fall 76–78% at depth 1 to 21–24% at depth 6
(Pre-TX 3.9%), so deep clades are scarcer and smaller — the trend is read within that.

**Still in flight at time of writing:** Pre-TX and Subclone `--maxd 6`
(`logs/43_soft_events_{Initial,Subclone}-1135651[89].out` → `results/soft_events_{arm}_d6.json`).
**Subclone is the decisive one** — its depth-4 partial fraction was 57.2%, and if that collapses at
depth 6 the inversion was clade resolution throughout.

**Sizing note:** the depth-4 soft runs took 35 s – 11 min against 6–16 h requests. Right-sized the
depth-6 reruns to 2 h.

⚠ **Process note.** Two commits this session lost their `notes/` and root-`CLAUDE.md` edits because
the shell was left in `analyses/2026-08_park-compatibility/` — a bare `notes/analysis_log.md` does
not exist there, and a bare `CLAUDE.md` resolves to the *analysis* file, so a `str.replace` keyed on
root-file text silently no-opped. **`cd` to the repo root before touching notes, and assert on every
replacement.**

**⚠⚠ MAX_D=6, final two arms — and they correct the reading above.** Pre-TX 15.1%→6.8% partial;
Subclone **85.2%→37.7%**, median $\hat\pi$ 0.507→0.947. The residual partial fraction at depth 6 is
**monotone in maximum clone size**: Mouse3 (210 cells) 5.3%, Pre-TX (127) 6.8%, Mouse1 (1,607) 19.8%,
Mouse2 (3,387) 26.0%, Subclone (10,996) 37.7%.

⇒ That is the recorder running out of resolution before the tree does — six sites resolve at most six
levels, which suffices only for small clones. Subclone is still steeply falling at the limit, with
clades a median of 12 cells inside clones of up to 10,996. **So the earlier "~20–26% genuine graded
residue" was wrong; extrapolating to arms where six levels suffice puts it at ~5–7%.**
⇒ **A per-tape absorbing state is the right and largely sufficient extension**; the
lineage-varying-rate case is weaker than the depth-4 numbers implied. Inference from a cross-arm
trend, not a direct measurement — a simulator with known ground truth would settle it.

**⚠⚠ Correction, and the calibration that forced it** (`src/45`). I claimed the soft statistic's
fitted-parameter advantage is "$\approx0.5$ nats, by Wilks". Both parts wrong. Wilks does not apply
($H_0$, per-cell $\tilde p_c$, is not nested in $H_1$, one shared $\pi$), and measuring
$\Lambda_{\rm soft}$ on **all 5,541,360 candidate (clade,tape) pairs** of Mouse3 with no threshold
gives a **null mean of −1.74**, not +0.5: the fitted-constant model is typically *worse*, because one
degree of freedom does not buy back the per-cell structure it discards. ⇒ **a positive
$\Lambda_{\rm soft}$ must overcome a deficit first; it is not a fitting artefact.** Enrichment over
the null: 1.8× at 2 nats, 5.6× at 10, 8.8× at 16, and 3,705 observed vs **0** null at 25.
A +1 sd fluctuation in a 20-cell clade buys 0.45 nats; the 10-nat threshold needs +4.6 sd.

**Next task specified: B = 1,000 permutations for proper p and q values.** Full spec — including the
correct statement of the permutation, measured per-scan costs, and the parallelisation design — is in
`analyses/2026-08_park-compatibility/CLAUDE.md` under "NEXT TASK".
*(Pointer updated 2026-09-03: that heading is now "B = 1,000 permutations → proper p and q values —
LAUNCHED 2026-09-03". The spec itself was followed; the one departure is the monotonisation
direction — see session 4 below.)*

⚠ **The permutation, stated correctly** (Justin asked, and the natural phrasing is wrong): we do not
shuffle *clone* labels for the sub-clone test. The rule is **permute the label being tested, blocked
by the label above it** — subclade membership within clone for the sub-clone tests, clone membership
within sample for the clone-level ones. Whole cell rows move; labels stay at their positions. The
within-clone version preserves each clone's own rate for every tape, which is why clone-wide losses
are invisible to the event catalogue and needed script `42` instead. The full scan is redone per
permutation, which is what makes $B=1{,}000$ expensive: measured per-scan costs give ~2 min for `42`
but ~11 h (`40`) and ~52 h (`43`, Subclone) serially, so `40`/`43` need ~20-way splitting over
permutation parts with a merge step.

---

## 2026-09-03, session 4 — B = 1,000 launched across everything; two new directions

**Implemented and submitted.** `40`/`43` gained `--permpart i/N`; `46_perm_merge.py`,
`47_submit_B1000.sh` and `48_collect_B1000.sh` are new. 15 permutation configurations (5 arms ×
{`40` hard sub-clone, `43` soft depth-4, `43` soft depth-6}) at 20 array tasks each, plus 5 direct
`42` runs. Design and the count-vector argument are in the analysis README under "B = 1,000
permutations — launched".

Three choices worth carrying forward:

- **Store counts, not candidates.** A part writes only `#{candidates ≥ t}` on a fixed 600-point
  grid. Subclone's soft scan yields ~876,000 candidates *per permutation*; 1,000 such lists would be
  ~7 GB and discarded. 600 integers is everything the FDR curve, the $q$'s and the global $p$ need.
- **Seed permutation $b$ from $(\text{SEED},b)$**, not from a stream advanced $b$ times. Slice
  boundaries then cannot change what any permutation *is*; parts merge in any order and a duplicate
  is detectable by index. `46` asserts the pooled set is complete and disjoint and names what is
  missing — a silently absent part would shrink the null and inflate every $q$'s denominator.
- **Keep the whole $(B\times G)$ count matrix**, not a mean: $q$ needs the mean null count, the
  global $p$ needs the distribution.

⚠ **Correction to the spec written last session.** It said the $q$ curve is monotonised "by a running
minimum from the top". It is the other way: an event at $\Lambda$ can be reported by any rejection
region $\{\Lambda'\ge t\}$ containing it, i.e. any $t\le\Lambda$, so $q=\min_{t\le\Lambda}q_{\rm raw}(t)$
— a running minimum **up** the grid, which is BH's direction and is automatically non-increasing in
$\Lambda$. Downward would invert that *and* drag the $q_{\rm raw}=1$ convention (assigned wherever the
observed count is zero, i.e. above the largest observed $\Lambda$) across the whole informative range.

**First results — the five `42` clone-wide runs, complete.** $p = 0.000999 = 1/1001$ in every arm at
both the chosen threshold and the 10-nat scan floor: **not one of 1,000 permutations came close.**
Loss counts unchanged from $B=200$ (Mouse 1 725, Mouse 2 410, Mouse 3 398). $B$ bought resolution on
$p$, not a different answer — which is the honest thing to report.

### ⚠⚠ Being a good lab citizen — a real error, caught by Justin

Submitted with the project's habitual `-p lesliec,cpu`, the 300 one-core 8 G tasks landed **194 CPUs
on `lesliec`, 76% of the lab's four private nodes**, one user holding them, while the general `cpu`
partition had ~9,700 CPUs idle. Worse than the raw share: those four nodes are the lab's **only** GPU
nodes, so pure-CPU work parked there can block a labmate's A100 job on CPUs with the GPUs free.
Cancelled and resubmitted `-p cpu`: `lesliec` back to 233/256 idle, this work 309 of 14,264 on `cpu`
(2.2%), every task `RUNNING` immediately — **no queue time lost by behaving well.** Recorded as a
standing rule in the root `CLAUDE.md`: the `lesliec,cpu` pairing is for jobs that are few, fat, or
GPU-bound; a wide array of small tasks goes to `cpu` alone.

### ⚠⚠ $B$ is not the detection floor — three knobs were being conflated

| knob | controls | what $B=1{,}000$ does |
|---|---|---|
| $B$ | resolution of $p$ | fixes it, $0.17\to0.001$ |
| FDR **threshold** | which candidates are *called* | lowers it where a 3-permutation null had pushed it up |
| scan **floor** `--lam` | which candidates are *collected at all* | **nothing — hard-fixed at 10 nats** |

The floor is a cutoff inside `scan()`; this run's grid starts at exactly 10 and can never see
beneath it. Mouse 3's soft threshold of 16.3 nats came from **three** permutations and should fall
toward 10 — but weaker events need a re-run at `--lam 4`, which *subsumes* the floor-10 run. Script
`45` already shows the terrain: 4.6× enrichment at $\ge4$ nats on Mouse 3. Deliberately not done by
restarting the live jobs — read this run's FDR curve at the floor first, then size it.

### ⚠⚠ Per-combo permutation p-values — asked, and they are the wrong tool at this floor

Justin proposed scoring each (clade, tape) combo against its own 1,000 permuted values. Well-defined
(prefix codes are never permuted, so a clade *slot* persists with its size fixed) and cheap, given an
exceedance **counter** per combo (~22 MB) rather than $B$ values per combo (~44 GB for Mouse 3 alone).
**But at the 10-nat floor it is strictly weaker than the pooled count.** Both spend $B=1{,}000$;
pooling buys $B\times N_{\rm combos}\approx10^{9}$ null draws (tail resolution $\sim10^{-7}$) against
1,000 ($10^{-3}$). Measured: **eight** permutations of Mouse 3's hard scan produced **zero** null
candidates above 10 nats across every combo, so every real event censors at $1/1001$ — the statistic
saturates exactly where the signal is strongest. Pooling's price is assuming $\Lambda$'s null is
exchangeable *between* combos, and the repair for that is stratification, not per-combo scoring.

⚠ **A second floor on $p$ that no $B$ can lift.** An $m$-cell clade in an $n_C$-cell clone has only
$\binom{n_C}{m}$ realisable permuted compositions, so $p\ge1/\binom{n_C}{m}$: a 4-cell clade in a
6-cell clone **cannot reach $p<0.05$ at any $B$** (15 compositions). ⇒ a structural mechanism for
script `41`'s "not one event in any clone under 20 cells", beyond $\gamma_{C,z}$ absorbing
everything. Worth saying aloud — it makes the power limit information-theoretic rather than a
choice of threshold.

⇒ **Folded into the `--lam 4` follow-up:** lower floor + **null stratified by clade size** (the real
repair for "weaker but real": a 500-cell clade has a far heavier $\Lambda$ tail than a 5-cell one and
one global FDR fits neither; per-stratum count vectors cost ~6× storage, i.e. nothing) + per-combo
counters, which stop being censored at $\Lambda\ge4$.

⚠ Terminology: the run already yields a per-event **$q$** (the FDR of the rejection region containing
it). It does not yield a per-event **$p$**. Different questions; $q$ is the more useful for a catalogue.

### ⭐ New direction — does the co-integrated symbol vanish when a tape is silenced?

Justin's idea, and it is the sharpest test available: every result so far infers silencing from
**missingness**, which is what technical dropout also looks like. pegRNA and TAPE share one cassette
and pegRNAs act in ***trans*** (§0), so silencing integration $z$ should remove symbol $s(z)$ from
**every other tape** in those cells — a channel transcript capture cannot reach. It measures what §0
currently asserts under "nasty coupling" and what row **A9** is built on.

⚠ Not a deduction: pegRNA is Pol III (U6), tape/mRFP Pol II (EF1α). Locus heterochromatin should take
both, but nothing forces it — which is exactly why it is worth measuring.

⚠ **The map $z\mapsto s(z)$ is unknown** (`TargetBC` 10-nt vs `NNNN` 4-nt, never linked) — the same
gap §"Open empirical question: is there a *cis*-preference" records for the Typewriter data.
**⇒ Recover the map instead of assuming it, and let its structure be the evidence.** Measured today:
all **166 TargetBCs are identical across all five arms**, so a recovered map has five independent
replicates; 166 draws from $4^4=256$ predicts **122** distinct symbols against **100–106** observed,
so near-injectivity at a *predicted* collision rate is checkable, and collisions predict partial
rather than complete drops. Also measured: $\xi$ is **smooth over a 570× range** (0.00008–0.0456),
not quantised by copy number — no cheap shortcut, but a large dynamic range for
$\mathrm{corr}(\beta_z,\xi_{s(z)})>0$ to live in.

**⚠ Step 0 is a power calculation and precedes the examples.** The tape is append-only, so symbols
written before the silencing stay; only post-loss insertions can show depletion. Tapes here are
~4.5–5 of 6 saturated — precisely the regime where most content is ancestral. Measure the fraction of
(tape, site) slots **polymorphic within a clone**, and within a called clade. That number decides
whether any of this works. It is also the same discipline §"cis-preference" already demands under
*phylogenetic non-independence*: count each edit once at the branch where it first appears, never
once per cell.

Then: 3–5 hand-inspected examples → screen the **clone-wide** layer first (strongest: 4,763 losses
over 1,188 Pre-TX clones, and the loss predates the clone founder so more of the clone's editing
postdates it) → the validation ladder → only then the sub-clone version.

**In flight at time of writing:** 80 array tasks, the depth-4 and depth-6 soft scans on Pre-TX and
Subclone (`perm43d{4,6}_{Initial,Subclone}`); `perm_collect` (11389784) pends on them and will pool
the parts and re-run each observed scan once with `--nullfile`.

### 2026-09-04 — B = 1,000 landed, all 15 configurations

`perm_collect` pooled 300 parts into 15 nulls of exactly 1,000 permutations, re-ran every observed
scan, 17 min, all completeness assertions passed. **$p=1/1001$ in all twenty runs** (hard, soft d4,
soft d6, clone-wide × five arms) at both the threshold and the scan floor.

**⭐ Quote the null maximum, not the p-value.** The most extreme of 1,000 random labellings produced
**111** candidates in Subclone where the real labelling produced **279,973** — 2,522×. Across arms
the ratio to the null *max* is 118× (Mouse 2) to 2,522× (Subclone). This answers Justin's own
framing ("where does the real labelling sit among the 1,000?") directly, makes the null explicit,
and is immune to "$1/1001$ is just the resolution floor".

⚠ The null is **right-skewed by 6–100×** (max/mean) — the same skew script `37` found for the
per-tape nulls. $B=3$ could never have seen that tail, and it is why the max rather than the mean is
the honest summary.

**Every hard-catalogue event clears $q\le1.9\times10^{-4}$**; the catalogue is significant as a
whole, with no weak tail. **Cell × tape entries inside called blocks**: Subclone 266,688 (19.09% of
all missing), Mouse 2 21,277 (6.43%), Pre-TX 22,256 (1.70%), Mouse 1 14,770 (3.64%), Mouse 3 2,916
(2.28%) — pooled shares are dilution artefacts (script `41`: nothing called below 20 cells).

⚠ **One genuinely marginal result to disclose rather than bury:** the clone-wide layer on Mouse 2,
FDR 4.36%, max $q$ **0.044**. Everything else is orders of magnitude clear.

#### ⚠⚠ CORRECTION to yesterday's prediction

I wrote that Mouse 3's soft threshold of 16.3 nats "came from three permutations and should fall
toward 10 once the null is properly estimated". **It did not** — 16.2 at depth 4, 15.5 at depth 6,
and every other configuration stayed exactly at the 10-nat floor where it already was.

⇒ **The three-permutation null was already unbiased in the mean; what it could not see was the
tail.** $B=1{,}000$ bought resolution on $p$ ($<0.17\to<0.001$) and the *shape* of the null, not any
new events and not one moved threshold. Stated plainly because the honest accounting matters: ~250
core-hours bought calibration and rhetoric.
⇒ It also sharpens the `--lam 4` case — the threshold now demonstrably sits **on** the scan floor in
14 of 15 configurations, so the floor is the only place more events can come from.

### ⚠⚠ 2026-09-04 — audit: how events are actually called (Justin asked; it matters)

Two thresholds, and **the arbitrary one is doing the work**:

- **scan floor `LAM0 = 10` nats** — hard-coded, fixed a priori. `scan()` only collects pairs with
  `L >= LAM0`; nothing below exists downstream.
- **calling threshold `LAM`** — adaptive, smallest grid point with monotonised FDR $\le5\%$.

**The adaptive step is degenerate in 13 of 15 configurations.** FDR at the floor is 0.002–0.019% for
the hard catalogue — **262× to 2,511× below the 5% target** — so `LAM` snaps to the floor. It binds
only for Mouse 3 soft (FDR at floor 10.03% d4 / 9.22% d6 ⇒ 16.2 / 15.5).
⇒ For the hard catalogue the operative threshold **is** the arbitrary 10 nats; the stated FDR rule
never bites.

**Consequence, and why the "how widespread" claim was not yet safe:** counts are set by the constant
and the *ranking between arms is not stable under it.* Pre-TX 1,783 events at $\Lambda\ge10$ but
**110** at $\ge20$ — 2nd place to 4th — because its events pile against the floor (78% within 5
nats, median 12.0, max 35, clades median 9 cells) while Subclone runs to $\Lambda=1{,}326$. The
dropout share moves likewise: Subclone 19.09% → 13.62% → 5.49% at $\ge10/20/100$; Pre-TX 1.70% →
0.20%. **⇒ report a curve or two thresholds, never a single number.**

⚑ **But the floor is conservatism, not rigour.** At FDR 0.002–0.019% the $\Lambda=10$–12 events are
almost certainly real — we are hundreds of times inside the cliff, so the counts are **loose lower
bounds**. This converts the `--lam 4` run from "fishing" into "applying the criterion we already
claim". Further structural lower-bound reasons: `MIN_CLADE=4`; clades exist only where an anchor tape
resolves them; $\gamma_{C,z}$ absorbs clone-wide by construction; nothing called under 20-cell clones.

**Checked and clean:** dedup keeps the highest-$\Lambda$ clade among overlapping ones for a given
(clone, tape), so **no cell × tape entry is double-counted**; and kept events are all-or-none —
`inside_rate` $\ge0.90$ at the 10th percentile in all five arms, median **exactly 1.000**, against an
expected 0.11–0.41.

⚠ **The 5% target is itself a choice not to inherit.** With $10^5$ candidates, 5% FDR admits ~14,000
false events in Subclone. For an inspectable catalogue an **absolute** criterion (expected false
events $\le10$) is more defensible than a rate. Decide this when setting the floor-4 threshold.

⇒ **Threshold-independent claims that should carry the talk:** the null comparison (118–2,522× beyond
the most extreme of 1,000 permutations), $q\le1.9\times10^{-4}$ for every hard-catalogue event, and
$\hat\pi$ median 1.000 against expected 0.11–0.41.

### ⭐ 2026-09-07 — the low-floor run landed, and it corrected me twice

`lf_collect` (40 min): 15 configurations, 300 parts, $B=1{,}000$ each, $p=1/1001$ throughout, all
merges complete. Then `52_collect_budget.sh` re-reported all 15 against the stored nulls under an
absolute criterion. **The criterion, not the floor, was doing almost all the work.**

**⚠⚠ The 5% target is unusable at a low floor, for a reason that did not exist at floor 10.** The
FDR is computed on *candidates*; the reported quantity is *events*, after the overlap collapse. At
floor 10 that gap was harmless (Mouse 1: 2.15 expected false candidates vs 320 events). At 5% it is
fatal — 5,515 expected false candidates vs 2,083 events — and the count-vector design cannot dedup
permuted sets to close it. **A candidate-level rate cannot bound an event-level error unless the
rate is tiny.** ⇒ added `--budget X`: per-stratum threshold = smallest $\Lambda$ with expected null
count $\le X$. An absolute budget bounds false *events* too, since dedup only reduces them.

**The defensible catalogue** (floor 2, stratified, expected false $\le2$/stratum): Subclone 6,597
events / 20.07% of missing; Pre-TX 2,019 / 2.46%; Mouse 2 452 / 7.02%; Mouse 1 424 / 4.34%;
Mouse 3 132 / 3.15%. Against floor 10 that is **1.1–1.8×, and 0.8× for Subclone** — not the 3–14×
the 5% target advertised. $\le11$ expected false candidates per arm; every event $q\le2\times10^{-3}$.

**⚠⚠ Correction 1 to 2026-09-04.** I wrote that the floor was "conservatism, not rigour" and the
counts "loose lower bounds". Direction right, **magnitude badly wrong**: I measured the distance to
the false-positive cliff in *FDR* units (262–2,511× below target) when the relevant units are
*nats*. The defensible per-stratum thresholds land at **5.9–12.4 nats** — essentially where the
arbitrary floor sat. The cliff is 0–4 nats below 10. ⇒ the floor-10 catalogue was accidentally
near-right; what this bought is a **justified** threshold rather than an inherited one, plus a
modest gain, plus knowing where the cliff is.

**⚠⚠ Correction 2, within the same session.** I computed each stratum's maximum achievable
$\Lambda=m\log((1-\varepsilon)/\bar p)$ with the **arm-median** $\bar p$, declared the 4–5 cell
stratum "DEAD in all five arms", and said so. **The budget catalogue calls 18–635 events there in
every arm.** $\tilde p$ is per (cell, tape), not per arm, and the events that clear are exactly
those where recovery was expected: called 4–5 events sit at $\tilde p$ 0.05–0.13 against arm medians
0.12–0.42, lifting their ceiling above the threshold. ⇒ the right statement is **not** "small clades
are undetectable" but "**a small clade is detectable only on a tape that should have been there**" —
capture-independence doing exactly the work it was built for. Small clades stay heavily suppressed
(34% of Mouse 1's candidates, 12% of its events) but the stratum is live.

**⚑ The `MAX_D=6` conclusion survives; my doubt about it was itself the 5% artefact.** I flagged
that the partial fraction rose at the lower floor (Mouse 3 d6 8.1%→18.0%) and that "the graded
appearance is very largely clade coarseness" needed re-reading. At the budget threshold: Mouse 1
15.9% (was 20.0%), Mouse 2 22.0% (29.6%), Mouse 3 8.9% (8.1%), Pre-TX 15.6% (13.7%), Subclone 46.3%
(50.3%) — unchanged, and *lower* in three arms. A partial loss scores lower $\Lambda$ by
construction, so **any** threshold drop inflates the partial fraction with no change in biology.
⚠ Always quote the partial fraction with its threshold.

**Still open:** `42` (clone-wide) is still unstratified, and clone sizes vary far more than clade
sizes — its Mouse 2 FDR of 4.36% is the last marginal number. An **event-level** FDR would need the
dedup run on permuted candidate sets, which count vectors cannot support; the absolute budget
sidesteps that rather than solving it.

## ⭐ 2026-09-08/09 — the detector rebuilt: nested families, the score test, $\hat\pi$ made non-circular

Justin asked to settle the "dropout is heritable" figure, and the thing never settled was **what a
silencing event is**. Two claims had been running together: **A**, dropout is non-exchangeable within
clones (established three ways, needs no event definition), and **B**, the losses are discrete and
complete — Dollo characters. Only B needs a definition, and B is what the SciPhy extension rests on.

**First: $\hat\pi$ and `inside_rate` are the same number** ($k/m$), so the completeness bar bolted onto
the per-combo route on 09-07 already was a $\hat\pi$ threshold. Hard and soft differ in exactly one
respect — whether completeness is fused into the statistic or reported separately. Made exact by
$\Lambda_{\rm hard}=\Lambda_{\rm soft}-m\,\mathrm{KL}(\hat\pi\|1-\varepsilon)$ (verified $<5\times10^{-13}$,
reproduces `43`'s docstring examples). Corollary: $\Lambda_{\rm hard}$ is maximised at the **largest
clade that is still complete** — the Dollo rule — which is why it was always the right attributor.

**⚠⚠ CORRECTION, and the Mouse 3 run caught it for 25 min of compute.** I proposed detecting on
$\Lambda_{\rm soft}$ ("$\Lambda_{\rm hard}$ without the completeness requirement") and reporting
$\hat\pi$. Wrong. It returned 340 events, median $\hat\pi$ 0.203 — of which **226 were clades of 50–199
cells with $k-E=+0.27$ cells and $z=0.08$**, i.e. no elevation whatsoever. Cause: $H_0$ (product of
*different* Bernoullis) is **not nested** in $H_1^{\rm soft}$ (product of *identical* ones), so they
differ in two respects at once — level and homogeneity — and
$\Lambda_{\rm soft}=m\mathrm{KL}(\hat\pi\|\bar e)+[\ell(\bar e\mathbf 1)-\ell(\tilde p)]$, the second term
being the dispersion of $\tilde p$ in the clade ($\mathbb{E}_{H_0}=-\sum_c\mathrm{KL}(\tilde p_c\|\bar e)$,
verified by MC). Real clades are more homogeneous in capture than random subsets of their clone, so that
term beats its own permutations with zero elevation. On the 226 artefacts the elevation term carried
**0.1%** of the median margin. ⚠ **Not a calibration failure** — the permutation null was exact; the
*alternative* was a mixture. $\Lambda_{\rm hard}$ escaped only by brute force ($-3.3m$ nats), i.e. right
by accident, and by the same term that blinds it to graded losses.

⚠⚠ **A second, earlier correction**: I first blamed the overlap collapse (63 ordered it by the detection
margin). Measured — collapsing by $\Lambda_{\rm hard}$ gives 341 events against 340. Coarse attribution
was ~15% of it, not the cause.

**The fix — a nested family.** $X_c\sim\mathrm{Bern}(\sigma(\eta_c+\delta))$ with
$\eta_c=\alpha_c+\beta_z+\gamma_{C,z}$: $\delta$ is one more additive term in the same logistic
regression, $\delta=0$ **is** $H_0$ in the interior, $\delta\to\infty$ is total loss, and the per-cell
heterogeneity is kept in both models so the misfit term cannot arise. Score and information at $\delta=0$
give $z=(k-E)/\sqrt V$ — simply (observed − expected)/(SD of the count), with $\mathbb{E}_{H_0}[S]=0$ and
$\mathrm{Var}_{H_0}[S]=V$ holding **exactly** for any $m$ and any heterogeneity. Heterogeneity enters the
scale, never the location; and by concavity of $p(1-p)$ a homogeneous clade gets the *largest* $V$, so the
artefact direction is suppressed.

⚠ Wald is unusable here (0.000 on a complete loss — complete separation), and $\hat\delta$ diverges for
complete losses, so $\hat\delta$ is the effect size for graded shifts and $\hat\pi$ for complete ones.
⚠ $\chi^2_1$ does not replace the permutation: a nominal 5% test rejects **9.07%** at $m=8$.
⚠ Score screens, exact LRT characterises the shortlist — purely cost (LRT is ~10× the scan; Pre-TX
2h → 20h per part).

**Calibration draws are now APPENDED, not held out.** Compare fresh permutations $B,B+1,\dots$ to the
stored $M=\max_{b<B}z_b$; same exchangeability, and NCAL stops being a launch-time decision (one scan
each, ~0.1% of the run). Validated on data generated under $H_0$ exactly: ties reproduced (0.0117 vs
0.0113, against 0.0053 continuous), quantiles match to 3 dp through the 99.99th percentile, LOO ratio
1.000–1.015, and **12 independent replicates** give ratios 0.998–1.218 with $t$ −0.19 to 1.75 — no
detectable anti-conservatism, bound ~20% at the smallest counts. ⚠ NCAL is the *only* source of precision
on the null count; HOLD=3 gave ±40% at the operating point.

**Mouse 3, all three detectors on the same 1,000 permutations.** $\Lambda_{\rm soft}$: 63,920 combos, of
which **45,916 are large clades with $z<0.5$** — 72% artefact. Score $z$: 19,844 combos, **zero** such.
$\Lambda_{\rm hard}$: 13,674, zero. ⇒ **127 events, 3.19% of all missing**, $\hat\pi$ median **1.000**
against $\bar e=0.376$ predicted, $\hat\pi\ge0.99$ in 60.6%, $z$ median 5.23 (min 2.26), LRT median 22.99,
$\hat\delta$ 60% diverged, smallest clade 4 cells, all beating every permutation. Converges on the
committed work from a different direction: 127/3.19% vs 128/3.09% (per-combo $\Lambda_{\rm hard}$, d4) vs
132/3.15% (floor-2 budget); 98 of 124 (clone,tape) pairs shared, 82 three-way.

⚠⚠ **Two numbers now superseded.** (i) The partial fraction is **29.9%**, not 8.9% — the hard route
cannot see partial losses, so its number was circular; this is the first honest one. (ii) **The
$\hat\pi$-by-depth gradient is gone** (1.000, 0.992, 1.000, 1.000, 0.992, 0.963 over depths 1–6). Under
$\Lambda_{\rm soft}$ it ran 0.155→1.000, so the "graded losses are largely clade coarseness" reading
recorded 09-03 and re-confirmed 09-07 **was an artefact of the detector**. Events are complete at every
depth; the partial tail sits in 6–9 cell clades where one present cell moves $\hat\pi$ a long way.

**Implementation checks.** `62` reproduces `53`'s observed $\Lambda_{\rm hard}$ **bitwise** (5,574,944
slots), before and after the score edit. Depth-6 prefix partitions **contain** depth-4 in all five arms,
so this is a superset of the committed run. Same seed ⇒ the two runs are paired permutation by
permutation. And $V\to0$ cannot blow up $z$: since $k\le m$, $z\le m(1-\bar e)/\sqrt V$, measured max
**0.06** across the 638,355 slots with $\bar e>0.999$ — the $\gamma$ margin working.

**Fig 5 built** (`66`, three standalone PNGs: the plane, the survival curves, the completeness histogram).
⚠⚠ **Panel a's orange "null ceiling" is misleading and must be replaced**: it is a maximum pooled over
clade sizes, so **58 of 127 called events sit below it**, every one clearing its *own* stratum threshold.
Root cause is structural — the decision variable (clade size) is on neither axis, so no curve in that
plane is the decision boundary. Fix by faceting on clade size, or by plotting margin − threshold(stratum).
Two other defects were caught on inspection and fixed: a filled contour drew a spurious *lower* null
boundary, and panel b was a per-bin histogram mislabelled as cumulative.

**Next:** the four other arms (same code; Pre-TX 24 G, Subclone 16 G), panel a refaceted, and `42`
(clone-wide) still unstratified.

## ⭐ 2026-09-10 — exact permutation moments, analytic p-values, and the machinery goes away

Justin pushed back: the 12 calibration draws are convoluted and hard to justify; why not p-values from
an LRT, accept that small clades have little power, and drop the bespoke apparatus? He was
substantially right, and chasing it exposed an error in the previous design.

**⚠⚠ CORRECTION — $V$ is the wrong scale, and my check of it was circular.** I asserted
$\mathrm{Var}[T]=V=\sum_{c\in S}\tilde p(1-\tilde p)$ *exactly* and verified it by Monte Carlo —
**simulating under the model null with $\tilde p$ FIXED**, i.e. assuming what was to be shown. Under
the **permutation** null with $\tilde p$ **fitted**, the SD of $z$ is 0.51–0.67, not 1 (measured from
the stored $s_1,s_2$). Calling was unaffected — the margin compares observed to permuted $z$ on the
same scale, so a constant factor cancels — but every reading of $z$ as "standard deviations" was
wrong. ⚠ The score is now written $T$; earlier notes used $S$ for both clade and score.

**The exact moments.** The $\gamma_{C,z}$ fit forces $\sum_{c\in C}(X-\tilde p)=0$ per (clone, tape)
(verified: max $2.1\times10^{-4}$), so under the within-clone permutation a clade is a simple random
sample **without replacement** from the clone's residuals. Hence $\mathbb{E}_\pi[T]=m\bar r$ and
$\mathrm{Var}_\pi[T]=m\frac{n-m}{n-1}\sigma^2$ — closed form, **no permutations required**. The factor
$(n-m)/(n-1)$ is exactly what $V$ lacked: at $m=n$ the clade is the whole population, $T\equiv0$, and
the variance correctly vanishes. ⚠ The population is block-conditional (the clone's cells reaching
that depth on that anchor), not the whole clone; the whole-clone form overstates the SD by 3–6%,
conservatively. **Verified against the stored permutations: predicted/empirical SD ratio 1.031–1.064
by stratum, 1.046 overall, correlation 0.93.**

**⚠⚠ CORRECTION 2 — Edgeworth was the wrong tool and was dropped before building.** Its skewness term
is a *central* correction: at $z=6$, skewness 1, the normal tail is $9.9\times10^{-10}$ and the
correction $3.6\times10^{-8}$ — 37× the leading term, so the expansion has diverged precisely where
decisions are made. Replaced by the **hypergeometric**, the exact permutation distribution in the
constant-$\alpha$ limit (Fisher's exact test on clade × missingness): discrete, skew-exact, no
expansion. Two p-values per combo, the conservative one used.

**Multiplicity collapses to one number.** $p\le\mathcal{E}/N_{\rm test}$ bounds expected false
*combos*; events are unions of combos and the collapse only merges or drops, so expected false
*events* $\le\mathcal{E}$ — **the candidate-vs-event mismatch that forced the absolute budget does not
arise**. ⚑ And the six clade-size strata are **gone**: $\mathrm{Var}_\pi$ depends explicitly on $m$
and $n$, so $p$ already conditions on clade size. The strata were a patch for an unstandardised
statistic.

**Testability is structural.** A combo is testable iff $\sigma^2>0$ and $m<n$ — 53.3% on Mouse3. The
other 47% are clades that *are* their whole block-conditional population, so $T\equiv0$ and no test
exists. That is the principled form of "do not call tiny clades": for many there is literally no test.

**Mouse 3.** At $\mathcal{E}=1$: 92 events, 2.97% of missing, $\hat\pi$ median 1.000 (expected 0.402),
$\ge0.99$ in 66.3%, partial 23.9%, smallest clade 6 cells. At $\mathcal{E}=12$ (matching the old
rule's ~12 expected false): 128 events, 3.79%. ⭐ **Stringency-matched the three routes converge —
128 (exact) / 127 (per-combo margin) / 132 (committed $\Lambda$ budget), with 85 (clone,tape) pairs
common to all three.** Validation on the 124,490 combos the permutations can resolve: median
$p_{\rm norm}/p_{\rm emp}=2.31$, $p_{\rm hyper}/p_{\rm emp}=13.0$, both conservative.
⚠ The hypergeometric's conservatism costs real power (normal alone clears 26,973 combos, the max
12,093); taking the max is a deliberate choice.

**Cost.** The merge runs in **9 s at 1.76 GB** against 10 min before, and needs **zero** permutations —
so the whole analysis becomes laptop-scale for future datasets. Scripts `67_exact_merge.py`,
`68_submit_exact.sh`; launched on all five arms.

⚠ Still owed: the normal tail is validated only to $p\approx10^{-3}$ while events sit near $10^{-9}$,
so beyond that it is extrapolation (the hypergeometric is the guard); Bonferroni over ~$10^6$
*dependent* tests over-corrects; and this is the **third** detector design in three days, so it is
held to the same standard — it reproduces the previous catalogues at matched stringency and both of
its approximations are checked against the stored permutations.

**Five arms landed (2026-09-10).** Cell × tape share of all missing: Subclone 25.48% · Mouse2 7.70% ·
Mouse1 6.74% · Mouse3 2.97% · Pre-TX 1.33%, against the committed $\Lambda$-budget catalogue's
20.07 / 7.02 / 4.34 / 3.15 / 2.46% — same ordering but for Mouse3/Pre-TX swapping. **(clone, tape)
membership agrees strongly** (1,240/1,382 Subclone, 774/817 Pre-TX, 154/163 Mouse2). $\hat\pi$ median
1.000 in 4/5 arms. Runtime 8 s – 12 min per arm.
⚠⚠ Two caveats recorded in full in the README: **Subclone's normal p is anti-conservative (0.46) with
no clean explanation** — the hypergeometric guard binds in 89.4% of its calls, so my earlier
suggestion to drop it for power is **withdrawn**; and **Subclone fragments at 4.2 events per (clone,
tape)**, reinforcing "quote cell × tape, never event counts".

## ⭐⭐ 2026-09-10 (cont.) — structured dropout without events: the variogram and the prediction task

Justin: events are hard to justify and probably useless for calibrating simulations; is there a
metric based on how similar related cells' dropout profiles are? Yes, and it is the best thing in the
project for that purpose. Scripts `73_dropout_variogram.py`, `74_profile_prediction.py`.

**The confound that would have invalidated it.** Lineage relatedness is read from the edit data and
dropout decides which edits are readable, so two cells that both lack tapes 1–50 look related
*because* their dropout matches. ⇒ **disjoint tape split**: relatedness from half A only, dropout
from half B only, repeated over random splits. Plus Pearson residuals (margins already removed) and
everything within clone.

**Pearson residual** $r^{*}=(X-\tilde p)/\sqrt{\tilde p(1-\tilde p)}$ — necessary because a
Bernoulli's variance depends on its mean, so a missing tape at $\tilde p=0.05$ is a 4.4 SD surprise
while at $\tilde p=0.5$ it is 1.0 SD, and raw products would be dominated by unreliable tapes.
⚠ SD floored so $|r^{*}|\le10.1$; **the floor binds on 16% of entries** and needs a sweep.

**Variogram result: monotone in 4/4 arms run.** obs − null from lowest to highest relatedness bin —
Subclone −0.0135→**+0.1523** ($|t|$ to 39), Mouse1 −0.0298→+0.0335, Mouse3 −0.0311→+0.0338, Mouse2
−0.0160→+0.0381. **The null is flat in every bin of every arm** (+0.0001 to +0.0009). Pre-TX still
running. ⚠ Pooling coverage varies a lot (Mouse3 97%, Pre-TX 88%, Mouse1 63%, Mouse2 13%,
Subclone 5%) because `--nsub` truncates the huge clones.

**⭐ The prediction task is the digestible one.** For cell $c$, take its $k$ nearest relatives *in the
same clone* from half A, **excluding itself**; form the neighbour signal over half-B tapes; fit ONE
scalar $w$ with the existing predictor as a fixed offset,
$\operatorname{logit}\Pr(X_{cz}=1)=\eta_{cz}+w\,u_{cz}$; evaluate on held-out cells.
Mouse2 c76: **+3.04 nats per cell** at $k=20$ against a null of −0.001 (sd 0.44 over 5 splits), and
the conditional table — holding the model's own prediction fixed — reads *where the model predicts
14%, the tape is missing in 9% of cells whose relatives all have it and **99%** of cells whose
relatives all lack it*.
⚠ Quote the table, not the odds ratio: $e^{w_f}=4.7$ badly understates a 9%→99% contrast because it
is a linear coefficient on a relationship that is flat then a cliff. ⚠ Extreme cells are thin
(35–345 entries); always print counts. ⚠ Quote observed − null only: the $\gamma$ fit forces
$\sum_c r_{cz}=0$, so a random clone-mate is negatively correlated with $c$ by $\approx-1/(n_C-1)$ —
~−5% in a 20-cell clone, comparable to the signal.

**⚠⚠ The organising fact: the effect is CONCENTRATED, not diffuse.** Most entries carry no lineage
signal, a minority carry an overwhelming one. That reconciles the variogram's 0.024, the event
route's $p=10^{-46}$, the 1.3–25% prevalence and the 9%→99% table — all correct, all different
questions. ⇒ **stratify, do not average**, when conveying magnitude.
⚠⚠ **Correction to my own claim**: aggregating over tapes does NOT inflate a correlation (it makes it
precise, ceiling $\sqrt{0.024}=0.155$); what aggregates is the likelihood gain.

**Also corrected today:** the completeness claim. $\hat\pi$ median 1.000 is partly the attribution
rule — the collapse is ordered by $\Lambda_{\rm hard}$, which is maximised at the largest COMPLETE
clade. Before the collapse, above-threshold combos have median $\hat\pi$ 0.88–0.98. Controlling for
detection power via slack $=(k-k_{\min})/(m-k_{\min})$: **strong in Pre-TX (0.87), Mouse2 (0.97) and
Subclone (0.93); equivocal in Mouse1 (0.57) and Mouse3 (0.56)**. Honest headline: *where there was
room not to be complete, losses are complete in three of five arms.*

**Decision recorded:** Direction 2 (the co-integrated symbol) is **parked, not discarded** — its value
is evidential, not modelling; the map is many-to-one ($\bar j\approx1.6$), depletion is partial
(6–47%), and the knock-on to $q$ is a fraction of a percent. The split-half result (87.3%, 33.9×
null) stands on its own if the biological interpretation ever needs defending.
**Also decided:** stop leading with event counts. Prevalence is a weak headline at 1.3–25% and the
collapse is a greedy heuristic that fragments 4.2× on Subclone; the excess curve (`72`) and now the
variogram/prediction pair replace it.

---

## 2026-09-11 — session 10: the owed list closed, and two defects found before running anything

Justin asked for a fuller report on the dropout-prediction results and for the outstanding
submissions. Reading `73`/`74` before sizing those jobs turned up two problems, so the sizing came
after the fixes.

**⚠⚠ Defect 1 — a cell could be its own nearest relative.** `74` marked self *and* cross-clone pairs
with $-\infty$; under `argsort(-rel)` those tie at $+\infty$ and are ordered by index, so the cell's
own row entered its top-$k$ whenever its clone held fewer than $k+1$ cells. Verified on a toy: a
4-cell clone at $k=5$ returns `[1 2 3 0 4]`. It is exactly the circularity step 2 exists to prevent,
and it does **not** cancel in observed − null, because the within-clone permutation replaces the self
term with a random clone-mate. Exposure at $k=50$: **Pre-TX 71.5%** of cells, Mouse3 29.7%, Mouse1
19.1%, Mouse2 9.3%, Subclone 0.1% — i.e. it would have bitten on the very next run in the queue.
⚑ The published Mouse2 c76 result had 0.0% exposure and reproduces to $10^{-15}$ after the fix.

**⚑ Defect 2 — 443× of the arithmetic was thrown away.** Both scripts built one $n\times n$ matrix
over all pooled cells and used only the within-clone blocks; the useful fraction is
$\sum_C n_C^2/n^2$, which is 0.22% on Pre-TX (549 clones of at most 127 cells). Blocking by clone
fixes defect 1 structurally *and* removes the waste: Pre-TX **68 min 53 s → 33 s** (per split
411 s → 2 s), and `--nsub` became unnecessary, so coverage went to 100% (Subclone 5%→100%, Mouse2
13%→100%, Mouse1 63%→100%).

**Validated, not asserted.** Both rewrites preserve the RNG draw order and pair ordering, so they
should reproduce the old numbers exactly — and do: `73` Mouse3 agrees to $7\times10^{-18}$ with
identical bin edges and counts; `73` Pre-TX reproduces all ten per-split values of the 68-minute run;
`74` Mouse2 c76 agrees to $9.8\times10^{-15}$ across 30 runs with a bit-identical conditional table.

**Results, all five arms.** Variogram top bin (obs − null): Pre-TX **+0.2618** ($t=95.9$) · Subclone
**+0.1406** ($t=31.8$, 127.0 M pairs) · Mouse2 +0.0419 · Mouse1 +0.0369 · Mouse3 +0.0190, null flat
in every bin of every arm. Prediction task, best $k$, nats/cell: largest clone Subclone c2 **+11.29**,
Pre-TX c7 **+10.64**, Mouse2 c76 +4.04, Mouse3 c110 +3.33, Mouse1 c36 +1.73; pooled Subclone
**+9.26**, Pre-TX **+2.86**, Mouse2 +2.65, Mouse1 +0.88, Mouse3 +0.62.

**⚠⚠ Correction to yesterday's reading.** "Subclone is not a gradient — it is a cliff" was an
artefact of the 5% subsample. At 100% coverage the curve is **strictly monotone across all 13 bins**
and the top bin moved +0.1523 → +0.1406. The lesson generalises: `--nsub` truncates the largest
clones hardest, so it **reshapes** the curve rather than merely adding noise.

**⭐ A new methodological finding: the best $k$ inverts with clone size.** Large clones rise with $k$
(Mouse2 c76 +1.48 → +3.41 → +4.04); small ones collapse (Pre-TX pooled +2.86 → +1.57 → +0.18). Same
$\gamma$ constraint as the null caveat: as $k\to n_C$ the neighbour mean tends to $-r_{cz}/(n_C-1)$,
an anti-signal. ⇒ **quote $k$ relative to $n_C$, never absolutely.**

**✅ The SD-floor sweep is answered: it does not matter** — a 20× change in the floor moves the top
bin by <0.5%. The reason is worth keeping: the floor binds on 3.8–23.4% of entries, but those are
entries where $\tilde p\approx0$ and the model is *right*, so $r^{*}\approx-0.001$ whatever the
floor; only the rare confident-and-wrong entries are rescaled, and they are too few to move a mean
over millions of pairs.

**✅ Fig 4e built** (`75_fig_variogram.py`, three standalone PNGs). Four label collisions were found
by rendering and inspecting; one annotation was deleted outright rather than repositioned, since
"minimal annotating text" is the standing preference and the point belongs in the README.

⇒ Owed items (1)–(4) are done, (4) superseded. **The simulator is the only one left**, and its
calibration target is now concrete: the variogram's slope and convexity, the top-bin value across
arms (+0.019 to +0.262, a 14× spread), and the nats-by-$k$ curve including its inversion — all from
the identical scripts, with no threshold, attribution or event definition anywhere in them.

---

## 2026-09-14 — session 12: the silencing evidence is called sufficient; figure programme closed

Justin's call, recorded: the four figures below are enough to demonstrate that heritable tape
silencing is real. Work now moves to the simulator.

**The four figures.** `fig6a` the ladder (held-out gain against parameters spent, per rung) ·
`fig6b` $Q$ (the lineage share of what the technical model left) · `fig6c` per-tape unanimity
(missing vs present, against the diagonal) · `fig6d` the companion (capture matched, relatedness
not). All five arms; 6c also at $k$=3.

**The numbers.** Ladder, floor 100, $k$=20, nats per held-out cell (obs − null): Subclone
**+8.31±0.09** ($Q$=**15.9%**), Pre-TX **+4.59±0.12** (7.1%), Mouse2 +2.47, Mouse3 +2.19,
Mouse1 +1.48 — **positive 5/5 at 16–90 se**, with one added parameter against 3,232–64,975 in the
technical model. Per-tape unanimity, $k$=5, common tapes, median fold relatives ÷ capture-matched:
Pre-TX **9.77× for missing against 1.22× for present**, Subclone **3.84× vs 1.15×**, mice ~1.2× vs
~1.03×.

**⚠⚠ Four corrections this thread paid for, none of them cosmetic.**
1. *The null must exclude self.* The $M_3$ score equation forces $\sum_{c\in C}r_{cz}=0$, so a
   self-excluding neighbour set carries $-1/(n_C-1)$ of the cell's own residual; `74`'s
   row-permutation null does not, so obs − null cannot cancel it. Pre-TX returned a **negative**
   lineage rung until fixed.
2. *A cell could be its own nearest relative* — self and cross-clone pairs both $-\infty$, tied
   under `argsort`, affecting 71.5% of Pre-TX cells at $k$=50.
3. *Stratifying on an estimate does not hold the truth fixed.* Simulating a world with **no** lineage
   structure, the "arbitrary clone-mates" curve still rose 14.9% → 28.2%, so the model-stratified
   figure's caption was wrong even though the underlying measurement was not.
4. *Capture matching must be by rank, not decile* — deciles inflated the mouse folds by up to 20%,
   and related sets really are worse captured (Mouse2 37.4 vs 54.5 A-tapes unmatched).

**⚑ And three framing rules.** Per tape, never pooled — the pooled unanimity figure was a mixture
dominated by near-dead tapes. Common-tape convention — all-$k$-missing needs a high rate and
all-$k$-present a low one, so the two panels otherwise sit on opposite ends of the range. And $k$=5
maximises the discriminator, with $k$=2 being the variogram in disguise.

**⇒ Next thread: the simulator**, to be discussed fresh. Calibration target per arm, through the
identical scripts: the ladder's 4-vector including the *sign* of the $M_2$ rung, the variogram's
slope and convexity, the nats-by-$k$ inversion, and the per-tape unanimity fold. These are **joint**
constraints on (tree, editing, dropout) rather than clean per-parameter ones — the variogram's
$x$-axis is the editing process — so tree and editing should be calibrated first against statistics
that do not involve dropout, and these used to falsify rather than to fit.

⚠ Parked with a known flaw: the $\beta_z$ skew test. Dispersion is hugely inflated ($z$ = 11–89 in
5/5 arms) but moment skew is not, because centring on the clone-weighted mean turns a one-directional
effect two-sided when there are few clones. Rebuild against a low quantile before concluding anything
about directionality.

---

## 2026-09-15 — session 13: scoping the simulator; SciPhy's Session 3 read; the tree model settled

**Scope corrected, twice, both times narrowing.**

⚠⚠ **Correction 1 — I proposed compute for a dead question.** I queued an Initial compatibility
rerun and a Subclone memory-pricing to serve fig 5 ("compatibility spread + homoplasy null").
Justin pushed back: max-compatible-set died in the 2026-09-01 pivot. Right. The pivot preserved
compatibility as a *diagnostic* ("the diagnostics are the deliverable"), but **that sentence is now
stale too** — the dropout diagnostic has been superseded by the ladder, $Q$, the variogram and
unanimity, all on five arms with proper nulls, and declared sufficient on 09-14. ⇒ **both runs
withdrawn.** Compatibility survives only as a free end-stage check: Mouse3 is complete on disk, so if
the finished simulator reproduces 63.56% / 94.66% / +31.09 through `src/14`, that costs nothing.

⚑ Found while checking: **Mouse3 needed no rerun** — `conflict_degrees_{as_absent,excluded}_Mouse3.npz`
have been on disk since Aug 31; the README note saying `src/14` saves only one convention is stale.
The "lost" concentration numbers, recomputed: **as-absent 91.3% of characters conflict, top 10% hold
27.1% of edges**; excluded 58.0% / 45.9% (reproducing the README exactly). The as-absent row
*strengthens* the finding — against homoplasy's top 1% carrying 65–73%, conflict is diffuse.
⚠ **Initial's compatibility run was budget-truncated and is biased**: 1,780 of 2,137 eligible clones,
**53% of characters**, and the largest clone processed had **41 cells against 127 eligible**.
`src/14` orders clones ascending and breaks on wall-clock. Its 80.60% / 91.91% / +11.31 is computed
on the small end, so the true spread is wider. Not worth fixing given the above — but do not quote it.

**⚠⚠ Correction 2 — the homoplasy null is narrower than I framed it.** For the mouse arms it is
largely unnecessary: script 06 already *measured* level-0 recurrences at 0.03–0.11 per node with
0.4–1.1% of nodes carrying any. It is load-bearing only for **Subclone** (median clade 933, 16.0
mean recurrences, 40.8% of nodes) — where the strongest dropout signals also live, so the two
processes are genuinely entangled — and for the **design sweep**, where $q$ is the sole channel from
$\xi$ to topology. And the currency should change: incompatibility is a perfect-phylogeny quantity,
while the likelihood route cares about *reconstruction accuracy*. ⇒ **fig 5 to be renamed and
rescoped** away from its skeleton-era title.

**Editing layer, measured (all five arms, from the cached `dropout_matrix_*.npz`).**
- ⚠ **Depth 0 is censored**: `10_build_characters.py` returns `ABSENT` when the prefix is empty, so
  "recovered but wholly unedited" is indistinguishable from dropout. Fits must be zero-truncated.
  At Initial's $\hat\Lambda\approx2.82$, $e^{-\Lambda}=5.9\%$ against a 25.07% `ABSENT` rate — **up
  to a quarter of Initial's apparent dropout may be unedited tape.** Negligible in Mouse3 (0.34%).
- **A single rate is decisively rejected**, $\chi^2/\mathrm{df}$ = 15,836–86,154 across arms
  ($\hat\Lambda$ = 2.87 Initial, 4.69 Subclone, 5.39/5.52/5.69 Mouse1/2/3).
- **Per-tape rates are real and large**: per-tape mean depth 1.33–5.54 (Mouse3), fitted $\Lambda_z$
  3.71–7.62. $\mathrm{corr}(\beta_z,\bar L_z)$ = −0.63/−0.59/−0.33/−0.59 (Mouse3/Mouse1/Subclone/
  Initial), **unchanged** in the best-captured decile (−0.647/−0.572/−0.331/−0.577). Not a
  cell-capture artefact. ⚑ SciPhy fits per-tape clock rates with shared $\xi$ — same decision,
  reached independently.
- ⭐ **Site 6 is mechanistically different, and it survives the per-tape fit.** obs/exp at depths
  4/5/6: Mouse3 0.99/**2.09**/**0.82**, Mouse1 1.04/1.85/0.84, Initial 1.02/1.34/0.66. A pile-up at
  five filled sites and a deficit at six, which tape heterogeneity cannot make. Independently, $q$ is
  elevated at site 6 in 5/5 arms. Two unrelated signals.
- ⚑ **The arms are a natural two-timepoint experiment**: Initial $\Lambda\approx2.87$ vs mice
  5.39–5.69, ratio $\approx1.9$. With the durations this identifies $\lambda$ and tests rate constancy.
- ⚠⚠ **Editing and dropout are not separately fittable.** Missingness is informative for depth
  ($\rho$ to −0.63), so fitting $\Lambda$ on clean reads conditions on a depth-dependent selection
  and is biased upward — visible as the per-tape fit over-predicting depths 1–3 by ~2× in Mouse3.
  This **revises** the "calibrate editing first, then dropout" ordering proposed earlier the same day.
- ⚠ Unresolvable in Park: "tape $z$ edits less" vs "tape $z$ is transcribed less and read less
  fully". The control conditions on *cell* capture, not per-tape read depth, and Park ships no
  molecule counts. dtt-mouse does (`n_loci`, `mean_dominance`).

**SciPhy Session 3 close read done — full write-up at `notes/sciphy_notes.md` §S3.**
Their recipe, the fact that they already simulate both dropout axes but filter rather than model
them, the 260-line Java simulator that takes the tree as an *input*, the language decision (Python
for simulation; BEAST only if SciPhy inference is wanted), and the birth–death sampling model derived
in full with the Riccati/pgf route. Everything verified by Monte Carlo.

**⭐⭐ Two results that settle the tree model.**
1. **Shape and time factorise.** Conditional on tip count, topology is independent of $b$ and
   $\delta$ — verified across turnover 0 → 0.75, root-split distributions agreeing with Yule–Harding
   to ≤2.0 se. ⇒ shape is parameter-free; all of $b,\delta,\rho$ act on the branching times.
   ⚠ "Random binary tree" is ambiguous — PDA ≠ Yule, and PDA is far more imbalanced. Clade sizes set $m$.
2. **A single homogeneous birth–death cannot produce Park's clone sizes.** The geometric tail
   predicts $2.8\times10^{-14}$ (Mouse3) to underflow (Mouse2) clones as large as the largest
   observed; one is observed in every arm. ⇒ **do not simulate clone sizes — take them from the data
   and condition the within-clone tree on observed size.** Removes the $(b,\delta,\rho)$ fit; with
   $\rho$ also non-identifiable, the tree contributes **one effective number**.

**⇒ NEXT:** derive and verify the conditional branching-time density (route (c) in §S3.3.5), checked
against forward+reject at small $n$; then the joint editing/dropout fit.


## 2026-09-20 — session 14: the forward tree simulator built and validated

**Scope reversal, recorded.** §S3.3.5 made direct BDS-density sampling "the one to build" with
forward+reject as validation only. On the stated priority — **purposes I and IV, IV leading** — that
is backwards: the design sweep has $n$ and $\rho$ as *design variables*, so it needs no conditioning,
no rejection and no density. The $n/\rho$ blow-up is entirely an artefact of Park's
$\rho\approx8\times10^{-4}$. ⇒ build forward, **measure** the cost, reach for the density only where
the measurement says so. Justin's call, after pushing back on why $P(\text{tree}\mid n)$ is needed at
all — the answer being that conditioning on $n$ is how we avoid modelling clone-size heterogeneity
(§S3.4), and that $T$ is not a free size knob because it is shared with the editing layer through
$\Lambda=\lambda T$.

**Built** (`analyses/2026-09_simulator/src/`): `01_bdtree.py` (Gillespie loop → genealogy →
$\rho$-sampling → reconstruction, plus a naive oracle), `02_validate_tree.py` (five checks),
`03_cost_curve.py` (cost measurement — **not yet run**).

**⭐ The two-pass design, and why it is exact.** The lineage each event lands on is uniform and
independent of everything else, so the law factorises as (count trajectory) × (lineage assignments |
trajectory), and the acceptance test reads only $N(T)$ and a $\mathrm{Binomial}(N(T),\rho)$ draw.
Pass 1 therefore gets $N(T)$ with no genealogy at all, vectorised; the $\approx2.7n$ rejected
attempts never allocate a node. Pass 2 replays the identical trajectory from the same seed.

**✅ VALIDATION PASSED at both $\rho=0.5$ and $\rho=8\times10^{-4}$** — 4,000 trees per $n\in\{4,5,8\}$,
zero structure faults in 24,000 trees, worst statistic at 0.64 / 0.71 of its own p99 critical value.
Pruning load (nodes traversed per branch point) **3.6–3.8 vs 22.0–25.4**, confirming the two runs
exercised genuinely different regimes.

**⚠⚠ THREE DEFECTS FOUND, ALL IN MY OWN TEST OR READING RULE, NONE IN THE SIMULATOR.**
1. **`clade_sizes` accumulated in the wrong order.** Descending node id is right for the full tree
   (children get larger ids) and exactly wrong for a `ReconTree` (slots are handed out from the root
   down, so a child has a *smaller* id). Parents were totalled before their children and counted
   only their direct leaf children. First run returned **FAIL at 63.25 se**. ⚑ Diagnosed by an
   internal contradiction: the labelled-history check *passed* in the same run, and root split is a
   function of the labelled history. ⚑ $n=4$ passed and masked it — the broken computation is
   coincidentally right for both 4-tip shapes. Fixed by ordering on node **time**, correct under
   either id convention. ⚠ **Check E saw nothing**, because it compares two implementations through
   the *same* measurement function, so a bug in the shared instrument cancels.
2. **The flat 3-se flag was miscalibrated.** Every check reports the MAXIMUM $|z|$ over its cells, so
   the threshold must grow with cell count. Null p99 of $\max|z|$ at $N=4000$: **2.56 ($m$=2), 3.03
   (4), 3.44 (18), 4.42 (180)**. The flat flag called `rho8e4` a FAILURE at 3.14 on a 180-cell check
   whose statistic sat at the **62nd percentile** of its own null. Now calibrated per check by
   multinomial Monte Carlo. ⚠ $\sqrt{2\ln 2m}$ is not an adequate substitute — it reads 2.04 where
   the null mean is 1.38, because cell counts at these $N$ are skewed, not normal.
3. **`03`'s cost null was turnover-blind** (caught before running). $2.7n^2/\rho$ is the $\theta=0$
   case; the turnover factor $(1+\theta)/(1-\theta)^2$ is **1.0 / 2.65 / 28.0** at $\theta$ = 0 /
   0.3 / 0.75, so the flat null would have printed a spurious VOID at high turnover. The exponent
   $\gamma$, which is the load-bearing test, was never affected — the factor is $n$-independent.

⚑ **Free check on the cost model**: analytic $P(k=n)$ gives 32.1 attempts per accepted tree at
$n=8,\rho=0.5$ against the envelope's $2.7n/(1-\alpha)=30.3$. The *rejection* half of the cost model
is confirmed; the population half is what `03` still has to measure.

**⚠ Owed:** B/C/D/E run at one turnover only (0.3) — §S3.3.3's claim is that topology is
turnover-*independent*, and one $\theta$ does not test that the simulator preserves it.
**⇒ NEXT: the `03` cost pilot**, then $R$ from the measured timing.


## 2026-09-20/21 — session 14 (cont.): validation completed, cost measured, density dropped

**✅ VALIDATION COMPLETE — PASS in every configuration.** 84,000 trees, **zero structure faults**,
across 3 turnovers × 3 $n$ × 2 capture fractions plus four large-$n$ structural points ($n$ = 210,
1,024, 3,387). Worst statistic 0.96× its own p99 critical value; all other runs 0.64–0.71×.
⭐ Pruning load spans **3.4 to 65.2** nodes per branch point (19×, monotone in turnover), so the
degree-2 suppression path is tested across the whole regime, not at one point.
⚠ Per-**check** multiplicity is still unhandled (~15 checks per job each at p99 ⇒ ~14% of correct
runs trip one); the 0.96 near-miss is exactly that rate, not evidence of a defect.

**✅ COST MEASURED ⇒ ⭐⭐ THE BDS DENSITY IS DROPPED.** Nothing in this project exceeds one Slurm
task. The decision quantity is **not** total core-hours (embarrassingly parallel over clones and
replicates) but the **atomic unit**: one tree of the arm's largest clone. Worst case in the whole
project is Subclone at Park $\rho$, **70.7 h**, inside the 7-day cap — and that is on an inflated
clone size (≈12 h on the recorded 10,997). ⇒ §S3.3.5 route (c) and the owed branching-time density
verification are **no longer needed**.
⭐ **The cost model validated OUT OF SAMPLE**: fitted on $n$ = 4, 32, 210, it predicted check F's
independent $n$=1,024 point to 0.91× (seconds) and 1.00× (peak live lineages).

**⚠⚠ FOUR MORE DEFECTS, ALL IN MY TEST / READING RULE / PROJECTION, NONE IN THE SIMULATOR.**
4. **`--f-grid` was word-split by `submit.sh`** (`$*` inside `--wrap` drops quoting) — both pass-2
   jobs died in 9 s. Separator changed to commas.
5. **Clones were grouped by ORGAN, not arm.** `clone_sizes` stripped only a trailing `_<digit>`,
   leaving M1_LL, M1_LN, M2_LV… Park clones span organs by construction (§D.4c), so this **split
   single clones and understated their size**. Fixed to `Sample.split("_")[0]`.
6. **A single mean $\tau$ was applied across the whole range.** Seconds per lineage-event ran
   $6.9\times10^{-6}$ ($n$=4) to $4.4\times10^{-8}$ ($n$=3,387) — **155×** — because at small $n$ an
   accepted tree is many cheap *attempts* (Python overhead) and at large $n$ the vectorised per-event
   cost dominates. The mean, set by the small-$n$ points, inflated the $n$=27,224 projection **~36×,
   enough to flip a verdict**. Replaced by a two-term fit, seconds = $c_1\cdot$attempts +
   $c_2\cdot$events. ⭐ The fix validates itself: $c_2$ = 3.47e-8 vs 3.83e-8 s fitted independently
   on regimes differing 625× in $\rho$.
7. **The verdict thresholds compared core-hours to a WALLTIME cap.** Total core-hours is a budget
   question; the atomic unit is the feasibility question. This is what changed the answer.

**⚑ Storage, measured.** ~**12 bytes per tip** compressed, so a full Park-scale library at $R$=100 is
~120 MB per configuration — a non-issue. Attempts are **independently seeded**, so an accepted tree
replays from (seed, attempt index) in one attempt: verified **bit-identical**, 55× cheaper at
$n$=1,024 and ~4,100× at Subclone scale. ⚠ Store **branch lengths, not node times**: float32 times
differenced give ~9% error on the shortest measured branch ($1.4\times10^{-6}T$) by cancellation,
while a branch length itself is safe in float32 to $3.6\times10^{-8}$.

**⇒ NEXT:** $R$ (now that timing exists), the tree library (proposed, unapproved), then the editing
layer. ⚠ Clone sizes still lack the paper's per-cell tape filter — quote the projection as an
over-estimate until it is applied.


## 2026-09-21/22 — session 15: the tree library built

**⭐ Design change, Justin's call: stock the RANGE of clone sizes, not Park's clone list.** Initial
alone has 2,544 clones, nearly all tiny. An arm is reassembled later by drawing, for each real clone
of size $n$, the stocked tree nearest to $n$. ⚠ Range is what we stock, **distribution is what we
draw** — the statistics this feeds are strongly clone-size dependent (best-$k$ inverts with clone
size; the clone floor was monotone 4/4), so averaging over the grid would destroy the comparison.
⚑ Refinement the instruction implied: **97.1% of Park clones are $\le$64 cells** (median 8), and
stocking every integer size there costs **1.6 core-h**, so substitution is now EXACT for 97% of
clones rather than up to 34% off.

**Built.** `04_tree_library.py` + `05_submit_library.sh`. Grid: 20 log-spaced sizes 2–11,081 pinning
the five observed arm maxima, plus every integer 2–64 · $\theta\in\{0.3\ldots0.7\}$ ·
$\rho\in\{0.25,0.10\}$ everywhere and $\{0.01,0.002\}$ capped at the mouse size range · $R$=20.
**29,640 trees, ZERO structure faults, 94 MB.** Each tree stores `parent` int32 + **branch lengths**
float32 (not node times — float32 differencing costs ~9% on the shortest branch), seed + accepting
attempt index (replay verified **bit-identical**), a hash, and edit-free summaries (tree length, B1,
LTT). ⚠ Prefix-clade-size-by-depth is deliberately NOT computed — it needs the editing layer.

**⭐⭐ THE FINDING: $\rho$ dominates $\theta$ by 3–5× on coalescent depth.** Median coalescent
depth (fraction of $T$ by which half the lineages exist; larger = shallower tree), at matched $n$:

| | span over $\rho$ 0.25→0.002 | span over $\theta$ 0.3→0.7 |
|---|---|---|
| $n$=63 | **+0.338** | +0.082 to +0.111 |
| $n$=3,387 | **+0.268** | +0.051 to +0.062 |

Monotone in the predicted direction in **20/20** $\theta$ rows and **5/5** $\rho$ comparisons — the
pull of the present, measured. ⇒ **The unknown I could not resolve from the paper (mouse $\rho$)
matters several times more than the one I derived from growth arithmetic ($\theta$).** That sharpens
the ask to Jihye Park: total viable cells per dissociation, or tumour mass, per organ at day 45.
⚑ The $\theta$ effect also SHRINKS with clone size (0.094 → 0.051 from $n$=63 to 3,387), so small
clones are the more informative ones about turnover. ⚠ On this one summary only; the real comparison
runs through prefix clades and inherits the editing layer's noise.

**⚠⚠ DEFECT 8 — the cost model was fitted at turnover 0.3 and applied to 0.7.** `03`'s `project()`
filters on `turnover == 0.3`. Measured on this run, actual/predicted rises monotonically with
turnover — **median 0.79 at $\theta$=0.3, 1.27 at 0.7, tail to 4.12×** — because total nodes
$\approx 2N(T)/(1-\theta)$ is 2.9× $N(T)$ at 0.3 but 6.7× at 0.7, and those extra nodes go through
the PYTHON genealogy loop (~1 µs/event) not the vectorised pass-1 (~38 ns/event). Aggregate was fine
(396 core-h actual vs 339 predicted, 1.17×); the **tail** hit the 12 h walltime on **4 of 120 tasks**,
losing 12 of 1,490 cells. Backfilled via a new `--only` flag at 17× walltime margin. ⇒ size
high-turnover work generously; the caveat is recorded at the top of `04`.

**⇒ NEXT:** the editing layer, then fit $\theta$ (and $\rho$) jointly against clade-size-by-depth.
