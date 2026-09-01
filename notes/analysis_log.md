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
