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
