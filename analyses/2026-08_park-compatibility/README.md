# park-compatibility

**Question.** Does Park's cross-tape character compatibility support the perfect-phylogeny route
(§D.4), and how big is the gap between the homoplasy prediction and the measurement?

**Answer.** Partly, and the diagnostics matter more than the answer. Compatibility is 94.66%
(missing-excluded) / 63.56% (missing-as-absent) on Mouse3, and only the pessimistic end is
constructible (§D.4d) — **dropout is the binding constraint, row A6**. The skeleton was then set
aside as a step sideways (see the strategic reassessment below); the diagnostics are the deliverable
and they motivate the likelihood route. Dropout is now characterised in its own right — see
"Fig 3 redesign".

## What was run

- `src/01_symbol_composition.py` — streams all five edit tables, reports shape, missingness, the
  insertion alphabet, and $q=\sum_i \xi_i^2$. Output: `results/symbol_composition.json`.
- `src/02_alphabet_checks.py` — rarefaction of $M_{\rm obs}$, and per-site $q$ as a test of row A3.
  Output: `results/alphabet_checks.json`.
- `src/03_missing_pattern.py` — separates tape-dropout from unedited sites.
- `src/04_xi_vectors.py` — full $\xi$ vectors, global and per site. Output: `results/xi_vectors.json`.
- `src/05_coupon_identifiability.py` — $\mathbb{E}[s\mid m]$, $\mathrm{sd}[s\mid m]$, and the inversion.
  Outputs `results/coupon_identifiability.json`, `figures/coupon_identifiability.png`.
- `src/06_estimate_m.py` — count-only $\hat m$ per prefix node. `results/m_estimates_*.tsv.gz`.
- `src/07_m_figure.py` — `figures/m_distribution.png`.
- `src/08_poisson_mle.py` — **the set-dependent MLE, with the full derivation in its header**;
  MC validation + per-node estimates. `results/m_mle_*.tsv.gz`, `results/mle_validation.json`.
- `src/09_mle_vs_moment.py` — estimator comparison and the (corrected) model check.
- `src/10_build_characters.py` — character sets per clone, with the determination mask $D$.
- `src/11,12_junk_*.py` — the two tests that promoted junk values to real symbols.
- `src/13_compatibility.py` — reference cross-tab engine; `src/14_compat_sparse.py` — production
  engine (two sparse products); `src/15_memory_bound.py` — the analytic sizing bound.
- `src/16_skeleton_C.py`, `src/17_skeleton_linear.py` — skeleton construction; `18,19` — the subclone
  ground-truth test and group relatedness.
- `src/20,21_fig{1,2}_*.py` — figures 1 and 2. `src/22_fig3_dropout.py` — **superseded**, kept as the
  record of the old Fig 3.
- `src/23_dropout_matrix.py` — caches the (cell × tape) recovery / depth / termination matrices;
  every Fig-3 script reads this so the CSVs are parsed once.
- `src/24_dropout_margins.py` — both margins, VIF and $\rho$, unfiltered and QC-filtered.
- `src/25_dropout_margins_followup.py` — the shelf, the common ≥100 cut, barcode loss vs $R_c$,
  per-tape rate reproducibility across libraries and arms.
- `src/26_dropout_depth.py` — depth coupling on both margins, each with its confound controlled.
- `src/29_rho_check.py` — $\rho$ verified three independent ways; Spearman shown ≡ Pearson for binary.
- `src/27,28,30,31` — Fig 3 panels a, b, c, and d+e. **No assembly script by design** — the panels
  are presented individually.
- `src/32_lineage_feasibility.py` — power/feasibility for the A9 test: usable cells per arm, clone
  size bands, informative absence characters, and the clone/sample purity that sets the batch confound.
- `src/33_collision_screen.py` — Park's cross-clone collision screen, vectorised, cached as per-cell
  flags. Needed because the Fig-3 cache carries no symbols.
- `src/34_lineage_rho.py` — **T0/T1**: the intraclass-correlation machinery, three centrings, the
  edit-feature calibration curve, $\rho_{\rm clone}$ on $R_c$, and the within-sample permutation null.
- `src/35_lineage_depth.py` — **T2**: the relatedness gradient below the clone, with the analytic
  within-clone null and its permutation verification.

## Findings

### Step 0 — the matrix is present and is what §D.4b needs

| table | cells | tapes | sites | missing | $M_{\rm obs}$ | edits | junk | $q$ | $1/q$ |
|---|---|---|---|---|---|---|---|---|---|
| Initial  | 37,810 | 166 | 6 | 62.2% | 244 | 13,892,618 | 2.47% | 0.016905 | 59.2 |
| Mouse1   | 12,232 | 166 | 6 | 51.2% | 197 |  5,809,280 | 2.36% | 0.017400 | 57.5 |
| Mouse2   |  6,899 | 166 | 6 | 50.0% | 184 |  3,339,410 | 2.72% | 0.017168 | 58.2 |
| Mouse3   |  2,904 | 166 | 6 | 53.2% | 167 |  1,322,277 | 2.29% | 0.017256 | 57.9 |
| Subclone | 39,606 | 166 | 6 | 42.7% | 218 | 22,116,803 | 2.22% | 0.016889 | 59.2 |

† *these `missing` figures pool tape-dropout with unedited sites — see the correction below.*

99,451 cells total. **166 tapes × 6 sites confirmed** against the notes' record of the design
(§1629 table), with a byte-identical tape-barcode set across all five tables — so the tables are
directly stackable and tape identity is consistent.

### ⚠⚠ $q$ is 4.3× the value assumed throughout the notes

$q\approx0.0170$, not the $\approx0.004$ that §D.4b and §1681 assume. The old figure took a flat
distribution over $M=256$; in fact only 167–244 symbols appear and the distribution is skewed
(top symbol `ATATGGA` = 4.4%), giving an **effective alphabet of ~58, not 256.**

This *inverts* the comparison to Mulberry: their "high diversity" benchmark is $q=1/64=0.0156$, so
Park is marginally **worse**, not "near homoplasy-free". Correction written into §D.4b.

**Consequence for this analysis:** the §D.4b null must be recomputed at $q=0.0170$ before the
predicted-vs-measured gap can be read as dropout/error. The residual framing is unaffected; the
subtrahend changes. Do this before running the compatibility check.

### Two data-handling facts the protocol has to absorb

1. **⚠ CORRECTED — the 43–62% "missingness" was two different things added together.**
   `None` marks both *tape not observed* and *site not yet edited*, and the Step-0 number pooled
   them. Decomposed on Mouse3 (482,064 cell×tape instances, `src/03_missing_pattern.py`):

   | | tape-instances | % of all entries |
   |---|---|---|
   | tape entirely absent (**true dropout**) | 212,802 (44.14%) | 44.14% |
   | observed tape, trailing unedited sites (**biology, not loss**) | — | 9.07% |
   | actual edits | — | 46.79% |

   Only 0.02% of tapes show an *internal* `None` (a gap with a symbol after it), so the sequential
   architecture is respected essentially perfectly — a strong sanity check on the whole dataset.

   ⇒ **True per-tape dropout is 44%, not 53%.** The other 9% is saturation, and it is informative
   data, not missing data. Step 4's missing-as-absent/missing-excluded split applies to the 44%
   only; treating unedited sites as "missing" would be a modelling error.

   ⚑ **Mean fill is 5.02 of 6 sites; 41.4% of observed tapes are completely full.** The recorder is
   heavily saturated in Mouse3 — relevant to dynamic range (§H.6.13) and to how much of the
   experiment the tapes actually witnessed.
2. **~2.3% of observed entries are not NNNNGGA** — 1,939 distinct malformed strings, lengths 1–132.
   The largest single class is 6-mers (21,244 in Mouse3), dominated by `CACGGA`, `GATGGA`,
   `GCGCGG` — consistent with truncation. **Undecided: drop these, or treat as a distinct
   observed state?** They are not missing data and silently coercing them to `None` would inflate
   the dropout estimate. Decide before Step 2.

### Sample structure

From `clonalbc_percell_hamming1_corrected.csv`'s `Sample` column, joined to the CellID suffixes:

| table | samples | reading |
|---|---|---|
| `Initial` | `Initial_1`, `Initial_2` | pre-transplant in vitro population, 2 replicates |
| `Mouse1` | `M1_LL`, `M1_LN`, `M1_TH`, `M1_LLL`, `M1_R`, `M1_IP` | one xenografted mouse, dissected by anatomical compartment |
| `Mouse2` | `M2_LV`, `M2_LL1`, `M2_LL2` | " |
| `Mouse3` | `M3_LL1`, `M3_LL2`, `M3_LV` | " |
| `Subclone` | `Subclone_1`, `Subclone_2` | the subclone arm, 2 replicates |

So the five tables are **experimental arms, not five clones** — `Mouse1/2/3` are the three
xenografted animals of Park, Chang et al., each split across harvest sites. The two-letter codes
are anatomical (LL/LLL lung, LN lymph node, LV liver, TH, R, IP) — **expansion inferred from
context, confirm against the paper before using in any writeup.**

⚑ Supporting evidence that `Initial` is the earliest timepoint: its site-6 occupancy is far lower
than the others ($n_6/n_1 = 5.1\%$, vs 38% Mouse1, 27% Subclone), i.e. the least tape saturation,
as expected for a pre-transplant sample.

### $M_{\rm obs}$ is mostly a sampling artefact — rarefaction

The tables differ 17-fold in edit count, so raw $M_{\rm obs}$ is not comparable.
Rarefied to Mouse3's 1,322,277 edits (`src/02_alphabet_checks.py`):

| table | edits | $M_{\rm obs}$ | $M$ @ rarefied |
|---|---|---|---|
| Initial | 13,892,618 | 244 | 201 |
| Mouse1 | 5,809,280 | 197 | 168 |
| Mouse2 | 3,339,410 | 184 | 163 |
| Mouse3 | 1,322,277 | 167 | 167 |
| Subclone | 22,116,803 | 218 | **157** |

The 167–244 spread narrows to 157–201, and the ordering changes — Subclone goes from second-highest
to lowest. **Do not read $M_{\rm obs}$ as a property of a sample.** $q$ is the robust statistic:
it varies by <3% across tables while $M_{\rm obs}$ varies by 46%.

### Row A3 partially fails across sites — site 6 has elevated $q$

$q$ per site (sites are ordered in time, so this is a direct test of A3):

| table | Site1 | Site2 | Site3 | Site4 | Site5 | Site6 |
|---|---|---|---|---|---|---|
| Initial | 0.01720 | 0.01688 | 0.01679 | 0.01662 | 0.01708 | **0.01773** |
| Mouse1 | 0.01807 | 0.01786 | 0.01748 | 0.01795 | 0.01874 | **0.02008** |
| Mouse2 | 0.01795 | 0.02024 | 0.01922 | 0.01848 | 0.01765 | **0.02290** |
| Mouse3 | 0.01752 | 0.01739 | 0.01726 | 0.01763 | 0.01740 | **0.01841** |
| Subclone | 0.01698 | 0.01869 | 0.01681 | 0.01806 | 0.01674 | **0.01806** |

Site 6 is the maximum in **all five** tables, 3–27% above the same table's minimum.

Not a small-sample artefact: the plug-in estimator has $\mathbb{E}[\hat q]=q+(1-q)/n$, so at
$n_6\approx10^5$–$2\times10^5$ the upward bias is $\approx5\times10^{-6}$ — three orders of
magnitude below the observed 0.001–0.005 elevation.

⇒ **$\xi$ is not constant across sites**, so the pooled Step-0 $\hat\xi$ is a
time-averaged quantity. The effect is modest (few %) and does not change the headline
$q\approx0.017$, but it means A3 is measurably violated in the direction that matters, and any
per-site homoplasy null should use the per-site $q$.

Edit counts fall monotonically $n_1>n_2>\dots>n_6$ in every table — the expected sequential-filling
signature, and a clean sanity check that the site ordering in the columns is real.

### Is $m$ recoverable from $s$? — yes, over the range that matters

$m$ (independent write events at a site within a prefix clade) is not observable; $s$ (distinct
symbols those events produced) is. They differ by exactly the homoplasy count, $m-s$, which is the
circularity. The way out is that the $m$ events are iid draws from $\xi$ (§1a.3), so

$$\mathbb{E}[s\mid m]=\sum_i\left(1-(1-\xi_i)^m\right)$$

is monotone in $m$ and invertible against the **measured** $\xi$. `src/05_coupon_identifiability.py`
computes this plus $\mathrm{sd}[s\mid m]$ (occupancy variance, incl. the negative cross-covariance)
and inverts the $\pm1$ sd band to get the recoverable interval $[m_{\rm lo},m_{\rm hi}]$.

**Both closed forms verified against Monte Carlo** (20,000 reps to $m$=2,574; 4,000 above) — $\mathbb{E}$
and sd agree to 3–4 significant figures at $m=5,11,74,313,2574,10506$.

| quantity | value |
|---|---|
| $m$ recoverable to better than 2× | $m \lesssim 1{,}199$ |
| $m$ has no upper bound above | $m \approx 5{,}997$ ($s$ saturated) |
| $m$ at half the alphabet | 89 |
| $m$ at 90% of the alphabet | 602 |
| $\mathbb{E}[s]$ at $m=74$ (clone size) | 44.9 |
| $\mathbb{E}[s]$ at $m=10{,}506$ (largest subclone) | 99.8 of 100 — saturated |

⇒ **Usable for the primary analysis unit.** §D.4b step 1 works within clones of ~74 cells, where
$m\le74$ sits comfortably below the resolution limit. It degrades in the largest subclones only,
where we fall back to the assumption-free bracket $s\le m\le n_L$.

### ⚑⚑ The real alphabet is ~100 symbols, not 256 — and the tail nearly fooled the estimator

Symbol counts in Mouse1 split cleanly:

| count band | # symbols | share of edits |
|---|---|---|
| >10,000 | 91 | 99.392% |
| 1,001–10,000 | 9 | 0.559% |
| 101–1,000 | 6 | 0.030% |
| 11–100 | 30 | 0.015% |
| 2–10 | 40 | 0.004% |
| 1 | 21 | 0.0004% |

**100 symbols carry 99.95% of all edits.** The other 97 carry <0.05% between them — rarest
$\xi=1.7\times10^{-7}$, i.e. seen *once* in 5.8 M edits. Those are sequencing artifacts that happen
to satisfy NNNNGGA, not real pegRNAs. So the design alphabet $M=256$ is not realised: the effective
pegRNA pool is ~100, with $1/q\approx57$ after skew.

⚠ **This tail is a trap for exactly this estimator.** Being rare, artifact symbols keep $s$ growing
long after the real alphabet is exhausted, so the untrimmed 197-symbol $\xi$ makes the inversion
*look* well-conditioned out to $m\sim10^6$. It isn't — that apparent resolution is entirely
manufactured by noise. Panel A of the figure shows both curves; the grey one is the illusion.
$q$ is unaffected (0.017400 untrimmed vs 0.017417 trimmed), which is again why $q$ is the robust
statistic and symbol *counts* are not.

![coupon identifiability](figures/coupon_identifiability.png)

### First pass: $\hat m$ per prefix node — homoplasy is rare, and concentrated

`src/06_estimate_m.py` walks the prefix trie within each clonal barcode, per tape, and records
$s$ and the clade size at every node; `src/07_m_figure.py` summarises. **1,567,321 nodes** with
≥3 carriers across all five tables.

One alphabet rule used consistently for both $\xi$ and $s$: a value is a symbol iff it matches
NNNNGGA **and** occurs ≥1,000 times in its table (~94–106 symbols, 99.7–99.97% of edits).
Anything else — scaffold read-through, indel variants, the rare tail — truncates the prefix.

![m distribution](figures/m_distribution.png)

**1. 95.8% of nodes sit below the birthday threshold** $m^\ast=\sqrt{2/q}=10.8$. Median $s$ is 1–3
almost everywhere: the carriers of a prefix overwhelmingly share *one* next symbol, which is the
single-origin behaviour a perfect phylogeny predicts.

**2. Homoplasy scales with clade size, as the mechanism predicts** — mean recurrences per node:

| clade size | 3 | 4–5 | 6–10 | 11–20 | 21–50 | 51–100 | 101–500 | >500 |
|---|---|---|---|---|---|---|---|---|
| mean recurrences | 0.01 | 0.02 | 0.07 | 0.20 | 0.53 | 1.65 | 5.75 | **50.49** |
| % of nodes with ≥1 | 0.0 | 0.0 | 0.0 | 2.0 | 18 | 41 | 49 | **55** |

⚠ This is *consistency*, not independent confirmation: $\hat m - s$ and $\binom{m}{2}q$ are related
by construction, since $\mathbb{E}[s\mid m]=m-\mathbb{E}[\text{collisions}]$. What is genuinely
informative is the **clade-size dependence**, which the estimator does not build in.

**3. ⚑⚑ It is concentrated, which is the answer §D.4b step 5 wanted.** Ranking nodes by recurrence
count: the **top 1% carry 73%** of all estimated homoplasy, the top 10% carry 92%. Per §D.4b's
decision rule that is the *easy* world — conflict that can be removed by deleting a small set of
characters, rather than conflict spread thin enough to make maximal-compatible-set genuinely hard.

**4. The arms differ enormously, and it is all clade size.** Level-0 nodes (site 1, whole clone):

| table | median clade | mean recurrences | % nodes with ≥1 |
|---|---|---|---|
| Mouse1/2/3 | 7–8 | 0.03–0.11 | 0.4–1.1% |
| Initial | 9 | 0.22 | 5.4% |
| Subclone | **933** | **16.0** | **40.8%** |

The subclone arm is where homoplasy lives, and it is there because subclones were *designed* to be
large. The metastasis mice are essentially homoplasy-free at the clade sizes they actually have.

### ⚠ Two data facts that complicate the per-clone protocol

**Clone sizes do not match the notes.** §1629 records "~74 cells/clone × 75 clones". The delivered
`ClonalBC` column has **3,294 distinct barcodes, median 7 cells, max 27,537**, with five clones
holding 50% of all cells. Whatever the "75 clones" refers to, it is not this column as delivered —
resolve before per-clone results are quoted against the paper.

**Clone-barcode dropout is substantial and uneven**: cells skipped for lacking a `ClonalBC` are
2,239/37,810 (5.9%) in Initial, 670/39,606 (1.7%) in Subclone, but **5,489/12,232 (44.9%) in
Mouse1**, 1,511/6,899 (21.9%) Mouse2, 970/2,904 (33.4%) Mouse3. The mouse arm loses a third to a
half of its cells before the analysis starts.

### Second pass: the Poissonised MLE — implemented, validated, run

The current estimator uses only $|A|=s$, discarding *which* symbols were seen. It shouldn't: the
likelihood of observing set $A$ in $m$ draws is
$P(A\mid m)=\sum_{B\subseteq A}(-1)^{|A|-|B|}W_B^{\,m}$ with $W_B=\sum_{i\in B}\xi_i$, which depends
on the mass of $A$, not just its size. A rare-symbol $A$ implies $\hat m\to s$; a common-symbol $A$
tolerates much larger $\hat m$.

Poissonising the draw count makes the symbols independent and the likelihood exact and tractable:

$$\ell(m)=\sum_{i\in A}\log\!\left(1-e^{-m\xi_i}\right)\;-\;m\!\!\sum_{i\notin A}\!\xi_i$$

The second term is the set-dependent penalty: large $m$ is punished in proportion to the mass
**not** observed. Full derivation — inclusion–exclusion, the Poissonisation step, concavity, the
saturation limit — is documented in the header of `src/08_poisson_mle.py`.

**Monte-Carlo validated against the true fixed-$m$ model** (`results/mle_validation.json`):

| true $m$ | mean $s$ | MLE bias / IQR | count-only bias / IQR |
|---|---|---|---|
| 5 | 4.8 | +4.1% / 2.3% | +3.7% / 0.0% |
| 50 | 34.8 | +0.6% / 14.2% | +0.7% / 16.8% |
| 200 | 71.5 | +0.5% / **17.6%** | −1.8% / 21.0% |
| 500 | 87.9 | +1.6% / **26.4%** | −0.1% / 35.0% |

Both estimators are near-unbiased, so Poissonisation costs ≤1.6% bias. The MLE's gain is
**variance**: 16% narrower IQR at $m=200$, 25% narrower at $m=500$. It is no better below $m\approx10$.

**On the data** (`src/09_mle_vs_moment.py`, all 1,567,321 nodes):

| $s$ | nodes | $\hat m_{\rm MLE}/\hat m_{\rm count}$ |
|---|---|---|
| 1–25 | 1,559,225 (99.5%) | 1.003–1.007 |
| 26–50 | 6,211 | 0.989 |
| 51–80 | 1,474 | 0.950 |
| 81–120 | 411 | **0.799** |

Identical where $s$ is small — as predicted, since $(1-W_A)\approx1$ there whatever was seen — and
divergent only near saturation, where $\hat m$ scales roughly as $1/(1-W_A)$ and is therefore acutely
sensitive to the observed mass that the count-only estimator throws away.

![mle vs moment](figures/mle_vs_moment.png)

⚠ **A trap worth recording.** The first version of the model check compared observed $W_A$ against
$\mathbb{E}[W_A\mid m]$ at the $m$ implied by $s$, and appeared to show a large-$s$ deficit — i.e.
model misspecification. That was an artifact: the data are selected on $s$, not $m$, and for fixed
$m$ a node reaching an unusually large $s$ got there by hitting unusually many *rare* symbols, which
depresses $W_A$. Conditioning the null on the realised $s$ instead (by simulation), **observed $W_A$
falls inside the 90% band at every $s$** — no misspecification detected. Panel B shows the corrected
comparison.

**Effect on the conclusions: none qualitatively.**

| | below $m^\ast$ | total recurrences | top 1% holds | top 10% |
|---|---|---|---|---|
| count-only | 95.82% | 629,428 | 73.1% | 91.6% |
| Poissonised MLE | 96.01% | 574,959 | 65.0% | 86.9% |

Homoplasy is 8.7% lower and slightly less concentrated, but the picture — overwhelmingly
sub-threshold, strongly concentrated — is unchanged. **Use the MLE numbers going forward**; it is
the more efficient estimator and costs nothing where the two agree.

### ⚠⚠ CORRECTION — junk values ARE heritable characters; the alphabet rule is wrong

Scripts 06–10 treat any non-NNNNGGA value as junk that truncates the prefix, on the grounds that
`CACGGA` occurs 47,287 times and would "manufacture false clades" if admitted. **That reasoning
assumed its own conclusion**, and the data contradict it (`src/11_junk_heritability.py`).

Two hypotheses, opposite predictions, no tree required. A *heritable* insertion at site $k$ of tape
$z$ arose once in an ancestor, so its carriers must share that ancestor's sites $1..k-1$ and sit in
one clone. A *readout artifact* lands on unrelated cells.

**Test 1 — preceding-prefix concordance** (median, Mouse1, matched on carrier count):

| carriers | 3–4 | 5–9 | 10–24 | 25–99 | 100–499 | 500+ |
|---|---|---|---|---|---|---|
| real symbols | 0.500 | 0.400 | 0.348 | 0.408 | 0.515 | 0.945 |
| junk values | **1.000** | **0.667** | **0.524** | **0.714** | **0.880** | **0.980** |

**Test 2 — clone restriction** (triples with ≥10 carriers):

| | n | top-clone share | top-sample share | clone/sample ratio |
|---|---|---|---|---|
| junk | 1,356 | **0.681** | 0.856 | **0.857** |
| real symbols | 50,813 | 0.456 | 0.783 | 0.625 |
| random cells (null) | — | 0.16–0.20 | — | — |

**Junk is *more* clone-restricted than real symbols**, and far above the random null. The
library-artifact confound is excluded: if junk arose per sequencing run it would be sample-restricted
but *not* clone-restricted, giving a low clone/sample ratio — instead junk's ratio (0.857) is
*higher* than real symbols' (0.625). The top junk triples carry 1,200–1,600 cells spread over ~50
clones with **92–96% in a single clone** — one ancestral event plus a small leakage tail.

Why *more* restricted than real symbols? Junk events are rarer per (tape, site), so fewer independent
origins, so each occurrence traces to one ancestor. That is what a rare heritable event looks like.
It is also mechanistically expected: pegRNA scaffold read-through physically writes scaffold sequence
into the tape, and once written it is as irreversible as any designed insertion.

⇒ **Change the alphabet rule to frequency alone: a value is a symbol iff it occurs ≥1,000 times in
its table.** Drop the NNNNGGA pattern requirement — frequency is the evidence of reality; conformance
to the design pattern is not required, and demanding it discards real heritable events. This is the
same threshold already applied to NNNNGGA symbols (which it uses to exclude ~97 rare artifacts).

Impact on Mouse1: 9 junk values promoted (`CACGGA` 25,757, `GATGGA` 12,294, `GCGCGG` 11,107, four
scaffold read-through variants, `TTTTGGGA`, `TTTTAGGA`), recovering **71.6% of junk edits** directly
plus the **3.03% of downstream symbols** previously lost to truncation. Nearly free in $q$:

| | alphabet | $q$ | $1/q$ |
|---|---|---|---|
| current rule | 100 | 0.017417 | 57.4 |
| frequency-only rule | 109 | 0.017104 | 58.5 |

**Confirmed independently on Subclone** (`src/12_junk_check_table.py`), against a null calibrated
per table — Subclone has 15 clones with one holding 28% of cells, so its random-cell baseline is
0.290 versus Mouse1's 0.238:

| table | | top-clone | null | **excess** | clone/sample |
|---|---|---|---|---|---|
| Subclone (15 clones) | junk | 0.727 | 0.290 | **+0.436** | 1.263 |
| | real | 0.583 | 0.286 | +0.295 | 1.049 |
| Mouse1 (296 clones) | junk | 0.681 | 0.238 | **+0.444** | 0.857 |
| | real | 0.456 | 0.239 | +0.222 | 0.625 |

Junk excess is +0.44 on both tables against different nulls and different clone structures.

**✅ APPLIED 2026-08-31.** Script 04 now counts every value and keeps those with ≥1,000 occurrences;
scripts 05–10 rerun. Alphabet 97–129 symbols per table (was 94–106), 3–23 promoted non-NNNNGGA
symbols carrying 1.6–2.3% of kept edits.

| | old rule | new rule |
|---|---|---|
| mean $q$ | 0.0172 | **0.016641** |
| birthday $m^\ast$ | 10.8 | 10.96 |
| prefix nodes | 1,567,321 | 1,629,739 |
| below $m^\ast$ (MLE) | 96.01% | **96.25%** |
| top 1% of nodes hold | 65.0% | **64.9%** |
| top 10% hold | 86.9% | **86.9%** |
| characters ≥3 carriers | 1,729,592 | **1,801,071** |
| $m$ recoverable to 2× | $m\lesssim1{,}199$ | **$m\lesssim2{,}757$** |

**Nothing qualitative moved.** $q$ falls 3%, so homoplasy becomes marginally *less* likely, and the
concentration result — the one the whole "conflict is removable" conclusion rests on — is unchanged
to within 0.1 points. We gain 4.1% more characters and materially better $m$-identifiability
(a larger real alphabet means $s$ saturates later). The corrected $W_A$ model check still finds
observed mass inside the null band at every $s$.

⚠ Scripts 01–03 predate this rule and still describe the NNNNGGA-only alphabet; they are kept as the
Step-0 record, and the $q\approx0.0174$ quoted in §D.4b CORRECTION 1 is the old-rule value.

### Step 4 — the compatibility check, first table (Mouse3, complete)

`src/13_compatibility.py` is the reference cross-tab engine (verified against explicit brute-force
set operations on real clones, both conventions). `src/14_compat_sparse.py` is the production
engine: the same computation recast as two sparse products,

$$n_{11}=(MM^\top)_{ij},\qquad n_{10}=(MN^\top)_{ij}-n_{11},\qquad n_{01}=(MN^\top)_{ji}-n_{11}$$

with $M$ = characters × cells membership and $N$ = characters × cells determination. Only the
non-zero entries of $MM^\top$ are materialised, so disjoint pairs are counted by subtraction and
never touched. Verified against script 13; **5× faster and it reaches the large clones script 13
could not** (91/91 Mouse3 clones in 560 s vs 63/91 in 915 s).

**Mouse3, all 91 clones, 55,875 characters, 1,856 cells:**

| convention | pairs | incompatible | compatibility |
|---|---|---|---|
| missing-as-absent | 27,368,822 | 9,972,354 | **63.56%** |
| missing-excluded | 27,368,822 | 1,462,105 | **94.66%** |
| | | | **spread +31.09 points** |

By clone size (missing-excluded / as-absent / spread):

| cells | clones | as-absent | excluded | spread |
|---|---|---|---|---|
| 3–4 | 24 | 87.32% | 99.81% | +12.49 |
| 5–10 | 28 | 69.74% | 98.58% | +28.84 |
| 11–20 | 18 | 53.94% | 96.36% | +42.43 |
| 21–300 | 21 | 63.62% | 93.52% | +29.90 |

**Reading 1 — dropout is the binding constraint.** A 31-point spread is the last row of §D.4b's
decision rule: *"dropout is the binding constraint, which is row **A6**, and argues for the SBI
route."* The spread widens with clone size, because more cells give more chances for a dropped-out
cell to fake a witness in $S_1\setminus S_2$. It also vindicates building $D$: the naive convention
reports 63.56% and would have condemned the perfect-phylogeny route outright.

**Reading 2 — 94.66% sits in the 80–95% band**, "live but needs conflict resolution", just under the
≳95% threshold for the easy verdict.

⚠ **Reading 3 — conflict was NOT concentrated, unlike homoplasy.** In the crashed session's run the
conflict graph had ~58% of characters carrying at least one conflict, with the top 10% holding only
~46% of edges — against homoplasy's top 1% holding 65%. If that holds, it is the *"spread thin"*
world of §D.4b step 5, where maximal-compatible-set is genuinely hard. **These degree numbers were
printed but never written to disk and are lost**; `src/14` now persists them
(`results/conflict_degrees_{table}.npz`) and Mouse3 must be rerun to recover them.

**Caveat on all of the above:** Mouse3 is the smallest table — median clone 4 cells, largest 210.
The size trend suggests larger clones will score lower. Mouse1, Mouse2, Initial, Subclone not yet run.

### ⚠⚠ $C$ is not yet measured — and the obvious route to it is invalid

Computing $C$ from the conflict-free characters of the missing-excluded run gave
$C/(n-1)=\mathbf{2.107}$ — impossible — with laminarity failing on 37,737 pairs. Diagnosis in
§D.4d: **missing-excluded compatibility is pair-specific** (each pair is laminar on its own
$D_1\cap D_2$), so it does not compose into the single laminar family a tree requires.

⇒ $C$ must come from the **missing-as-absent** graph, where compatibility does imply laminarity.
`src/14` currently saves degrees only for the excluded convention; it needs to save both, then rerun.
`src/16_skeleton_C.py` is correct machinery pointed at the wrong input.

⇒ Reframes the spread: 94.66% is what we *would* see with complete data; **63.56% is what a skeleton
can actually be built from.** Only the pessimistic end is constructible.

### Verdict so far, against §D.4b's decision rule

| signal | Mouse3 | reading |
|---|---|---|
| compatibility (excluded) | 94.66% | 80–95% band: "live but needs conflict resolution" |
| compatibility (as-absent) | 63.56% | below the <80% "different project" line |
| spread | **+31.09 pts** | **dropout is the binding constraint — row A6** |
| conflict concentration | top 10% hold 45.9% | **"spread thin"** — max-compatible-set is genuinely hard |
| $C$ | **not measured** | the question is not answered until it is |

⚠ Mouse3 is the *most favourable* table: compatibility falls monotonically with clone size (99.81%
at 3–4 cells → 93.52% at 21–300), and its largest clone is 210 cells against Subclone's 10,997.

### Compatibility, second table: Initial

| table | clones | cells | characters | as-absent | excluded | spread |
|---|---|---|---|---|---|---|
| Mouse3 | 91 | 1,856 | 55,875 | 63.56% | 94.66% | **+31.09 pts** |
| Initial | 1,780 | 17,521 | 534,328 | 80.60% | 91.91% | **+11.31 pts** |

Initial's clones are tiny (median 6 cells, max 127), which is why its as-absent number is much
higher and its spread much smaller — fewer cells means fewer chances for a dropped-out cell to fake
a witness. The spread scaling with clone size, established within Mouse3, holds across tables.

**Mouse1 and Mouse2 never completed** — five OOM failures between them (128 G, 96 G, 48 G, 32 G),
all on their single largest clone. See the strategic note below: this is now moot, because the
pairwise matrix is not what a skeleton needs.

---

## ⚑⚑ STRATEGIC REASSESSMENT (2026-09-01) — read this before continuing

**The skeleton is a step sideways, not forward.** Three reasons, in order of seriousness:

**1. The output is arbitrary.** Three reasonable algorithms on identical Mouse3 data gave
$C = 689$, $783$, $1{,}003$ — a **45% spread**. Which skeleton you get depends on tie-breaking
order. Using one as a hard constraint injects a coin flip that *looks* like structure, which is the
opposite of honest uncertainty. (Partial salvage: treat randomised restarts as a sampler and use
per-clade survival frequency as a support measure. Heuristic, not calibrated.)

**2. The premise barely holds for this dataset.**

| | |
|---|---|
| clones >1,000 cells | **9 of 2,547 (0.35%)** |
| of those, in Subclone (an *in vitro* validation arm) | **7** |
| in the actual metastasis experiment | **2** |

SciPhy's practical limit is ~1,000 tips, so **382 of 384 mouse clones are already within reach of
likelihood inference**. Median clone size is 4 cells. The search-space reduction the skeleton exists
to provide is needed for two clones of scientific interest.

**3. The mismatch is structural.** The skeleton is a hard combinatorial device on soft, incomplete
data, and every difficulty traced to that: the missing-as-absent/excluded fork (§D.4d), false
clades, the restarts. Felsenstein pruning has none of them — it sums over unobserved states, so a
cell's ~120 observed tapes determine its position while its missing tape contributes a marginal.
**That is exactly "use other tapes probabilistically to mitigate dropout", and it is not a method to
invent — it is SciPhy's existing likelihood.**

⇒ **The diagnostics are the deliverable.** $q$, the homoplasy quantification, the dropout spread —
these justify the likelihood route and do not need the skeleton to be a production tool. The
skeleton's residual role is soft decomposition or initialisation for the ~9 large clones, never a
constraint.

### Ground-truth test (subclone colonies) — the method passed, better than expected

Pooled 11 ClonalBC groups (200 cells each), built characters blind to labels, built the skeleton,
scored clades against colony identity. After implementing Park's cross-clone collision screen
(3.77% of cells pruned) the median clade purity was **1.000**, with 94% of clades within five cells
of pure.

⚑ The large apparent "false clades" turned out to be **correct**: `ClonalBC` over-splits colonies,
because a founder with multiple barcode integrations yields several groups. Group-consensus
similarity is bimodal — off-diagonal median **0.048**, but **2.18–2.54** for four specific pairs:

- `{ATCCATACGA, CGGATAGTGG, CTATGGTAAG}` = one colony split three ways
- `{GCCTTATCAC, TAAACTAAGC}` = one colony split two ways

Collapsing those gives **exactly 8 colonies — the paper's number**, recovered from tape content
alone. The skeleton was rejoining split colonies; the ground-truth labels were wrong. A depth-5
clade spanning three groups cannot be homoplasy ($q^5\approx1.4\times10^{-9}$).

⇒ **Reusable QC finding: `ClonalBC` over-splits colonies, and tape-consensus similarity recovers the
true grouping.** This affects any per-clone analysis, including Park's own and our compatibility runs.

## Caveats

- $\hat m$ is capped at the clade size (correct: $m\le n_{\rm next}$) and floored at $s$. The cap
  is not binding at the largest nodes observed, so it is not distorting the tail.
- The $m$-identifiability result uses Mouse1's $\xi$; the other tables' $q$ agree to <3% so it
  should transfer, but it has not been recomputed per table.
- $q$ pooled over all **tapes** (per-tape breakdown not yet done — a per-tape $\xi$ would test
  the rest of row A3, and bears on the *cis*-preference question at §1312).
- Clone structure not yet joined in — all numbers above are per *table*, not per clone.

---

## Fig 3 redesign — characterising dropout (2026-09-01)

The old Fig 3 described *our handling* of `None`. The redesign measures the assay: is dropout
random, is it correlated within cells, is it correlated with edit depth? Informative dropout biases
any imputation or likelihood correction, so this is a prerequisite for the observation model that
row **A6** calls for.

`src/23_dropout_matrix.py` caches the (cell × tape) recovery / depth / termination matrices once
(`results/dropout_matrix_{arm}.npz`, gitignored); scripts 24–27 read that cache.

### Panel a — dropout is not a coin flip on entries

Null: every (cell, tape) instance recovered independently with one common $p$, so
$R_c=\sum_t Y_{ct}\sim\mathrm{Bin}(k,p)$, $k=166$, $\mathrm{Var}(R_c)=kp(1-p)$. The variance
inflation factor is the ratio of observed to null variance; giving each cell its own propensity
$\pi_c$ (mean $p$, variance $\sigma^2$) and applying the law of total variance gives
$\mathrm{Var}(R_c)=kp(1-p)+k(k-1)\sigma^2$, hence $\mathrm{VIF}=1+(k-1)\rho$ with
$\rho=\sigma^2/[p(1-p)]$ = the correlation between two dropout indicators **in the same cell**.

| arm | $\bar p$ | sd$(R_c)$ | null sd | VIF | $\rho_{\rm cell}$ | $\rho_{\rm tape}$ |
|---|---|---|---|---|---|---|
| Mouse1 | 0.606 | 29.9 | 6.30 | **22.6** | 0.131 | 0.250 |
| Mouse2 | 0.601 | 32.3 | 6.31 | 26.2 | 0.153 | 0.262 |
| Mouse3 | 0.559 | 35.5 | 6.40 | 30.8 | 0.180 | 0.202 |
| Initial | 0.749 | 19.0 | 5.58 | 11.5 | 0.064 | 0.231 |
| Subclone | 0.779 | 13.1 | 5.35 | 6.0 | 0.030 | 0.258 |

For Mouse1, $\sigma=0.177$: cells' recovery propensities scatter by ±18 points around 61%, where
the null allows none. Equivalently $P(\text{tape }t'\text{ missing})=39.4\%$ rises to **47.3%**
given another tape is missing in the same cell.

⚑ **$R_c$ is bimodal, not merely overdispersed** — a mode near 120 plus a flat shelf from the QC cut
at 20 up to ~100 holding **38%** of Mouse1 cells (Mouse2 37%, Mouse3 46%).

### The shelf is the whole per-cell effect, and it is a QC artefact

Re-measured at a common ≥100-tape cut, $\rho_{\rm cell}$ collapses in every arm and the mice become
the *most* homogeneous:

| | raw | at ≥100 |
|---|---|---|
| Mouse1 / 2 / 3 | 0.131 / 0.153 / 0.180 | **0.016 / 0.012 / 0.014** |
| Initial / Subclone | 0.064 / 0.030 | 0.021 / 0.024 |

Initial and Subclone were admitted at ≥100 tapes, the mice at ≥20 — so the other arms' QC deleted
exactly this population. This sharpens §D.4c's "part of the dropout asymmetry is a QC threshold"
to: **all of it is.** (Applying the filter ourselves drops 4.66% / 0.76% / 0.11% of cells,
reproducing §D.4c's 2.1% overall exactly.)

### The shelf and clonal-barcode loss are one phenomenon — informative selection

$P(\text{no ClonalBC}\mid R_c)$ is monotone in the mice and flat elsewhere:

| | 20–40 | 60–80 | 100–120 | 140+ |
|---|---|---|---|---|
| Mouse1 | **64%** | 52% | 44% | **17%** |
| Mouse3 | 52% | 47% | 29% | — |
| Initial | 3% | 3% | 8% | 3% |

A cell enters tree reconstruction only if it has a `ClonalBC`, so the analysis population is selected
on capture quality *twice* — the 44.9% of Mouse1 cells discarded for lacking a barcode are the
low-recovery cells. Censoring at the level of which cells exist, before per-entry dropout.

### ⚑⚑ The tape axis is large, reproducible, and estimable

Per-tape recovery rates span **0.006–0.962** in Mouse1 (deciles 0.195, 0.828), and
$\rho_{\rm tape}>\rho_{\rm cell}$ in every arm. Correlation of the 166 rates:

- between replicate libraries of the same population: **0.997** (Initial_1/2), **0.999** (Subclone_1/2)
- across arms — different animals, different preps: **0.76–0.96**

⇒ *which* tapes are badly recovered is a fixed property, not noise. A per-tape recovery probability
$\beta_t$ is measurable to three decimals and transfers across experiments, so it can enter the tip
emission $P(\mathrm{obs}\mid\mathrm{true})$ as a known constant (§1c.2 — the line where A6 enters).

### ⚑⚑ Dropout is informative on the tape axis, not the cell axis

Recovery vs mean depth given recovered (`src/26_dropout_depth.py`), Spearman over the 166 tapes,
and over cells:

| arm | per tape, raw | per tape, controlled† | worst decile → best decile | per cell, controlled‡ |
|---|---|---|---|---|
| Mouse1 | +0.344 | +0.281 | 2.92 → 4.95 sites (**+2.03**) | +0.046 |
| Mouse2 | +0.121 | +0.180 | 3.86 → 4.83 (+0.97) | −0.008 |
| Mouse3 | +0.326 | +0.02 ‡‡ | 3.55 → 5.06 (+1.51) | +0.077 |
| Initial | +0.324 | **+0.306** | 2.52 → 3.12 (+0.60) | +0.109 |
| Subclone | +0.139 | +0.147 | 3.86 → 4.45 (+0.58) | −0.155 |

† scored only on cells with $R_c\ge140$, so every tape is measured in a comparable cell population.
‡ mean depth over a fixed reference set of the 30 easiest tapes, identical for every cell.
‡‡ unpowered — Mouse3 has 30 such cells, Mouse2 76, Mouse1 296; Initial has 5,133 and is the one to
trust, where the control moves $\rho$ by 0.02.

**The tapes we mostly cannot see are the tapes that recorded least** — 2 sites of 6 in Mouse1, a
third of the recorder's dynamic range. Consistent with a shared per-locus latent (a closed
integration site is both poorly recovered and poorly edited): row **A9**'s mechanism, measured.

Two confounds excluded by design: the amplicon-length story (more editing → longer amplicon → worse
recovery) predicts the **opposite sign**; the coverage story (poorly captured cells lose terminal
sites, masquerading as unedited) predicts depth rising with $R_c$ within a tape, and it does not
(Mouse1: 4.892 → 4.986 across $R_c$ 20→140).

⚠ **Ceiling caveat.** The mouse arms sit at 4.8–5.1 of 6 sites filled (37–43% of tapes full), so
depth has little room to vary and their ≈0 per-cell coupling is partly attenuation. Initial, the
arm with dynamic range left (depth ~3), gives +0.109.

⚠ **Correction to a test proposed and then dropped.** An earlier draft argued that 0.02% internal
`None` proves trailing `None` is biology rather than dropout. That only excludes *random* site loss:
**terminal 3′-biased truncation produces no internal gaps at all**, and is exactly what the per-tape
depth deficit would look like. The flat per-cell depth profile is what actually excludes it.

⇒ **The two axes behave oppositely, in the favourable arrangement.** The informative axis (tape) is
the one that can be measured once and fixed; the axis that cannot be measured per cell is not
informative about editing, so it can be marginalised — which is what Felsenstein pruning already
does. This is the core argument of the redesigned figure.

**Still to measure:** whether dropout is lineage-correlated (do related cells lose the same tapes?
row **A9** proper). Mulberry & Stadler name this as their reason for punting on dropout (§1c.2) and
it is unmeasured in the literature. Attenuated here because shelf cells largely lack barcodes;
likely wants the simulator's null and its own figure.

### Panel b — the tape axis, and why the $\rho$'s are the comparable statistic

Panel a summed rows of the (cell × tape) matrix; panel b sums columns:
$\hat\beta_t=R_t/n$. Same null, roles swapped: $R_t\sim\mathrm{Bin}(n,p)$ so
$\mathrm{Var}(\hat\beta_t)=p(1-p)/n$ — for Mouse1 an sd of **0.0044**, i.e. the null puts all 166
tapes within ±1.3 points of 60.6%. Observed sd across tapes is **0.244**, 55× wider; rates run
**0.006 to 0.962** (deciles 0.195, 0.828).

⚠ **Do not compare the two margins' VIFs.** $\mathrm{VIF}=1+(m-1)\rho$ and $m$ is the number of
items summed — 166 tapes per cell, but 12,232 cells per tape — so $\mathrm{VIF_{tape}}=3{,}052$ vs
$\mathrm{VIF_{cell}}=22.6$ reflects the shape of the experiment, not the strength of the effect.
Only the $\rho$'s compare: **$\rho_{\rm tape}=0.250$ vs $\rho_{\rm cell}=0.131$ — the tape is the
larger axis.**

$\rho$ is a genuine correlation between two 0/1 entries sharing a unit. With $\pi_c$ the cell's
propensity (mean $p$, variance $\sigma^2$): $\mathrm{Var}(Y_{ct})=p(1-p)$ (a mixture of Bernoullis
is still Bernoulli marginally, so heterogeneity is invisible in one entry), while
$E[Y_{ct}Y_{ct'}]=E[\pi^2]=p^2+\sigma^2$ gives $\mathrm{Cov}=\sigma^2$ and hence
$\mathrm{corr}=\sigma^2/[p(1-p)]=\rho$. **Between-unit variance and within-unit covariance are the
same number.** Three readings of the magnitude:

1. $\rho=R^2$ of a one-way ANOVA of the indicator on that factor — 25% of the variance in "was this
   entry recovered?" is explained by which tape it is, 13% by which cell.
2. $P(\text{miss}\mid\text{miss in the same unit})=(1-p)+\rho p$: from 39.4% to **47.3%** (same
   cell) or **54.5%** (same tape).
3. $\sigma=\sqrt{\rho p(1-p)}$ = 0.177 (cells), 0.244 (tapes), in probability units.

**The spread is real by an analytic argument, not just by replication.**
$\mathrm{Var}_t(\hat\beta_t)=\mathrm{Var}(\beta_t)+E[\beta_t(1-\beta_t)]/n$, and the noise term is
**0.033% of the observed variance**, giving $\sigma_\beta=0.2441$ from 0.2441. On Initial the
naive spread (0.2081) and the covariance-based $\sigma_\beta$ from two independent libraries
(0.2081) agree to four digits — the check, not the argument.

Reliability measured anyway (`results/fig3b_tape_axis.json`), reported in text rather than plotted
since it was a foregone conclusion: Initial_1 vs Initial_2 $r=\mathbf{0.9968}$ against 0.9998
predicted from counting noise alone. Regressing out the global library-depth difference
($p$ 0.759 vs 0.741, slope 1.003) leaves a residual sd of 0.0168, of which counting noise supplies
0.0045 — so **1.6 points of library-specific wobble per tape against a 21-point real spread, i.e.
99.4% of the between-tape variance is a property of the tape.** Cross-arm $r$ vs Mouse1: Mouse3
0.959, Mouse2 0.875, Initial 0.846, Subclone 0.760 — largely transferable, not perfectly.

⇒ $\hat\beta_t$ has a signal-to-noise of 55 and enters the tip emission
$P(\mathrm{obs}\mid\mathrm{true})$ as a **known constant** — no prior, no extra parameter, no cost
to the pruning recursion. Estimate it per arm (replicates say one library suffices; cross-arm
correlations say do not import it from another animal). Mild circularity to state: $\beta_t$ comes
from the same data, but from the missingness marginal only and never from the tree — a plug-in
step, not double-counting.

⚠ **Plot caveat.** The four grey curves are each sorted on their own rates, so rank 40 is a
different tape in each arm: they show the *shape* is universal, not that the same tapes are bad
everywhere. That claim rests on the cross-arm correlations. Plotting the grey arms in Mouse1's
ordering was tried and is unreadable.

**Palette note (panel b).** Five arms now plot as five colours. Categorical slots 1,3,4,5,7 of the
`dataviz` reference palette — blue / aqua / yellow / magenta / violet, assigned Mouse1→Mouse3,
Pre-TX, Subclone, so Mouse 1 stays blue across panels. Orange (slot 2) is skipped: it is reserved
for the coin-flip null, and it fails the normal-vision floor against magenta ($\Delta E$ 12.9) and
yellow (13.7); red (slot 8) fails it worse (7.1). With orange unavailable to the series, **panel b's
null band is drawn neutral** — there it is reference furniture, not a competing series — so orange
never carries two meanings in the figure. Validated on the adjacent pairlist (the documented one for
line charts): CVD $\Delta E$ 9.1, normal-vision 19.6, lightness and chroma pass; aqua/yellow/magenta
fall below 3:1 on the light surface, so the relief rule applies and is met by the legend labels plus
the per-arm tables here.

⚠ The shipped `validate_palette.js` will not run on the cluster (node v10.24 cannot parse its ESM or
`??=`). Ported faithfully to Python — same thresholds, same Machado–Oliveira–Fernandes 2009
severity-1.0 matrices, same OKLab $\Delta E\times100$ — in the session scratchpad; re-port if needed.

⚠ The null band in panel b is **Mouse 1's** ($p=0.606$); each arm has its own $p$. The claim is the
band's *width* (±1.3 points, sd 0.0044–0.0053 across arms), not its position, and the annotation
says so.

### Panel c — the shelf is a QC choice, not a difference between arms

All five arms on panel a's axis, with the two admission thresholds the paper actually used marked:
**≥20** recovered tapes for Mouse1–3, **≥100** for Initial and Subclone. The mice carry a broad
plateau across 20–100 tapes; Pre-TX and Subclone are empty there, and their emptiness begins exactly
at their own cut.

| share of cells with 20–100 tapes | $\rho_{\rm cell}$ raw | at a common ≥100 cut |
|---|---|---|
| Mouse 1 **38%** · Mouse 2 37% · Mouse 3 46% | 0.131 / 0.153 / 0.180 | **0.016 / 0.012 / 0.014** |
| Pre-TX 4% · Subclone 1% | 0.064 / 0.030 | 0.021 / 0.024 |

**Diagnostic run first, and it could have overturned the claim.** If the shelf concentrated in
particular harvest sites it would be a dissection/prep batch effect. Using the same variance
machinery — $\rho_{\rm sample}$ = fraction of $\mathrm{Var}(1[20\le R_c<100])$ explained by sample
identity — gives **0.015 / 0.013 / 0.050** for Mice 1–3, pooled **0.022** over 12 harvest samples,
against $\rho_{\rm cell}$ of 0.13–0.18. Harvest site explains about a seventh of what cell identity
does, so the shelf is a per-cell property. (Mouse3's 0.050 is the largest and is driven by M3_LV,
$n=287$, 12.9% shelf against its littermates' 46–53%.)

⚠ **What this does not say.** The mouse analysis population *is* admitted at ≥20, so
$\rho_{\rm cell}=0.13$–0.18 is the real number any inference on those data must handle. The ≥100
comparison establishes what *kind* of thing it is — a data-quality population with no biological
content, hence marginalisable — as against the tape axis, which is informative and must be modelled.

### Panels d and e — dropout is informative on the tape axis, not the cell axis

Plotted on a **shared $y$-axis** (mean edit depth of recovered tapes, in sites of 6), because the
contrast is the argument and it disappears if they are scaled independently. Mouse 1.

**(d) per tape.** $\rho=+0.344$ raw, $+0.281$ scored only on high-capture cells ($R_c\ge140$), so
the confound that low-recovery tapes are seen mostly in good cells explains little. Decile means run
**2.92 → 4.95 sites**.

**⚑ The decile trend shows the coupling is a THRESHOLD, not a gradient** — a fact the correlation
coefficient hid. Depth climbs steeply to $\hat\beta_t\approx0.3$ and is flat above it:

| | tapes | mean depth | $\rho$ within |
|---|---|---|---|
| $\hat\beta_t < 0.3$ | 25 (15%) | **3.58** | — |
| $\hat\beta_t \ge 0.3$ | 141 (85%) | **4.86** | **+0.119** |

⇒ **the informative part of dropout is confined to a removable minority of tapes.** Drop the worst
15% and the coupling largely goes with them — structurally the same result as the homoplasy
concentration finding, and the same remedy. It also fits row **A9**'s mechanism better than a graded
effect would: a *subset* of integration sites is closed or silenced, and those loci are
simultaneously unreadable and unedited.

**(e) per cell.** $\rho=+0.046$, binned depth **4.88 → 4.98 sites** across $R_c$ 20→145 — flat.
Scored on a fixed reference set of the 30 easiest tapes, identical for every cell, without which a
low-$R_c$ cell would be graded on its easy tapes only.

⚠ **Ceiling caveat** (on the panel): Mouse 1 tapes average 4.9 of 6 sites, so depth has little room
to vary and (e)'s flatness is partly attenuation. Pre-TX, at ~3 of 6, is the arm with dynamic range,
and it gives $\rho=+0.109$ — still small beside (d).

⇒ **The axis that is informative is the one you can measure and, if need be, delete; the axis you
cannot measure per cell carries no information about editing and can be marginalised.** That is the
figure's conclusion and the argument for the likelihood route.

---

## Row A9 — is tape loss heritable? (2026-09-02)

**Question.** Mulberry & Stadler punt on dropout because it may be non-random: *"if DNA Tapes are
simultaneously lost for groups of related cells"* (§1c.2). Nobody has measured it. Do related cells
lose the **same** tapes?

**Why the answer changes the model.** Fig 3's emission $P(\mathrm{obs}\mid\mathrm{true})$ built from
$\alpha_c,\beta_z$ is valid only if missingness is conditionally independent across tips given the
tree. If loss is heritable: (i) the tip factorisation that pruning relies on fails in the
*observation* layer, and computes a wrong number with nothing to flag it; (ii) we are discarding a
Dollo character — and in a chromosomally unstable cancer, copy-number loss at an integration site is
the obvious mechanism; (iii) clustered *cell* loss breaks the birth–death uniform-sampling prior.
The remedies differ: (i)/(ii) want a per-tape irreversible loss process — **one absorbing state in
the pruning recursion, a constant factor, not an explosion**; only (iii) is a genuine SBI argument.

### The statistic

For cell $c$ and feature $f$, $r_{cf}=Y_{cf}-\hat p_{cf}$ is the **surprise** — outcome minus what
the cell and the feature alone predict, with $\hat p=\sigma(\alpha_c+\beta_z)$ so that a clone of
uniformly poor cells leaves no residual and only **tape-specific** structure survives. Products of
two cells' surprises, averaged over pairs sharing a group and normalised by
$v=\hat p(1-\hat p)$, give an intraclass correlation. Pairs are never enumerated:

$$\sum_{c \neq c' \in g} r_c r_{c'} = \Big(\sum_{c\in g} r_c\Big)^2 - \sum_{c\in g} r_c^2$$

so two running totals per (group, feature) deliver every pair — $O(nk)$ for all 166 tapes. Computed
at clone and sample level, so between-clone comes free (same-sample pairs minus same-clone pairs).

**T0 is a calibration curve, not a number.** $\rho$ depends on a feature's marginal frequency, so the
control is the *same statistic on features known to be heritable* — "has tape $z$ reached site $L$",
$L=2\ldots6$, whose marginals sweep 0.97 down to 0.35 and so overlap the missingness marginals.
Missingness is read against that curve; the within-sample permutation supplies the zero line.

⚠ **$\rho_{\rm between}$ is NOT the comparator.** The two-way fit has no per-sample term, so the mean
residual within a harvest sample is not exactly zero; squared over group sums of thousands of cells
that offset dominates and inflates the between term. The **permutation null** carries the same
offset by construction and is the sound contrast. Retained in the output for the record only.

### T1 — clone level (`src/34_lineage_rho.py`)

Excess = observed $\rho_{\rm within}$ minus the within-sample permutation null (200 draws).

| arm | missing marg | excess | matched T0 | T0 marg | T0 excess | ratio |
|---|---|---|---|---|---|---|
| Mouse 1 | 0.363 | **+0.1662** | depth≥6 | 0.377 | +0.3191 | 0.52 |
| Mouse 3 | 0.400 | **+0.1634** | depth≥6 | 0.400 | +0.5363 | 0.30 |
| Pre-TX | 0.232 | **+0.1608** | depth≥4 | 0.346 | +0.1208 | **1.33** |
| Subclone | 0.218 | **+0.2630** | depth≥6 | 0.272 | +0.2842 | 0.93 |
| Mouse 2 | 0.371 | +0.0092 | depth≥6 | 0.350 | +0.0290 | 0.32 |

**Heritable in every arm**, at 30–130% of a known-heritable character of matched frequency. In
Pre-TX missingness is *more* clone-clustered than the edit features.

⚠ **Mouse 2's small value is a power artefact, not biology.** One clone holds 3,387 of its 5,382
cells, so a permuted group of that size drawn from the same samples is nearly the real clone and the
test has no leverage there — its null sits at +0.030 where the other arms' sit at ~0.002. Mouse 2 is
the floor, not the signal.

Most of the *raw* within-clone agreement is tape marginals, not lineage: Mouse 2 $\rho_{\rm within}$
falls 0.406 (grand mean) → 0.032 (tape-centred) → 0.040 (two-way).

### T3 — is capture quality itself clone-clustered?

$\rho_{\rm clone}$ on $R_c$, the birth–death uniform-sampling question, kept as a finding in its own
right rather than only as a nuisance to condition away:

| arm | within clone | between clones, same sample |
|---|---|---|
| Pre-TX | **+0.1534** | +0.0047 |
| Mouse 1 | **+0.0899** | +0.0059 |
| Subclone | +0.0383 (screened +0.0446) | −0.0192 |
| Mouse 3 | +0.0323 | +0.0212 |
| Mouse 2 | +0.0036 | −0.0024 |

Related cells do vary together in overall capture, strongly in Pre-TX and Mouse 1 — the uniform
sampling assumption failing directly. ⚠ Not yet callable as biological: clone-mates could share a
capture level through co-encapsulation or sub-lane structure the `Sample` label does not resolve.
**Needs a control before it goes in a figure.**

### T2 — finer than clonal: the relatedness gradient (`src/35_lineage_depth.py`)

Clone is coarse (Subclone's median clone is 997 cells) and, in Subclone, colony ≡ culture batch. Both
problems resolve by looking below the clone. A subclade = cells sharing the depth-$d$ prefix of one
**anchor** tape; missingness is measured on every tape **except** the anchor (cross-tape only, per
§D.4 fact 1), averaged over ~160 anchors.

**The null is analytic.** Under random assignment of a clone's cells to subgroups of the observed
sizes, every unordered pair is equally likely to land together, so

$$E\Big[\sum_{c\neq c'\in g} r_c r_{c'}\Big] = \frac{m(m-1)}{n(n-1)}\sum_{c\neq c'\in C} r_c r_{c'}$$

— the permutation mean is a size-weighted rescaling of the clone's own pair sum, no simulation.
**Verified** against explicit within-clone permutations (8 anchors × 40 draws): ratios 0.9992–1.0006,
**mean 1.0000**, each within its permutation sd.

Excess = agreement *beyond* the clone. A flat batch effect predicts zero at every depth.

| arm | d1 | d2 | d3 | d4 | rise | d4 ÷ T0 |
|---|---|---|---|---|---|---|
| Pre-TX | +0.0628 | +0.0973 | +0.1320 | **+0.1547** | 2.5× | 0.58 |
| Mouse 3 | +0.0192 | +0.0369 | +0.0511 | **+0.0587** | 3.1× | 0.37 |
| Mouse 1 | +0.0072 | +0.0105 | +0.0161 | **+0.0195** | 2.7× | 0.39 |
| Subclone (screened) | +0.0070 | +0.0079 | +0.0091 | **+0.0122** | 1.7× | 0.63 |
| Mouse 2 | +0.0004 | +0.0011 | +0.0026 | **+0.0038** | 9.5× | 0.47 |

⚑⚑ **Monotone in five of five arms** (six of six with Subclone unscreened). Closer relatives agree
more about which tapes they have lost, every time — the Dollo prediction, not the batch prediction.

⚑ **The depth test did the job it was added for.** Subclone has the *largest* clone-level signal
(+0.263) and the *flattest* gradient, with sub-clone structure only ~5% of its clone-level number;
Pre-TX's sub-clone excess (+0.155) nearly equals its clone-level one (+0.161). So Subclone's headline
is substantially colony-≡-well batch, while Pre-TX's and the mice's is genuinely tree-structured. The
colony confound is now bounded rather than argued about. (Clone-level and sub-clone excesses use
different nulls and cell subsets, so the ratio is indicative, not exact.)

### ⚑ An unexpected result about the collision screen

At clone level the screen barely matters and moves things the predicted way (Subclone +0.263 →
+0.273; misassigned cells dilute). **At subclade level it halves the signal** (+0.0189 → +0.0122 at
d4). In hindsight this is right: a misassigned cell carries a foreign edit prefix *and* foreign
dropout, so it forms a spurious subgroup that is internally coherent — it *manufactures* subclade
agreement rather than diluting it. ⇒ **the screen matters more for the fine test than the coarse one,
and in the opposite direction.** Use screened numbers for the depth sweep, unscreened for clone level.
Mouse 2 flags 11/5,382 cells (0.20%); Subclone 1,304/38,652 (3.37%), against `src/18`'s 3.77% on
pooled 200-cell groups.

### Where this leaves the model

Tape loss is heritable, and heritable **below** the clone — the property that makes it a Dollo
character on the tree rather than a per-clone nuisance. That points at the tractable remedy: a
per-tape irreversible loss process, one absorbing state integrated along branches in the pruning
recursion. **A9 does not, on this evidence, force SBI** — contrary to the §1c.3 table. The
uniform-sampling half of A9 (clustered *cell* loss) is untouched by that fix and remains the genuine
SBI argument.

**Not yet done:** figure panels; a control for the Pre-TX capture-ICC; depth sweep beyond d=4.

### Dollo or a heritable propensity? (`src/36_dollo_test.py`, 2026-09-03)

The depth gradient shows concordance follows the topology, but **two mechanisms predict a gradient**
and they need different models: a single irreversible **loss** in an ancestor (⇒ one absorbing state
per tape in the pruning, a constant factor) versus a heritable **propensity** — a lineage with a
lower recovery rate at that locus (⇒ a per-lineage rate multiplier, a latent field over the tree, a
real SBI case). The discriminator is the distribution of

$$a = \frac{\text{cells of subclade } S \text{ missing tape } z}{|S|},\qquad z \neq \text{anchor},$$

over subclades whose **own clone still carries the tape** (so the absence is a loss *below* the
clone, not one inherited from it).

⚠⚠ **The first version of this test used the wrong null and overstated the result.** Comparing $a$
against Bernoulli draws from $\hat p=\sigma(\alpha_c+\beta_z)$ assumes the additive-logit fit is
correct; **any** cell×tape interaction, lineage-related or not, makes $a$ U-shaped against it — and
the enrichment we saw at $a\approx0$ was itself evidence the fit is too smooth. The honest
comparator shuffles **whole cell rows within a clone**, preserving each cell's entire missingness
profile and each tape's marginal while destroying only the alignment between cells and subclades.
The permutation null is far higher than the Bernoulli one (Mouse1, depth 3, 10–20 cells: 20,065 vs
12,057), so much of the first-pass signal was cell×tape structure. **All numbers below are against
the permutation null.**

**Complete absence in a subclade whose clone still carries the tape**, observed ÷ permuted:

| arm | 10–20 cells | 20–50 | 50–200 | 200+ | well-recovered tapes ($\beta\ge0.5$) |
|---|---|---|---|---|---|
| Mouse 1 | 1.7–2.2× | 2.5–3.7× | 3.2–5.5× | 3.5–10.6× | 4.5–25× |
| Mouse 3 | 1.8–2.3× | 2.2–2.5× | **6.3–11.8×** | — | 3.5–9.8× |
| Mouse 2 | 1.4–1.8× | 1.6–3.2× | 2.0–5.0× | 1.5–4.8× | 2.4–7.1× |
| Pre-TX | 1.9–2.3× | 2.0–3.0× | 2.9–3.8× | — | 3.0–8.5× |
| Subclone | **0.9–1.6×** | 1.1–1.8× | 1.1–1.9× | 2.0–2.2× | 0.8–2.2× |

⚑ **Subclone is the weak arm here** — at shallow depths and small subclades it sits *at or below* the
null. Consistent with the depth sweep: its large clone-level signal is substantially colony-≡-well
batch, and genuine sub-clone loss is modest.

**⚑⚑ The shape is the real evidence, and it is consistent in 5/5 arms.** Mass does not shift smoothly
toward higher $a$; it moves out of "nearly all missing" and into "**exactly** all missing":

| arm | $a\in[0.85,0.90)$ | $[0.90,0.95)$ | $[0.95,1.00]$ | middle $[0.25,0.85)$ |
|---|---|---|---|---|
| Mouse 1 | 0.70× | 0.76× | **1.46×** | 0.96× |
| Mouse 2 | 0.70× | 0.80× | **1.34×** | 0.92× |
| Mouse 3 | 0.78× | 0.91× | **1.36×** | 1.01× |
| Pre-TX | 0.91× | 1.00× | **1.50×** | 0.99× |
| Subclone | 0.70× | 0.90× | **1.50×** | 0.96× |

A graded propensity broadens the middle and fills the near-complete bins; a discrete irreversible
loss empties them into the complete bin. **The sharper the criterion, the larger the excess** —
$a\ge0.95$ gives 1.3–1.5×, $a=1$ exactly gives 2–4× in the same cells — which is itself the
signature of a discrete state rather than a rate.

⇒ **Reading: dropout is a MIXTURE** — a technical component well described by $\alpha_c\beta_z$, plus
a heritable discrete loss. The modelling consequence is the cheap one: a per-tape irreversible loss
process alongside the editing process, one absorbing state in the pruning. Effect sizes over the
honest null are modest (1.2–1.5× on shape, 1.4–4× on counts), so this is a real but not overwhelming
component — it should be *modelled*, not treated as the dominant term.

⚠ **A consequence for Fig 3b that needs following up.** If part of per-tape missingness is heritable
locus loss rather than per-observation failure, then $\beta_z$ is not purely a technical constant —
treating it as a fixed emission probability would model as noise something that is a tree-structured
state. The cross-library reproducibility ($r=0.997$) does not settle this, because all arms descend
from one engineered line, so an ancestral loss is shared by construction. **Open.**

---

## Methods, in full: the within-clone correlation $\hat\rho_{\rm within}$

*Written out completely (2026-09-02) so it can be re-read before presenting. Nothing here is new
analysis; it is the derivation behind every number in the A9 section above.*

### 1. The data

One binary number per cell–tape pair,

$$Y_{cz} = 1 \ \text{if tape } z \text{ is entirely unrecovered in cell } c,\qquad 0 \text{ otherwise},$$

i.e. a 6,737 × 166 matrix for Mouse 1. (Flipping the 0/1 convention on both members of a pair leaves
every correlation below unchanged; only consistency matters.)

### 2. Why a model is needed at all

We want to know whether two cells of a clone agree **more than expected**, so "expected" must
already contain the two nuisances Fig 3 measured: some cells are badly captured ($\rho_{\rm
cell}=0.13$), some tapes are badly recovered ($\rho_{\rm tape}=0.25$). Otherwise sisters look
concordant merely for sharing a bad tape.

⚑ **And each $(c,z)$ pair is observed exactly once.** There is no repeat measurement from which to
estimate $E[Y_{cz}]$ locally. A model is the only device that can borrow strength from the whole row
and the whole column to produce an expectation for a single unrepeated observation. That is the sole
reason $\alpha$ and $\beta$ exist.

### 3. Why not just use the observed fractions?

Let $a_c$ = fraction of tapes missing in cell $c$, $b_z$ = fraction of cells missing tape $z$,
$\bar p$ = overall. A combining rule must (i) reproduce the observed margins and (ii) return a
number in $[0,1]$.

**Multiplicative**, $\hat p_{cz}=a_cb_z/\bar p$, satisfies (i) —
$\sum_z a_cb_z/\bar p=(a_c/\bar p)\,K\bar p=Ka_c$ ✓ — but **fails (ii) on our real data**: Mouse 1
has $\bar p=0.394$, a real cell with 20/166 tapes recovered has $a_c=0.88$, and a real tape recovered
in 0.6% of cells has $b_z=0.994$, giving $\hat p = 0.88\times0.994/0.394 = \mathbf{2.22}$. Then
$v=\hat p(1-\hat p)$ is negative and the whole construction collapses.

**Additive on the probability scale**, $\hat p = a_c+b_z-\bar p$, also matches the margins and also
fails: $0.88+0.994-0.394=1.48$.

Both linear rules break precisely at the badly recovered cells and tapes — which is where heritable
loss would live.

### 4. The model: $\alpha$, $\beta$, and $\hat p$

Add on the **log-odds** scale, which cannot leave $(0,1)$ and saturates correctly:

$$\mathrm{logit}(p)=\log\frac{p}{1-p},\qquad \mathrm{logit}^{-1}(x)=\frac{1}{1+e^{-x}},\qquad
\hat p_{cz}=\mathrm{logit}^{-1}(\alpha_c+\beta_z).$$

* $\alpha_c$ — cell $c$'s propensity to be missing, in log-odds. Large = poorly captured cell.
* $\beta_z$ — tape $z$'s propensity to be missing, in log-odds. Large = badly recovered tape.

Worked: $\alpha_c=0.4,\ \beta_z=0 \Rightarrow \hat p=0.599$; same cell on a good tape
$\beta_z=-2 \Rightarrow \hat p=0.168$. The additive-log-odds form means a bad cell multiplies the
*odds* of missingness by the same factor on every tape — exactly what a capture failure does.

⚠ **Notation.** Earlier drafts used $\sigma$ for both the logistic function and a variance. Here
$\mathrm{logit}^{-1}$ is the link; $\sigma^2$ only ever means a variance.

### 5. How $\alpha,\beta$ are fitted — and why Newton is forced, not chosen

The maximum-likelihood score equations for this model, with row and column indicators as covariates,
are $\sum (Y-\hat p)\times(\text{covariate})=0$, i.e.

$$\sum_z\big(Y_{cz}-\hat p_{cz}\big)=0\ \ \forall c,\qquad \sum_c\big(Y_{cz}-\hat p_{cz}\big)=0\ \ \forall z .$$

**The ML fit reproduces every observed row and column total** — precisely the property wanted from
the fractions in §3. The difference is that $\mathrm{logit}^{-1}$ is nonlinear, so the $\alpha,\beta$
achieving those totals are *not* the fractions themselves and must be solved for: alternate
one-dimensional Newton steps on $\alpha$ given $\beta$ and vice versa. Adding a constant to all
$\alpha$ and subtracting from all $\beta$ changes nothing, so fix $\sum_c\alpha_c=0$.

⇒ **The fractions are what the fit matches; $\alpha,\beta$ are whatever numbers land on them.**
(Standard object: the Rasch model of item-response theory — "person ability + item difficulty" on the
logit scale; the margin-matching is the logistic analogue of iterative proportional fitting.)

⚑ **This is the crux of the test.** Because the fit reproduces every row and column total, any
purely per-cell or purely per-tape structure is gone. Whatever survives in the residuals is *by
construction* an interaction — a particular cell missing a particular tape beyond what its own
quality and that tape's quality imply.

### 6. $r$ and $v$ are not extra concepts

$r_{cz}=Y_{cz}-\hat p_{cz}$ is just "$Y$ minus its expectation", which the definition of correlation
already requires; $\hat p_{cz}$ *is* the model's $E[Y_{cz}]$. And $Y_{cz}$ is a coin flip with
probability $\hat p_{cz}$, so its variance is $v_{cz}=\hat p_{cz}(1-\hat p_{cz})$ — the yardstick we
divide by. Written with neither symbol:

$$\hat\rho_{\rm within}=\frac{\sum_C\sum_{c\neq c'\in C}\big(Y_{cz}-\hat p_{cz}\big)\big(Y_{c'z}-\hat p_{c'z}\big)}
{\sum_C\sum_{c\neq c'\in C}\sqrt{\hat p_{cz}(1-\hat p_{cz})}\ \sqrt{\hat p_{c'z}(1-\hat p_{c'z})}}$$

**The modelling choice is *which* expectation to centre on**, and it is the whole design:

| centre on | question answered | Mouse 2 |
|---|---|---|
| grand mean $\bar p$ | do sisters agree? | +0.41 — dominated by "some tapes are bad" |
| tape marginal $\beta_z$ | …beyond tape quality? | +0.03 |
| $\alpha_c+\beta_z$ | …beyond cell **and** tape quality? (**partial correlation**) | +0.04 |

### 7. Collapsing the pair sums

For any list $x_1\ldots x_m$, $\big(\sum_i x_i\big)^2=\sum_i x_i^2+\sum_{i\neq j}x_ix_j$, hence
$\sum_{i\neq j}x_ix_j=\big(\sum_i x_i\big)^2-\sum_i x_i^2$. Applying it to the numerator with $x=r$
and the denominator with $x=\sqrt v$ (so $(\sqrt v)^2=v$):

$$\hat\rho_{\rm within}=\frac{\sum_C\Big[\big(\sum_{c\in C}r_{cz}\big)^2-\sum_{c\in C}r_{cz}^2\Big]}
{\sum_C\Big[\big(\sum_{c\in C}\sqrt{v_{cz}}\big)^2-\sum_{c\in C}v_{cz}\Big]}$$

Four running totals per (clone, tape); every pair accounted for; $O(nk)$ for all 166 tapes. Ordered
pairs throughout, so the factor 2 cancels.

**Worked example — one clone of three cells, one tape.** $\hat p=(0.6,0.6,0.3)$, so
$v=(0.24,0.24,0.21)$, $\sqrt v=(0.490,0.490,0.458)$, and denominator
$(1.438)^2-0.69=1.378$ either way.

*All three missing*, $Y=(1,1,1)$, $r=(0.4,0.4,0.7)$: numerator $(1.5)^2-0.81=1.44$; brute force
$2[0.16+0.28+0.28]=1.44$ ✓; $\hat\rho=+1.05$.
*Mixed*, $Y=(1,0,1)$, $r=(0.4,-0.6,0.7)$: numerator $(0.5)^2-1.01=-0.76$; $\hat\rho=-0.55$.

Agreement positive, disagreement negative, weighted by how surprising it was. A single extreme clone
can exceed 1 — this is a ratio-of-moments estimator with a fitted mean; the $\rho\le1$ bound applies
to the pooled quantity.

### 8. Why the marginal frequency must be matched

With $\pi_C$ the clone's own propensity (mean $p$, variance $\sigma^2$), $\pi^2\le\pi$ gives
$\sigma^2=E[\pi^2]-p^2\le p-p^2=p(1-p)$, so $\rho=\sigma^2/[p(1-p)]\le1$. The useful consequence is
at the extremes: a feature present in 97% of cells has $p(1-p)=0.029$ — almost no room for clones to
differ — against 0.24 at $p=0.4$. Visible in our own data: the `depth≥2` control, marginal 0.97,
gives the smallest excess in **every** arm (Mouse 1 +0.219 vs +0.319 for `depth≥6`; Mouse 3 +0.084 vs
+0.536). ⇒ controls must be compared at matched marginal, which is why panel a is a curve.

### 9. Why the permutation null is subtracted rather than assumed zero

$$\text{excess}(f)=\hat\rho_{\rm within}(f)-E_\pi\big[\hat\rho_{\rm within}(f)\big]$$

with $\pi$ shuffling cells **within harvest sample**, preserving each clone's size and its sample
composition. One expects $E_\pi\approx0$, and in three arms it is — but not in all:

| arm | worst between-sample per-tape $r$ | permutation null |
|---|---|---|
| Subclone | 0.999 | **−0.0000** |
| Pre-TX | 0.997 | **−0.0000** |
| Mouse 1 | 0.900 | +0.0021 |
| Mouse 3 | 0.867 | +0.0216 |
| Mouse 2 | 0.815 | +0.0303 |

⚑ **Perfectly monotone in how much the samples disagree about tape quality.** The fit has no
per-sample term, so where samples differ in their $\beta_z$ profile the residuals are not mean-zero
within a sample; squared over group sums of thousands of cells, that offset appears even for random
groups. The permutation carries the identical offset by construction, so subtracting it makes the
estimator robust to a misspecification we can *measure*. Assuming zero would have inflated Mouse 2 by
**4.3×** (+0.0395 vs +0.0092) and Mouse 3 by 13%.

⇒ The null value is a **diagnostic**, not a formality, and should be reported alongside the excess.

### Per-tape permutation nulls: p-values, z-scores, and what they are worth (`src/37_pertape_null.py`)

$\hat\rho_{\rm within}$ is computed per tape, and the per-tape values are wildly heterogeneous, so
the pooled number hides the structure a Dollo mechanism predicts. Script 37 keeps the whole
per-feature vector on every permutation draw (B = 5,000 Mouse1; 2,000 elsewhere), giving each tape

$$z = \frac{\rho_{\rm obs}-\overline{\rho_{\rm perm}}}{\mathrm{sd}(\rho_{\rm perm})},\qquad
p = \frac{1+\#\{b:\rho_{\rm perm}^{(b)}\ge\rho_{\rm obs}\}}{B+1}$$

one-sided upper (nothing predicts sisters agreeing *less* than strangers), with BH-FDR across the
166 tapes. `depth>=6` is run identically as a known-heritable positive control.

**✅ The null calibration check passes cleanly.** Tapes significantly *negative*: **0 out of 166, in
all five arms, for both feature families** — ten independent opportunities for a misspecified null to
produce spurious anti-concordance, and it produced none.

| arm | missing: sig. at FDR 5% | median $z$ | control: sig. | control median $z$ | null skew |
|---|---|---|---|---|---|
| Pre-TX | 166/166 | +100.0 | 165/166 | +54.1 | +0.14 |
| Subclone | 166/166 | +422.2 | 165/166 | +2063.5 | +1.49 |
| Mouse 1 | 165/166 | +47.8 | 160/166 | +95.8 | +1.58 |
| Mouse 3 | 154/166 | +12.9 | 158/166 | +43.9 | +0.75 |
| Mouse 2 | 128/166 | +7.8 | 147/166 | +19.3 | +0.32 |

⚠⚠ **Significance is saturated and must not be the headline.** $\rho$ pools over millions of
within-clone pairs — Mouse 1's largest clone alone contributes $1607\times1606\approx2.6$M ordered
pairs — so the permutation sd is minute and the test detects arbitrarily small departures. The
p-value answers *"is there any detectable excess"* (yes, on 77–100% of tapes) and not *"is it
large"*. **Effect size is the informative quantity.** Relatedly, the per-tape null is right-skewed
(median skew +0.14 to +1.58), so $z$ is an **ordering statistic**, not something to convert to a
Gaussian p-value.

⇒ The reading is **pervasive but unequal**: with hundreds of clones per arm, even a low per-locus
loss rate means nearly every tape is lost in *some* lineage and so shows a detectable excess, while a
minority of loci carry most of the magnitude.

### ⚑⚑ What the per-tape nulls actually buy: dropout heritability is *more concentrated* than edit heritability

With per-tape nulls the excess $\rho_{\rm obs}-\overline{\rho_{\rm perm}}$ is comparable across
tapes, and the concentration can be compared against the control measured identically:

| arm | missing: median excess | max | **top 10% of tapes hold** | top 25% | control top 10% |
|---|---|---|---|---|---|
| Mouse 1 | +0.049 | +2.08 | **58%** | 77% | 23% |
| Mouse 3 | +0.058 | +1.00 | **43%** | 71% | 18% |
| Subclone | +0.034 | +1.18 | **44%** | 81% | 24% |
| Mouse 2 | +0.002 | +0.10 | 57% | 82% | 59% |
| Pre-TX | +0.146 | +0.59 | 20% | 41% | 16% |

**In three of five arms dropout heritability is 2–2.5× more concentrated than edit-depth
heritability measured the same way.** That contrast is the point: edit depth accumulates along
lineages on *every* tape alike, so its heritability is spread evenly (control top decile 16–24%);
dropout heritability piles into a minority of loci (top decile 43–58%). Discrete loss events at
*particular* integration sites predict exactly that asymmetry, and a diffuse technical effect does
not. Reported as a **comparison** rather than an absolute, which is what makes it defensible.

Mouse 1's spread is 100-fold across tapes: deciles +0.011, +0.019, +0.049 (median), +0.103, +0.258,
95th +0.534, 99th +1.364, max +2.078.

⚠ Pre-TX is again the exception — diffuse (20% vs 16% control), as it was in the raw per-tape
distribution. Real arm difference, not yet explained.

**Link back to Fig 3, now null-controlled.** Spearman(tape recovery rate, per-tape excess) is
negative in every arm: Subclone −0.70, Mouse 3 −0.44, Mouse 2 −0.37, Mouse 1 −0.33, Pre-TX −0.11.
**Poorly recovered tapes are the heritably lost ones** — the mechanism behind Fig 3's "informative
dropout sits in a removable 15% of tapes", and direct evidence bearing on whether $\beta_z$ is a
technical constant. ⚠ Treat as suggestive: subtracting the per-tape null removes the baseline but
not the residual marginal dependence of the achievable range, and recovery rate *is* the marginal.

### ⚑⚑ The pooled $\hat\rho$ was hiding the signal, not creating it (`src/38_pool_robustness.py`)

$\hat\rho=\sum_C \mathrm{num}_C/\sum_C \mathrm{den}_C$ is a **ratio of sums**, so each clone enters
weighted by its pair count $n_C(n_C-1)$. Share of the pooled weight held by the single largest clone:

| arm | clones | largest clone | its share of the weight |
|---|---|---|---|
| Mouse 2 | 216 | 3,387 cells | **98.9%** |
| Mouse 1 | 295 | 1,607 | **84.0%** |
| Subclone | 15 | 10,996 | 47.6% (top 3: 83.9%) |
| Mouse 3 | 149 | 210 | 28.8% |
| Pre-TX | 2,946 | 127 | **1.6%** |

So "Mouse 2's $\rho$" was essentially *one clone's* $\rho$, and "5/5 arms replicate" overstated the
independence of the replication. Decomposing per clone:

| arm | pooled | drop largest clone | **equal weight per clone** | clones ≥5 cells, % positive |
|---|---|---|---|---|
| Mouse 2 | +0.009 | **+0.326** | **+0.246** | 78, 99% |
| Mouse 1 | +0.166 | +0.226 | **+0.222** | 137, 100% |
| Mouse 3 | +0.163 | +0.148 | **+0.237** | 67, 100% |
| Pre-TX | +0.161 | +0.163 | **+0.214** | 1,685, 100% |
| Subclone | +0.263 | +0.342 | **+0.353** | 13, 100% |

⚑ **Mouse 2's anomaly is fully explained** — its dominant clone has near-zero excess and carries
98.9% of the weight. Drop it and Mouse 2 is an ordinary arm. Equal-weighted, the five arms collapse
to a strikingly stable **+0.21 to +0.35**, where pooled they ranged 0.009–0.263.

⇒ Two consequences. The replication claim is far stronger than "5/5 arms": it is **~1,980 individual
clones, essentially all positive**. And **large clones are the wrong unit** — a clone spanning much
of the tree dilutes within-clone concordance, since a loss partway down does not make the whole clone
concordant. That is the same conclusion the depth sweep reached from the other direction. **Use the
equal-clone-weighted estimator from here on.**

### Panels a and b (`src/39_fig4ab.py`)

**`fig4a_heritable.png` — related cells lose the same tapes.** $x$ = the feature's marginal
frequency, $y$ = excess over the permutation null, equal weight per clone. Each arm contributes a
control *curve* (the five `depth≥L` features, $L=2\ldots6$, marginals sweeping 0.97→0.35) and one
filled diamond for missingness, with a dashed connector to the control point at the nearest
frequency. The curve is necessary because $\rho\le\sigma^2/[p(1-p)]$ bounds what is achievable at
extreme marginals — visible in the plot as every arm's control collapsing near $p=0.97$.

Read-off: **every diamond sits far above the null and inside its arm's control band** —
Mouse 1 40%, Mouse 3 39%, Mouse 2 31%, Subclone 82%, **Pre-TX 126%** of the matched control.

**`fig4b_perclone.png` — …in essentially every clone, not one big one.** One point per clone of ≥5
cells, sized by cell count, colour = missingness and grey = control, medians barred, and **each arm's
largest clone ringed**. The rings make the pooling result visible: Mouse 2's 3,387-cell clone sits
*on the null line* while its 77 smaller clones sit well above it.

⚠ Caveat carried on the panel: the control here is `depth≥6`, whose marginal matches missingness in
the mice (0.35–0.40 vs 0.36–0.40) but **not in Pre-TX** (0.05 vs 0.23). Pre-TX's 126% uses the
nearest available control (`depth≥4`) in panel a; its per-clone comparison in panel b is against
`depth≥6` and so is not frequency-matched.

---

## Panels c/d groundwork: an event catalogue (2026-09-03)

### ✅ Construct and readout verified from the paper — the mechanism is not speculative

`refs/metastasis_lineage_recording.pdf`, Methods p39–42:

- The recorder is **one piggyBac cassette, `PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC`** — so
  **pegRNA and tape are co-integrated in Park too**, as in the Typewriter and mouse-embryo
  lineage configurations. Prime editor separate (`LSL-PEmax-P2A-mClover3`, Cre-activatable).
- **The readout is RNA.** Cells sorted mRFP⁺, 10x 3′ v4, *"Feature cDNA Primer 3 replaced the
  standard cDNA primers so that the Read 2N primer was included and **TAPE cDNA co-amplified**"*,
  split by amplicon size into transcriptome / TAPE+TargetBC / ClonalBC libraries.

⇒ **A tape is recovered only if its integration is transcribed.** So epigenetic silencing of an
integration (i) removes the tape from the readout and (ii) removes the co-integrated pegRNA's symbol
from the cell's writing pool — both heritable. The mouse Typewriter paper already attributes low
per-tape recovery to *"epigenetic silencing of a subset of circTAPE-encoding integrations"*, and
chromatin-dependence of prime editing is measured (Li et al. 2024, *Cell* 187:2411). What is new
here is that it is **lineage-resolved**.

Three consequences:
- **Fig 3b is re-explained.** Per-tape recovery (0.006–0.96, $r=0.997$ across libraries) is the
  **integration site's expression level**, not primer efficiency — which is why it transfers across
  arms (one founder line, same insertion sites).
- **Fig 3d is re-explained.** Poorly recovered tapes having ~2 fewer filled sites is one chromatin
  state suppressing transcription *and* editing at the same locus — row A9's shared latent cause.
- **No selection confound**: mRFP comes from all 166 integrations, so silencing one leaves the cell
  mRFP-high and in the data.

⚑ 166 integrations each drawing `NNNN` from 256 predicts $256(1-e^{-166/256})=122$ distinct symbols;
measured alphabet is 97–129. So the tape↔symbol map is **1:1 by construction**, ~1.4 integrations per
symbol. Why the pairing is unknown: the library was cloned in two independent degenerate steps
(`NNNN` then `N10`), randomly paired, and **no read spans both** — the pegRNA is a separate U6
transcript at the far 5′ end. Recovering it needs long-range sequencing of the intact integration,
which Park did not do (the mouse paper did, for 10 of 11).

### The statistic: a log-likelihood ratio for a Dollo loss

For clade $S$ ($m$ cells, $k$ missing tape $z$), against the fitted per-cell $\tilde p_{cz}$:

$$\Lambda(S,z) = \sum_{\text{missing}} \log\frac{1-\varepsilon}{\tilde p_{cz}} + \sum_{\text{present}} \log\frac{\varepsilon}{1-\tilde p_{cz}}$$

Two properties, both wanted: **capture-independence is built in** (a missing cell earns
$-\log\tilde p$, so a well-captured cell missing a reliable tape dominates: +2.98 nats at
$\tilde p=0.05$ vs +0.09 at 0.90), and **it demands all-or-none** (a present cell costs up to
$-4.55$). ⚑ The penalty is *larger* for a good cell ($-4.55$ at $\tilde p=0.05$) than a bad one
($-2.30$ at 0.90): a well-captured cell would have shown the tape if it were there. $\Lambda$ is in
nats and is the log-likelihood a per-tape absorbing state would gain.

Clades: cells sharing the depth-$d$ prefix of an **anchor** tape *and* a clone, scored on every tape
$z\neq a$ (cross-tape), deduplicated by cell set, and collapsed so overlapping passing clades for one
(clone, tape) give **one** row. ⚠ Only clades some anchor resolves are visible, so the catalogue is a
**lower bound**.

### ⚠⚠ CORRECTION — the first version had no third margin and was useless

Scored against $\mathrm{logit}^{-1}(\alpha_c+\beta_z)$ alone, a tape lost **clone-wide** makes every
subset of that clone look spectacular — in real *and* permuted data alike. Those candidates flooded
both counts and **the FDR sat at 64–73% at every threshold**. Fixed by adding a per-(clone, tape)
offset fitted to the clone's own margin,

$$\sum_{c\in C}\mathrm{logit}^{-1}(\alpha_c+\beta_z+\gamma_{C,z}) = k_{C,z}$$

so $\Lambda$ measures only **within-clone** deviation. $\gamma$ fixes *how many* cells of the clone
lack the tape; it says nothing about *which* — and that is exactly what $\Lambda$ tests. A clone-wide
loss now gives $\tilde p\to1$ and $\Lambda\approx0$. **FDR fell from ~65% to ≈0%.**

### Results — all five arms, FDR ≈ 0

| arm | events | tapes | clones | median clade | inside rate vs expected | $\Lambda$ total | $R_c$ all / event |
|---|---|---|---|---|---|---|---|
| Subclone | 8,220 | 166 | 11 | 9 | 1.00 vs 0.11 | 258,295 | 131 / 131 |
| Pre-TX | 1,783 | 161 | 266 | 9 | 1.00 vs 0.25 | 23,794 | 128 / 131 |
| Mouse 2 | 388 | 126 | 12 | 14 | 1.00 vs 0.22 | 15,637 | 116 / 117 |
| Mouse 1 | 320 | 128 | 35 | 16 | 1.00 vs 0.34 | 7,097 | 114 / 116 |
| Mouse 3 | 73 | 51 | 9 | 23 | 1.00 vs 0.41 | 1,567 | 113 / 118 |

**Inside rate is 1.00 in every arm**, and **event cells are never worse captured** — the
capture-independence result, clean.

### ⚑⚑ Detection power is the whole story (`src/41_event_strat.py`)

Share of a clone-size band's missing entries inside a called event:

| clone size | Mouse 1 | Mouse 2 | Mouse 3 | Pre-TX | Subclone |
|---|---|---|---|---|---|
| 5–10 | 0.00% | 0.00% | 0.00% | 0.00% | — |
| 10–20 | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| 20–50 | 0.30% | 0.19% | 0.24% | 1.22% | 0.00% |
| 50–200 | 2.73% | 3.47% | 4.12% | 7.54% | — |
| 200+ | **8.56%** | **9.23%** | **6.52%** | — | **19.12%** |

**Not one event in any clone under 20 cells, in any arm** — $\gamma_{C,z}$ is fitted from that
clone's own cells and absorbs everything when there are few. So pooled percentages are dilution
artefacts. **Where detection is possible, 6.5–19% of missing entries sit inside a confident
inherited loss.**

⚠ **"Fraction of missingness explained" was the wrong phrase and is retired.** Because $\gamma$
matches each clone×tape total exactly, $\sum_{c\in C}(\text{miss}-\tilde p)=0$ within every
clone-tape cell — the excess inside an event is balanced by a deficit elsewhere in the same clone.
The quantity measures **concentration**, not an additive share. Report instead: missing entries
inside called blocks (3.64% of all missing in Mouse 1; 14,770 of 14,897 slots, i.e. 99.1%, confirming
completeness) and $\Lambda$ in nats.

### The clone-wide layer (`src/42_clonewide.py`)

Same statistic one level up: clade = the whole clone, scored against the **two**-margin fit
($\gamma$ is what we now want to measure rather than absorb); null permutes cells **across clones
within sample**.

| arm | losses | tapes | clones | inside vs expected | **% of all missing** | $\Lambda$ | FDR |
|---|---|---|---|---|---|---|---|
| Subclone | 155 | 71 | 12 | 0.94 vs 0.34 | **32.90%** | 263,004 | 0.1% |
| Mouse 2 | 410 | 130 | 69 | 1.00 vs 0.28 | **21.68%** | 25,427 | 4.4% |
| Mouse 1 | 725 | 154 | 118 | 0.97 vs 0.30 | **14.34%** | 40,712 | 0.0% |
| Mouse 3 | 398 | 133 | 61 | 1.00 vs 0.32 | **13.20%** | 13,715 | 0.9% |
| Pre-TX | 4,763 | 161 | 1,188 | 0.94 vs 0.17 | **7.10%** | 114,923 | 0.0% |

⚑ **The two layers are one phenomenon at two epochs.** Founder-39 was monoclonal and sorted on high
mRFP, so all 166 integrations were active at cloning; ClonalBC transduction was day −11 and the
bottleneck to ~8,000 clones day 3. So **clone-wide losses are silencing events before the bottleneck;
sub-clone events are after it.** Same mechanism, different time.

⚠ **The layers do NOT add to shares of a total.** Each fit matches its own margins, so every layer
sums to zero against the layer beneath. Report each against the one below — $\Lambda$ in nats and
entries inside called blocks — never as a partition.

### Soft variant: complete losses, or graded shifts too? (`src/43_soft_events.py`)

The hard test pins $H_1$ at $P(\text{missing})=1-\varepsilon$ with **no free parameter**, so a real
0.30→0.80 shift in a 20-cell clade earns only +2.11 nats against +23.88 for 0.30→1.00. The soft
variant fits the clade's own rate, $\hat\pi = k/m$ (the MLE: maximising
$k\log\pi+(m-k)\log(1-\pi)$ gives $k/\pi=(m-k)/(1-\pi)$, i.e. the observed fraction missing). Same
examples: **+10.68** and **+24.08**. $\Lambda_{\rm soft}\ge\Lambda_{\rm hard}$ always, since
$H_1^{\rm hard}$ is the point $\pi=1-\varepsilon$ of $H_1^{\rm soft}$.

⚠ **Correction: Wilks does not apply.** $H_0$ (per-cell $\tilde p$, no free parameters) is nested in
$H_1^{\rm soft}$ only if every $\tilde p$ in the clade is equal, which it is not. The models are
non-nested, $\Lambda_{\rm soft}$ can be negative, and there is no $\chi^2_1$ reference — the
permutation is the only calibration. One-sided ($\hat\pi>$ expected); the opposite direction is a
calibration check.

| arm | soft threshold | FDR | soft events | $\hat\pi$ median (expected) | $\hat\pi\ge0.99$ | $\hat\pi<0.90$ **partial** | hard keeps |
|---|---|---|---|---|---|---|---|
| Mouse 1 | 10.0 | 3.4% | 37,160 | 0.989 (0.409) | 48.5% | 20.3% | 76.3% |
| Mouse 2 | 10.0 | 3.1% | 22,646 | 1.000 (0.258) | 51.8% | 31.7% | 69.8% |
| Mouse 3 | 16.3 | 5.0% | 6,476 | 0.993 (0.512) | 53.6% | 7.9% | 48.3% |
| Pre-TX | 10.0 | 0.0% | 67,628 | 1.000 (0.254) | 67.3% | 13.8% | 81.4% |
| Subclone | 10.0 | 0.7% | 516,549 | **0.857** (0.115) | 26.4% | **57.2%** | 54.2% |

⇒ **The hard test is strict but not badly so.** In four arms the typical event is complete
($\hat\pi$ median 0.989–1.000 against 0.25–0.51 expected), about half exceed 0.99, and
$\hat\pi<0.90$ accounts for 8–32%. So a per-tape absorbing state covers most of it, with a real
minority that would want a lineage-varying rate instead.

⚠ **Subclone inverts this (57.2% "partial") and it is probably an artefact of clade resolution, not
biology.** Two different things are pooled under $\hat\pi<0.90$: genuine graded silencing, and **a
complete loss on a smaller clade than the one tested**. Prefix-defined clades stop at depth 4, so in
an arm whose clones run 200–11,000 cells they are coarse relative to the true tree, and a complete
loss on a sub-subclade reads as a partial loss on the clade we happened to score. Subclone has by far
the largest clones and by far the lowest expected rate (0.115), so both effects push the same way.
⇒ **8–32% is an upper bound on genuine partial silencing in the mice; 57% overstates it in
Subclone.**

⭐ **Cheap decisive test, not yet run.** Rebuild the prefix cache with `MAX_D = 6` instead of 4
(`35_lineage_depth.py` sets it; the tapes have six sites so the depth exists) and re-run. If the
"partial" events are complete losses on finer clades they will **resolve into complete events at
depth 5–6**; if they stay partial at maximum resolution the graded component is real. That
distinguishes *absorbing state* from *lineage-varying rate* event by event — the modelling question.

⚠ **Correction: "fitting $\varepsilon$" was the wrong idea.** Any estimate from called events is
selected on having few present cells, so it is biased downward. What is honest is a **sensitivity
sweep**, and it shows $\varepsilon$ is not load-bearing: across a 10× range (0.005→0.05) the share of
soft events the hard test keeps moves only 73.5→77.7% (Mouse 1), 68.5→74.0% (Mouse 2), 46.9→53.6%
(Mouse 3).

### ⭐ The `MAX_D = 6` test — run (2026-09-03)

`src/44_prefix_codes6.py` rebuilds prefix codes to the full six sites
(`prefix_codes6_{arm}.npz`, a new filename so depth-4 results stay reproducible; cell filter and row
order kept bit-identical and asserted against the depth-4 cache). `src/43_soft_events.py --maxd 6`
then reruns the soft catalogue at depths 1–6 and stratifies $\hat\pi$ **by the depth of the clade**.

The logic: if "partial" events ($\hat\pi<0.90$) are really complete losses on clades we scored too
coarsely, then **$\hat\pi$ must rise and the partial fraction must fall as clades get finer.** If
they persist at maximum resolution, the graded component is real.

**Result — both are true, in different proportions:**

| arm | median $\hat\pi$ d1 → d6 | $\hat\pi\ge0.99$ d1 → d6 | **$\hat\pi<0.90$ d1 → d6** |
|---|---|---|---|
| Mouse 3 | 0.969 → **1.000** | 43.8% → **70.5%** | 12.7% → **5.3%** |
| Mouse 2 | 0.969 → **1.000** | 46.1% → 56.6% | 37.5% → **26.0%** |
| Mouse 1 | 0.984 → **1.000** | 44.9% → 52.0% | 22.9% → 19.8% |

⇒ **Clade resolution explains a substantial part of the "partial" fraction but not all of it.**
Finer clades give more complete losses in every arm — Mouse 3 falls to 5.3% partial, essentially all
artefact. But Mouse 1 (19.8%) and Mouse 2 (26.0%) still carry a partial component at **maximum
recorder resolution**, so a genuine graded element survives.

⚠ Depth is not free: determined (cell, tape) pairs fall from 76–78% at depth 1 to 21–24% at depth 6
(Pre-TX to 3.9%), so deep clades are scarcer and smaller. The monotone trend is read within that.

**⇒ Modelling reading.** A per-tape absorbing state covers the majority — complete losses are ~50–70%
of events at full resolution and rise with resolution — but a ~20–26% graded residue in the mice is
not an artefact of coarse clades and would want a lineage-varying rate. So the answer to "absorbing
state or latent field" is **mostly the former, with a real minority of the latter**, and that
minority is now bounded rather than assumed.

**All five arms, $\hat\pi$ stratified by clade depth:**

| arm | max clone | $\hat\pi\ge0.99$ d1 → d6 | **$\hat\pi<0.90$ d1 → d6** | median $\hat\pi$ at d6 | median clade at d6 |
|---|---|---|---|---|---|
| Mouse 3 | 210 | 43.8% → 70.5% | 12.7% → **5.3%** | 1.000 | 58 |
| Pre-TX | 127 | 59.7% → 85.7% | 15.1% → **6.8%** | 1.000 | 7 |
| Mouse 1 | 1,607 | 44.9% → 52.0% | 22.9% → **19.8%** | 1.000 | 20 |
| Mouse 2 | 3,387 | 46.1% → 56.6% | 37.5% → **26.0%** | 1.000 | 18 |
| Subclone | 10,996 | **5.0% → 44.7%** | **85.2% → 37.7%** | 0.947 | 12 |

⚠⚠ **CORRECTION to the reading written before the last two arms landed.** I concluded there that a
"~20–26% graded residue is not an artefact of coarse clades". The completed set says otherwise:
**the residual partial fraction at depth 6 is monotone in maximum clone size** — Mouse 3 (210 cells)
5.3%, Pre-TX (127) 6.8%, Mouse 1 (1,607) 19.8%, Mouse 2 (3,387) 26.0%, Subclone (10,996) 37.7%.

That is the signature of the *recorder running out of resolution before the tree does*. Six sites
resolve at most six levels; in arms whose clones are small, six levels suffice and the partial
fraction converges to **5–7%**. In arms with thousand-cell clones the trees are far deeper than six
levels, so clades stay coarse and partial events persist. Subclone makes it plainest: $\hat\pi$ runs
0.507 → 0.947 and the partial fraction falls 85.2% → 37.7%, **still steeply falling at the resolution
limit**, with clades still a median of 12 cells inside clones of up to 10,996.

⇒ **Revised reading: the graded appearance is very largely clade coarseness, not a graded mechanism.**
Extrapolating to the arms where six levels *are* enough, the genuine partial component looks like
**~5–7%, not 20–26%.** For the model that strengthens the earlier conclusion rather than qualifying
it: **a per-tape irreversible absorbing state is the right and largely sufficient extension**, and
the case for a lineage-varying rate is weaker than the depth-4 numbers implied.

⚠ This is an inference from a trend across arms, not a direct measurement — the recorder cannot be
pushed past six levels. A simulator with known ground truth could settle it, which is another use for
the one Figs 5–7 need.
