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

### Is $\Lambda_{\rm soft}$ just "we fitted a parameter"? — measured (`src/45_lrt_calibration.py`)

$\hat\pi = k/m$ matches the clade's observed count exactly, so the natural worry is that $H_1$ is
guaranteed to beat $H_0$ and the statistic is uninformative. **Measured on all 5,541,360 candidate
(clade, tape) pairs in Mouse 3, with no threshold:**

| | mean $\Lambda_{\rm soft}$ | sd |
|---|---|---|
| observed | **−0.824** | 2.79 |
| within-clone permutation | **−1.740** | 3.90 |
| Wilks reference for a 1-parameter fit | +0.500 | — |

⚠⚠ **CORRECTION.** I earlier wrote that the free fitting advantage is "$\approx0.5$ nats, by Wilks".
**Wrong on both counts.** Wilks does not apply — $H_0$ is not nested in $H_1$, since $H_0$ gives every
cell its own $\tilde p_c$ (varying through $\alpha_c$) while $H_1$ imposes one shared $\pi$. And the
measured null mean is **−1.74**, not +0.5: on a typical clade the fitted-constant model is *worse*,
because one degree of freedom does not buy back the per-cell structure it discarded. **A positive
$\Lambda_{\rm soft}$ therefore has to overcome a deficit first — it is not a fitting artefact.**

Enrichment over the permutation null, by threshold:

| $\Lambda_{\rm soft}\ge$ | 2 | 4 | 10 | 16 | 25 |
|---|---|---|---|---|---|
| observed | 150,463 | 52,169 | 16,301 | 8,954 | **3,705** |
| null | 83,185 | 11,246 | 2,929 | 1,015 | **0** |
| enrichment | 1.8× | 4.6× | 5.6× | 8.8× | — |

**What 10 nats costs**, for a 20-cell clade at a clone rate of 0.30 (≈6 missing expected, sd 2.05):
$k=8$ (+1 sd) → $\Lambda=0.45$; $k=14$ → 6.79; $k=15$ → 8.59; $k=16$ → **10.68**. A one-sd
fluctuation buys under half a nat; the threshold needs a **+4.6 sd** departure.

**On "why not just a likelihood ratio test":** $\Lambda_{\rm soft}$ *is* one —
$\log[L(H_1)/L(H_0)]$. What does not apply is the *standard* LRT's $\chi^2$ p-value, for the
non-nesting reason above; hence permutation. A genuinely simpler alternative is a **2×2 test**
(clade vs rest-of-clone × missing vs present), which would rank events similarly — but it weights
every cell equally and so cannot distinguish "the missing cells are the well-captured ones" from
"the missing cells are the poorly-captured ones". The per-cell $\tilde p_c$ inside $\Lambda$ is what
makes this a test of *inherited loss* rather than of *dropout rate*, and that is the whole
discriminator against technical dropout.

### B = 1,000 permutations — launched (2026-09-03)

Every catalogue so far ran on $B=3$–200 permutations, so its permutation $p$-value was floored at
$1/(B+1)$: the *formal* claim for the event catalogue was only $p<0.17$, however overwhelming the
counts (28,367 candidates against ~1.4 null). This raises $B$ to 1,000 for **all three layers × all
five arms**, and for the soft catalogue at **both** clade resolutions (depth 4 and the full
recorder depth 6) — 15 configurations, plus the five clone-wide runs.

**Cost, and why it needs splitting.** The entire scan — every anchor × depth × clade × tape — is
redone per permutation. Measured per-scan cost (elapsed ÷ (nperm+1), from a completed `sacct`):
`42` 0.02–0.12 s, `40` 2–40 s, `43` 9 s–3.2 min. At $B=1{,}000$ that is ~2 min for `42` but ~11 h
for `40` on Initial and ~52 h for `43` on Subclone, serially.

**The design (`src/47_submit_B1000.sh`).**

| piece | what it does |
|---|---|
| `--permpart i/N` in `40`/`43` | runs permutations $[\,(i{-}1)B/N,\;iB/N)$ and writes **only a count vector**, then exits without doing the observed scan at all |
| a **fixed** threshold grid | 600 geometric points, 10 → $10^7$ nats, ratio 1.023, *independent of the data* — count vectors from different array tasks are addable only if every task scored the same thresholds |
| seed $(\text{SEED}, b)$ | permutation $b$ is seeded from its own index, not from a stream advanced $b$ times, so slice boundaries do not change what any permutation is, parts merge in any order, and a duplicated part is detectable |
| `46_perm_merge.py` | pools the parts, **asserts the set is complete and disjoint**, writes `permnull_{script}_{arm}{tag}.npz` |
| `48_collect_B1000.sh` | re-runs each observed scan **once** with `--nullfile`, attaching $q$ and the global $p$ |

⚠ **Storing counts, not candidate lists, is what makes this affordable.** Subclone's soft scan
produces ~876,000 candidates per permutation; 1,000 of those lists would be ~7 GB and useless
afterwards. The count vector is 600 integers and is everything the FDR curve, the $q$-values and the
global $p$ need.

**The two quantities, and why the whole per-permutation matrix is kept rather than a mean.**

- **Global $p$** — $p=\bigl(1+\#\{b: C_b\ge C_{\rm obs}\}\bigr)/(B+1)$ with $C$ the total candidate
  count above a threshold. Needs the *distribution* over $b$, so the $(B\times G)$ matrix is stored.
  With $B=1{,}000$ and zero exceedances this licenses $p<0.001$.
- **Per-event $q$** — $q(\Lambda)=\overline{\text{null}}(\ge\Lambda)\,/\,\text{obs}(\ge\Lambda)$,
  monotonised, then attached as a column of `events_{arm}.tsv.gz` and `clonewide_{arm}.tsv.gz`.
  Needs the *mean*.

⚠ **Correction to the task spec: the monotonising running minimum goes UP the grid, not down.**
The spec said "made monotone by a running minimum from the top". Implemented the other way, and the
reason is BH's: an event at $\Lambda$ can be reported by *any* rejection region $\{\Lambda'\ge t\}$
containing it, i.e. any $t\le\Lambda$, so its $q$ is $\min_{t\le\Lambda} q_{\rm raw}(t)$ — a running
minimum from the low-$\Lambda$ end. That is automatically non-increasing in $\Lambda$, which is the
direction a $q$-value must run. A running minimum from the high end would do the opposite and would
also drag the $q_{\rm raw}=1$ convention — which the code assigns wherever the observed count is
zero, i.e. at every threshold above the largest observed $\Lambda$ — back down over the entire
informative range.

**⚠ Partition: `cpu` only, not the project default `-p lesliec,cpu`.** Listed against both, 300
one-core tasks put **194 CPUs on `lesliec` — 76% of the lab's four private nodes, held by one
user**. Those four nodes are also the lab's *only* GPU nodes, so pure-CPU work parked there can
block a labmate's A100 job on CPUs while the GPUs sit idle. The general partition has 239 nodes and
~9,700 idle CPUs, 30× what this needs. Resubmitted `-p cpu`: `lesliec` went back to 233 of 256 CPUs
idle, and this work now occupies 309 of 14,264 CPUs on `cpu` (2.2%). **Rule for the future: a job
that needs neither a GPU nor a ~1 TB node has no business on `lesliec`**, and a wide array is
exactly where the default job-pair listing does the most damage.

**Sizing.** 1 core, 8 G, 20 array tasks per configuration. Peak RSS across every `40`/`43` run to
date is 935 MB, so 8 G is ~8× headroom; the inner loop is `bincount`/`exp`, not BLAS, so extra cores
buy nothing and a 1-core 8 G task backfills into gaps a fat one cannot. Walltimes 6 h (16 h for
`43` on Initial/Subclone) against a predicted worst case of 2.7 h. **All 300 array tasks entered
`RUNNING` immediately.**

**Verified before launch** (`logs/smoke_B1000-11388862.out`, Mouse 3, small $B$): parts → merge →
observed-with-`--nullfile` reproduces the previously published Mouse 3 numbers exactly (73 events,
$\Lambda$ total 1,567 nats, 398 clone-wide losses, 13.20% of missing), and the merge **refuses** an
incomplete set, naming the missing permutation indices.


---

## Two directions from the B = 1,000 run (recorded 2026-09-03, to pick up later)

### Direction 1 — how widespread is heritable silencing?

**The idea, as put:** with a properly calibrated null there should be more (clade, tape) combos
that clear significance but sat below the earlier effect-size cutoff, surfacing more *partial*
silencing; and we can then count how many **cell × tape** entries sit inside a called event, and
still separate complete from partial via $\hat\pi$.

**This is the right question, and the quantity to lead with is the cell × tape count, not the event
count.** An "event" is an artefact of how the clade search happens to carve the tree — the
catalogue is explicitly a lower bound on events. Cell × tape entries inside called blocks is the
quantity that says *what a per-tape absorbing state would actually buy the likelihood*, it is
already partly measured (6.5–19% of missing entries where clones are large enough to detect
anything), and it is comparable across layers and arms. ⚠ It is a **count inside called blocks**,
not a share of a total — the retired "fraction of missingness explained" phrasing is wrong for the
reason recorded above ($\gamma$ forces $\sum_{c\in C}(\text{miss}-\tilde p)=0$ within every
clone × tape cell).

**⚠⚠ But $B$ alone will not surface the weaker events, and this run cannot answer below 10 nats.**
Three things are being conflated, and separating them is the whole content of this note:

| knob | what it controls | what $B=1{,}000$ does to it |
|---|---|---|
| $B$ | resolution of the $p$-value; precision of the null count at each threshold | fixes it — $p$ floor $0.17\to0.001$ |
| the FDR **threshold** $\Lambda$ | which candidates are *called* | lowers it where a noisy small-$B$ null had pushed it up |
| the scan **floor** `--lam` | which candidates are *collected at all* | **nothing — it is fixed at 10 nats** |

The scan floor is a hard cutoff inside `scan()`: nothing below `LAM0` is ever put in the candidate
list, and this run's grid starts at exactly 10. So the run will move the *threshold* — Mouse 3's
soft threshold of 16.3 nats came from **three** permutations and should fall toward 10 once the null
is properly estimated — but it can never look beneath 10.

**⇒ The follow-up this implies: re-run `40`/`43` with `--lam 4`.** $B=1{,}000$ is precisely what
makes a low floor *interpretable*: at $B=3$ the null count at 4 nats is far too noisy to divide by.
Script `45` already measured the terrain on Mouse 3 — observed vs null 52,169 / 11,246 at
$\Lambda\ge4$ (4.6× enrichment) against 16,301 / 2,929 at $\ge10$ (5.6×) — so there is real signal
below the current floor, at a worse but possibly tolerable FDR. A floor-4 run **subsumes** the
floor-10 one (its grid contains 10), so it is a replacement, not an addition, and costs ~2–3× the
current run.

**Decision taken: do not restart the running jobs to do this.** The current run delivers the
headline $p$-values tonight, its FDR curve at the floor is exactly what tells us how much room lies
below 10, and 250 core-hours is not worth protecting on a 14,264-CPU partition. Read the curve
first, then launch floor-4 with a sizing measured from this run.

⚠ **Expect a chunk of the new low-$\Lambda$ events to be clade coarseness, not partial silencing.**
That is the settled reading of the `MAX_D=6` test: the residual partial fraction is monotone in
maximum clone size (5.3% at 210 cells → 37.7% at 10,996), i.e. the recorder runs out of resolution
before the tree does. Any claim that a lower floor reveals *graded* silencing has to survive the
$\hat\pi$-by-clade-depth stratification, which is already built into `43`.

### Direction 2 — ⭐ does the co-integrated symbol disappear when a tape is silenced?

**This is the sharper of the two, and it is the first genuinely orthogonal test of the mechanism.**
Everything so far infers silencing from *missingness*, which is also what a technical dropout looks
like; capture-independence argues against the technical reading but does not exclude it. This
predicts a signal in a **completely different channel — the symbols written into other tapes** —
where transcript capture cannot reach.

**The mechanism, and why the prediction is not tautological.** The cassette is
`PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC`: one integration carries both a pegRNA (bearing a
4-nt insert barcode `NNNN`) and a tape (labelled by a 10-nt `TargetBC`). pegRNAs act in ***trans***
(§0, four independent grounds), so integration $z$'s symbol is written into **every** tape in the
cell, not its own. Therefore, if silencing an integration kills the whole locus:

> losing tape $z$ from the readout ⇒ symbol $s(z)$ should also stop appearing **at all other tapes**
> in the same cells.

⚠ **The two promoters are different polymerases** — pegRNA on U6 (Pol III), tape/mRFP on EF1α
(Pol II). Locus-level heterochromatin should take both, but nothing *forces* it, so the coupling is
a hypothesis, not a deduction. That is exactly why the test is worth running: it measures the thing
§0 currently asserts ("nasty coupling", row A9) and Fig 3b interprets ($\beta_z$ = the integration's
expression level).

**⚠ The blocker: the $z \mapsto s(z)$ map is not known.** `TargetBC` (10-nt, tape) and `NNNN`
(4-nt, symbol) are different barcode spaces, and nothing in the delivered tables links them —
the same gap §"Open empirical question" records for the Typewriter lineage data. Recovering it
would need sequencing of the intact integrations, which we do not have.

**⇒ Invert it: do not assume the map, *recover* it.** Rather than testing a known pair, ask of each
silencing event *which* symbol is depleted. The map becomes the output, and its internal structure
becomes the test — because a spurious map has no reason to be consistent, injective, or reproducible.

**The validation ladder — this is what makes it convincing, not the per-event p-values:**

1. **Cross-arm reproducibility (the killer).** Verified 2026-09-03: **all 166 TargetBCs are
   identical across all five arms** — same engineered line, same integrations. So tape $z$ must name
   the *same* symbol in Mouse 1, Mouse 2, Mouse 3, Pre-TX and Subclone, recovered independently.
   Agreement by chance among ~100 symbols is 1%. Five-way concordance is essentially unfakeable.
2. **Injectivity and the collision structure.** 166 integrations drawing `NNNN` uniformly from
   $4^4=256$ predicts $256(1-(1-1/256)^{166}) = \mathbf{122}$ distinct symbols. We observe
   **100–106** design-conforming symbols carrying essentially all edits (103 Pre-TX, 100 Mouse 1,
   106 Subclone). The recovered map should be near-injective with *exactly* that collision rate —
   and collisions are not a nuisance but a prediction: a symbol carried by two integrations should
   show only a **partial** drop when one is silenced.
3. **Dose.** A symbol named by $j$ tapes should carry ~$j\times$ the base $\xi$.
4. **$\mathrm{corr}(\beta_z,\ \xi_{s(z)})>0$** across integrations — tape recovery rate and symbol
   frequency are two readouts of one locus's expression. ⚠ Measured today: $\xi$ is **smooth over a
   570× range** (0.00008–0.0456), *not* quantised into copy-number multiples, so per-integration
   expression varies enormously. That kills quantisation as a cheap shortcut but strengthens the
   premise of test 4 — there is a large dynamic range for the two readouts to correlate over.

**⚠ Step 0 first, because it decides whether any of this has power: how much of the edit content
postdates the loss?** The symbols already written into a tape before an integration was silenced
stay there — the tape is an append-only record. So only insertions laid down *after* the silencing
can show the depletion, and everything inherited from before dilutes it. Tapes here are close to
saturation (mean ~4.5–5 of 6 sites edited), which is precisely the regime where most content is
ancestral. **The measurement:** the fraction of (tape, site) slots that are **polymorphic within a
clone** (⇒ written after the clone founder) and, one level down, polymorphic within a called clade
(⇒ written after its MRCA). That single number is the power calculation for the whole idea, it is
cheap, and it should be computed before any examples are pulled.

**Recommended order** — and the user's instinct to look at examples first is right, with Step 0
slotted ahead of it so we know what we are looking at:

| step | what | why here |
|---|---|---|
| 0 | polymorphic-slot fraction, within clone and within clade, per arm | decides feasibility; cheap; also tells us which layer to use |
| 1 | **hand-inspect 3–5 examples**: highest-$\Lambda$, largest-clade events in Subclone/Mouse 2 big clones. Just print the symbol-frequency table inside vs outside | what was asked for; catches design errors no aggregate would |
| 2 | screen on the **clone-wide** layer first, not the sub-clone one | it is the strongest layer — Pre-TX alone has 4,763 losses over 1,188 clones, ~30 losing clones per tape — and the loss predates the clone founder, so a larger share of the clone's editing postdates it |
| 3 | the validation ladder above | this is the actual evidence |
| 4 | only then the sub-clone version, restricted to post-MRCA insertions | the genuinely *lineage* claim, but power-limited by Step 0 |

**The statistic for steps 2/4.** For tape $z$ and symbol $s$, contrast $s$'s share of insertions in
cells/clones that lost $z$ against those that did not, on **post-loss insertions only** (within a
clade: positions beyond the clade's common prefix on each tape; within a clone: clone-polymorphic
slots). Two controls fall out of the design: the rest of the symbol vector is an internal control
against a general compositional shift (removing $\xi_s\approx1\%$ renormalises everything else up by
~1%, which must be divided out), and the tape being scored is **never** the tape supplying the
insertions — the same cross-tape rule that keeps $\Lambda$ honest.

⚠ **A global alternative worth remembering:** build $D[z,s]$ = depletion of $s$ in cells lacking $z$
and solve it as an **assignment problem** (Hungarian) rather than tape-by-tape. One statistic, all
166 assignments at once, compared against permuted $D$. More powerful than 166 separate tests, and
the natural form of the question — but do it after the per-tape version, so the examples stay
inspectable.

### ⚠⚠ Per-combo permutation p-values: well-defined, cheap, and the wrong tool at this floor

**The proposal.** Rather than a pooled count, give every (clade, tape) combo its own permutation
p-value: $p(S,z) = \bigl(1+\#\{b:\Lambda_b(S,z)\ge\Lambda_{\rm obs}(S,z)\}\bigr)/(B+1)$, then take a
threshold much stricter than 0.05.

**It is well-defined.** Prefix codes are never permuted, so a clade *slot* persists across
permutations with its size fixed; $\Lambda_b(S,z)$ is the score that slot gets when its cells are
random members of the clone. That is exactly the right conditional null.

**It is cheap.** Not by storing $B$ values per combo — Mouse 3 alone has 5,541,360 combos, so 1,000
each is ~44 GB — but with **one exceedance counter per combo**, incremented in place: ~22 MB.

**⚠ But at the 10-nat floor it is strictly weaker than the pooled count, and the arithmetic is not
close.** Both spend the same 1,000 permutations; they differ in how many null draws that buys:

| | null draws | resolvable tail probability |
|---|---|---|
| per-combo | $B=1{,}000$ | $10^{-3}$ |
| pooled count | $B\times N_{\rm combos}\approx10^{9}$ | $\sim10^{-7}$ |

Measured: **eight** permutations of Mouse 3's hard scan produced **zero** null candidates above 10
nats, across every combo. So the per-combo exceedance probability is below $\sim10^{-6}$ and 1,000
draws return zero essentially always — every real event comes back at exactly $1/1001$, censored,
unrankable. "The best examples exceed every random permutation" is therefore the *failure* mode of
the statistic, not its payoff: it saturates precisely where the signal is strongest. The pooled
count estimates the same tail to ~3% relative precision by borrowing strength across combos — at
the cost of assuming $\Lambda$'s null is exchangeable *between* combos, which is what stratification
below repairs.

**⚠⚠ A second floor on $p$, and this one no $B$ can lift.** A within-clone permutation of an
$m$-cell clade inside an $n_C$-cell clone can only realise $\binom{n_C}{m}$ distinct compositions,
so $p \ge 1/\binom{n_C}{m}$ **however large $B$ is**:

| clone $n_C$ | clade $m$ | compositions | best achievable $p$ |
|---|---|---|---|
| 6 | 4 | 15 | 0.067 — **cannot reach 0.05** |
| 10 | 4 | 210 | 0.0048 |
| 20 | 9 | 167,960 | $\sim10^{-5}$ |

⇒ This is a mechanism for "**not one event in any clone under 20 cells**" (script `41`) beyond
$\gamma_{C,z}$ absorbing everything: the permutation null of a small clone has too few states to
certify anything. Worth saying in the talk — it makes the power limit structural rather than a
choice of threshold.

**⇒ What to actually do, all three in the `--lam 4` follow-up, not in the run now finishing:**

1. **Lower the scan floor to 4 nats** — still the binding constraint.
2. **Stratify the pooled null by clade size** (and expected rate). This is the real repair for
   "weaker but real". Pooling assumes one null distribution for all combos, and a 500-cell clade has
   a far heavier $\Lambda$ tail than a 5-cell one, so a single global FDR judges both by an average
   fitting neither. Per-stratum count vectors cost ~6× the storage — nothing — and keep almost all
   the resolution while conditioning on what makes a combo easy or hard.
3. **Then add the per-combo exceedance counter, where it earns its keep.** At $\Lambda\ge4$ the null
   is substantial (Mouse 3 soft: 11,246 null vs 52,169 observed), so per-combo p-values there are
   *not* censored. That is exactly the "beaten by a few permutations but still significant" regime.

**⚠ Terminology, to keep straight in the talk.** The current run already attaches a per-event
number — a **$q$-value**, the FDR of the rejection region containing that event. What it does not
give is a per-event **$p$**. They answer different questions ("what fraction of calls at this
threshold are null?" vs "how exceptional is this one combo?"), and for a catalogue the $q$ is the
more useful of the two.

### B = 1,000 — results, all 15 configurations (2026-09-04)

`perm_collect` (11389784) pooled 300 parts into 15 nulls of exactly 1,000 permutations each and
re-ran every observed scan against them. 17 min. Every merge passed the completeness/disjointness
assertion.

**$p = 0.000999 = 1/1001$ in all twenty runs** — hard, soft (both depths), and clone-wide, in every
arm, at both the chosen threshold and the 10-nat scan floor. Not one of 1,000 permuted labellings
produced as many candidates as the real one, anywhere.

**⭐ The number to say out loud is the null *maximum*, not the p-value.** "Where does the real
labelling sit among the 1,000 random ones?" has a much better answer than "$p<0.001$":

| arm | observed | null mean | **null MAX over 1,000** | obs / null max |
|---|---|---|---|---|
| Subclone | 279,973 | 43.2 | **111** | **2,522×** |
| Pre-TX | 55,076 | 10.5 | 104 | 530× |
| Mouse 3 | 8,538 | 0.17 | 17 | 502× |
| Mouse 1 | 28,367 | 2.15 | 60 | 473× |
| Mouse 2 | 15,801 | 2.17 | 134 | 118× |

⇒ **the most extreme of a thousand random labellings reached 111 candidates where the real one
reached 279,973.** That makes the null explicit, which a p-value hides, and it is immune to the
objection that $1/1001$ is just the floor of what $B$ can resolve.

⚠ Note the null is **strongly right-skewed** — max/mean runs 6× (Pre-TX) to 100× (Mouse 3). That is
the same skew flagged for the per-tape nulls in script `37`, and it is why $B=3$ could never have
characterised this tail. It is also why the *max* is the honest summary: quoting the mean alone
understates what a lucky permutation can do.

**Every event in the hard catalogue clears $q\le1.9\times10^{-4}$** (max $q$ over kept events:
Mouse 3 $2.0\times10^{-5}$, Mouse 1 $7.6\times10^{-5}$, Mouse 2 $1.4\times10^{-4}$, Subclone
$1.5\times10^{-4}$, Pre-TX $1.9\times10^{-4}$). So the catalogue is significant *as a whole*, not
just at its top — there is no weak tail to defend.

**Cell × tape entries inside called blocks** — the quantity to lead with rather than event counts:

| arm | events | **cell × tape entries** | of all missing | cells touched |
|---|---|---|---|---|
| Subclone | 8,220 | **266,688** | 19.09% | 37,785 |
| Pre-TX | 1,783 | 22,256 | 1.70% | 8,233 |
| Mouse 2 | 388 | 21,277 | 6.43% | 4,013 |
| Mouse 1 | 320 | 14,770 | 3.64% | 3,842 |
| Mouse 3 | 73 | 2,916 | 2.28% | 701 |

⚠ Read these with script `41`'s stratification: the pooled share is a **dilution artefact**, since no
event is called in any clone under 20 cells. Where detection is possible it is 6.5–19%.

**Full table, all fifteen configurations:**

| config | arm | threshold | FDR | soft/hard events | null mean | null max |
|---|---|---|---|---|---|---|
| `40` hard | Mouse 1 / 2 / 3 / Pre-TX / Subclone | 10.0 all | 0.00–0.02% | 320 / 388 / 73 / 1,783 / 8,220 | 0.17–43.2 | 17–134 |
| `43` soft d4 | Mouse 1 / 2 / Pre-TX / Subclone | 10.0 | 0.02–3.3% | 37,160 / 22,646 / 67,628 / 516,549 | 15–3,931 | 153–4,860 |
| `43` soft d4 | Mouse 3 | 16.2 | 4.8% | 6,500 | 313 | 1,123 |
| `43` soft d6 | Mouse 1 / 2 / Pre-TX / Subclone | 10.0 | 0.02–3.0% | 48,794 / 37,083 / 69,185 / 876,479 | 16–4,422 | 154–5,462 |
| `43` soft d6 | Mouse 3 | 15.5 | 5.0% | 9,097 | 452 | 1,603 |
| `42` clone-wide | all five | 10.0 | 0.01–4.4% | 725 / 410 / 398 / 4,763 / 155 | 0.11–17.9 | — |

Soft-event counts rise at depth 6 in every arm (Mouse 2 22,646 → 37,083; Subclone 516,549 →
876,479), consistent with finer clades resolving losses the coarse ones smeared.

⚠ **The clone-wide layer on Mouse 2 is the one genuinely marginal result**: FDR 4.36%, null mean
17.9, max $q$ **0.044**. Every other configuration sits orders of magnitude clear. Say so rather
than letting it ride on the shared "$p<0.001$".

#### ⚠⚠ Correction: I predicted the thresholds would fall, and they did not

Written yesterday: *"Mouse 3's soft threshold of 16.3 nats came from **three** permutations and
should fall toward 10 once the null is properly estimated."* **Wrong.** At $B=1{,}000$ it is
**16.2** at depth 4 and 15.5 at depth 6. Every other configuration stayed at exactly the 10-nat
scan floor, where it already was.

⇒ **The three-permutation null was already unbiased in the mean; what it could not see was the
tail.** $B=1{,}000$ bought (i) resolution on $p$, from $<0.17$ to $<0.001$, and (ii) the *shape* of
the null — the 6–100× right skew above, which is new information and changes how the result should
be quoted. It did **not** move a single detection threshold. That is worth saying plainly: the
extra 250 core-hours bought rhetoric and calibration, not any new events.

⇒ **And it sharpens the case for the `--lam 4` follow-up.** If more events are wanted, the floor is
now demonstrably the only place they can come from — the threshold is already sitting on it in
14 of 15 configurations.

### ⚠⚠ How events are actually called — nail this down before any "how widespread" claim

Justin asked (2026-09-04) whether the 10-nat cutoff is arbitrary or adaptive. **It is both, and the
arbitrary one is doing the work.** Two distinct thresholds:

| | what it is | value |
|---|---|---|
| **scan floor** `LAM0` | hard-coded, fixed before seeing data. `scan()` only puts pairs with `L >= LAM0` in the candidate list — nothing below it exists downstream **at all** | **10 nats** |
| **calling threshold** `LAM` | adaptive: smallest grid point where the monotonised FDR curve reaches $\le5\%$ | chosen per config |

**The adaptive step is degenerate in 13 of 15 configurations** — `LAM` snaps to the floor because the
FDR is *already* far below target there:

| config | FDR at floor | × below the 5% target | chosen |
|---|---|---|---|
| `40` hard × 5 arms | 0.002–0.019% | **262× – 2,511×** | 10.0 |
| `43` soft, Pre-TX / Subclone / Mouse 1 / Mouse 2 (both depths) | 0.02–3.3% | 1.5× – 219× | 10.0 |
| `43` soft **Mouse 3** d4 / d6 | **10.03% / 9.22%** | 0.5× | **16.2 / 15.5** ← only binding case |

⇒ **For the hard catalogue the operative threshold is the arbitrary 10 nats**, and the stated 5%
FDR rule never bites. We are hundreds to thousands of times stricter than the criterion we claim.

#### Consequence 1 — event counts are set by the constant, and the ranking is not stable

| $\Lambda\ge$ | Subclone | Pre-TX | Mouse 2 | Mouse 1 | Mouse 3 |
|---|---|---|---|---|---|
| 10 | 8,220 | **1,783** | 388 | 320 | 73 |
| 15 | 4,667 | 390 | 230 | 164 | 48 |
| 20 | 3,102 | **110** | 174 | 114 | 37 |

**Pre-TX falls from 2nd to 4th between $\Lambda\ge10$ and $\ge20$.** Its events pile against the
floor — 78% within 5 nats of it, median $\Lambda$ 12.0, max 35 — because its clades are small
(median 9 cells). Subclone runs out to $\Lambda=1{,}326$. So "which arm has the most silencing" is
threshold-dependent and must not be asserted from one cut.

#### Consequence 2 — the "fraction of dropout" number is a curve, not a number

Cell × tape entries inside called blocks, as % of all missing:

| $\Lambda\ge$ | Subclone | Mouse 2 | Mouse 1 | Mouse 3 | Pre-TX |
|---|---|---|---|---|---|
| 10 | 19.09% | 6.43% | 3.64% | 2.28% | 1.70% |
| 20 | 13.62% | 5.38% | 2.23% | 1.70% | 0.20% |
| 100 | 5.49% | 2.24% | 0.61% | — | — |

⇒ **Never quote a single figure.** Quote the curve, or two thresholds, and say which.

#### ⚑ But the floor is *conservatism*, not rigour — which cuts the other way

At FDR 0.002–0.019%, the events at $\Lambda=10$–12 are almost certainly real; we are nowhere near a
false-positive cliff. The reported counts are therefore **lower bounds, and demonstrably loose
ones.** That converts the `--lam 4` run from "fishing" into "actually applying the criterion we
already state".

Four further reasons the counts are lower bounds, independent of the floor:
`MIN_CLADE = 4`; clades exist only where some **anchor tape resolves them** at depth $\le4$ (or 6);
$\gamma_{C,z}$ absorbs anything **clone-wide** by construction (that is the `42` layer); and nothing
is called in clones under 20 cells (script `41`).

#### What *is* clean

- **No double counting.** Dedup collapses overlapping clades for one (clone, tape) to the
  highest-$\Lambda$ one, so a cell × tape entry cannot be counted twice; non-overlapping ones are
  genuinely separate losses. ⚠ It keeps the *coarsest* passing clade, which is the conservative
  choice for event count and the liberal one for entries per event.
- **Kept events are all-or-none, not marginal enrichments.** `inside_rate` $\ge0.90$ at the 10th
  percentile in every arm, median **exactly 1.000**, against an expected rate of 0.11–0.41.

#### ⚠ And the 5% target itself is a choice that should not be inherited

With candidate counts in the $10^5$, **a 5% FDR admits thousands of false events** — 5% of
Subclone's 279,973 candidates is ~14,000. For a catalogue meant to be inspected event by event, an
**absolute** criterion ("threshold such that the expected number of false events $\le10$") is more
honest than a rate, and at the current floor we already satisfy something far stricter. Decide this
deliberately when setting the floor-4 threshold rather than defaulting to 0.05.

**⇒ Claims that survive any threshold choice**, and which should therefore carry the talk: the null
comparison (118–2,522× beyond the most extreme of 1,000 permutations), $q\le1.9\times10^{-4}$ for
every event in the hard catalogue, and $\hat\pi$ median 1.000 against an expected 0.11–0.41.

### The low-floor, size-stratified run — launched 2026-09-04 (`src/49_submit_lowfloor.sh`)

Answering "should we proceed with `--lam 4`": **yes, but the sizing exposed two things that changed
the design.**

#### ⚑⚑ Finding 1 — a fixed $\Lambda$ floor is a CLADE-SIZE FILTER in disguise

Extrapolating the candidate curve below the floor gave Pre-TX a log-log slope of **−4.4** against
−0.95 to −1.9 for every other arm. That is not noise, it is a cap. A fully-missing clade of $m$
cells at expected rate $p$ scores about $m\log\bigl((1-\varepsilon)/p\bigr)$ nats:

| arm | max clone | median clade | **cap for a median clade** | max $\Lambda$ seen |
|---|---|---|---|---|
| Pre-TX | 75 | 9 | **12.4 nats** | 35 |
| Mouse 1 | 1,131 | 16 | 17.1 | 150 |
| Subclone | 4,967 | 9 | 19.6 | **1,326** |
| Mouse 3 | 151 | 23 | 20.1 | 45 |
| Mouse 2 | 2,588 | 14 | 21.2 | 634 |

⇒ **a 9-cell Pre-TX clade cannot exceed ~12.4 nats however complete the loss.** A single 10-nat
threshold therefore excludes small clades *by construction, regardless of how real the loss is* —
which is exactly why Pre-TX collapses 1,783 → 110 events between 10 and 20 nats and falls from
second place to fourth. **So clade-size stratification is not a refinement; without it a lower floor
would simply pile up small-clade candidates and be judged against a null dominated by large ones.**

Strata: `[4,6) [6,10) [10,20) [20,50) [50,200) [200,∞)`. Counts are kept **per stratum**, the
threshold is chosen **per stratum**, and each event's $q$ comes from **its own stratum's curve**.

#### ⚑ Finding 2 — the scan had to be rewritten to stream, and that made the low floor free

The old `scan()` concatenated every candidate's $\Lambda$, and the observed pass built one dict plus
a cell array per candidate. At floor 2, Pre-TX yields millions of candidates per scan; 1,000
permutations of that is unmaterialisable, and the observed pass alone would have needed GBs.

Now each `(depth, anchor)` block's survivors are histogrammed into the fixed grid and discarded, so
memory is $O(\text{strata}\times|\text{grid}|)$ whatever the floor, and the observed run is
**two-pass**: pass 1 counts only → per-stratum thresholds chosen from the whole curve → pass 2
collects only what clears its own stratum's threshold. ⚑ **Cost is essentially unchanged** — the
prefilter is the same comparison the old code already did, and only survivors are binned:

| | floor 10 | floor 2 |
|---|---|---|
| Mouse 3 hard, s/scan | 1.1 | **1.0** |
| Mouse 3 soft d6, s/scan (floor 4) | 5.3 | 6.5 |

#### Smoke test (Mouse 3, small $B$) — the effect is large

| | floor 10, unstratified | floor 2, stratified |
|---|---|---|
| hard candidates | 8,538 | 39,964 |
| **hard events** | **73** | **419** |
| median clade | 23 cells | **6 cells** |
| threshold | 10 nats (one, arbitrary) | **2.5–4.1 nats, per stratum** |

Per-stratum thresholds are *not* monotone in size (4-5 cells → 4.1 nats; 20-49 → 2.5; 50-199 → 3.2),
which no hand-set rule would have guessed and which is the point of estimating them.

⚠ The soft partial fraction **rises** at the lower floor (Mouse 3 d6: 7.1% → 18.0%), as weaker
events are more often partial. The depth gradient survives (22.5% → 13.7% over depths 1–6), but
**the "graded appearance is very largely clade coarseness" conclusion must be re-read against the
low-floor catalogue** rather than carried over.

#### Configuration

Floors: **hard 2 nats, soft 4 nats.** Soft stops at 4 because script `45` measured
$\Lambda_{\rm soft}$'s null at mean −1.74, sd 3.90 — below ~4 nats there is nothing worth
calibrating. Outputs are tagged `_lam2` / `_lam4`, so **the floor-10 catalogue is not overwritten**
and the tables above it stay reproducible.

300 array tasks, all `RUNNING` immediately, all on `-p cpu` (zero on `lesliec`). Collector
`lf_collect` (11476489) pools and re-runs the observed scans at 64 G — the parts are counts-only
and need 8–12 G, but the observed pass does hold per-candidate arrays at a low floor.

⚠ **Not yet done, and deliberately out of scope here:** the clone-wide layer (`42`) still uses a
single un-stratified floor, and clone sizes vary far more than clade sizes do, so the same
criticism applies to it with more force. Its Mouse 2 result (FDR 4.36%) is the one already-marginal
number in the whole programme. **Stratify `42` by clone size next.**

⚠ **The 5% target is still not decided.** The smoke run's strata all land at 4.5–4.9% FDR, i.e.
right on the target, so at 5% roughly 1 event in 20 is false — ~21 of Mouse 3's 419. The run stores
the **full per-stratum curves**, so any criterion (a stricter rate, or an absolute "expected false
events ≤ 10") can be applied afterwards at zero cost. Decide it when the curves land.

### ⭐ The low-floor run — results, and two corrections to my own reasoning (2026-09-07)

`lf_collect` completed in 40 min: 15 configurations, 300 parts, $B=1{,}000$ each, $p=1/1001$
throughout. Then `52_collect_budget.sh` re-reported all 15 against the **already-stored** nulls under
an absolute criterion. The headline is that **the criterion, not the floor, was doing almost all of
the work** — and that the floor-10 catalogue was closer to right than I claimed.

#### At the 5% target the gains look dramatic. They are mostly false positives.

| arm | events @ 5% | @ floor 10 | **expected FALSE candidates** | events reported |
|---|---|---|---|---|
| Pre-TX | 25,274 | 1,783 | **47,898** | 25,274 |
| Subclone | 24,684 | 8,220 | **35,784** | 24,684 |
| Mouse 1 | 2,083 | 320 | **5,515** | 2,083 |

⚠⚠ **The FDR is computed on *candidates*; the reported quantity is *events*, after the overlap
collapse.** At floor 10 that gap was harmless — Mouse 1 had 2.15 expected false candidates against
320 events, so even total survival of every false candidate left the catalogue clean. At 5% it is
fatal: 5,515 expected false candidates against 2,083 events, and the count-vector design cannot
dedup the permuted sets to close the gap. **A candidate-level rate cannot bound an event-level
error when the rate is not tiny.**

⇒ `--budget X` added to `40`/`43`: each stratum's threshold is the smallest $\Lambda$ whose
**expected null count** is $\le X$. An absolute budget bounds the false *event* count too, because
dedup can only ever reduce it.

#### The defensible catalogue — floor 2, stratified, expected false $\le2$ per stratum

| arm | **events** | floor 10 | ratio | **cell × tape** | floor 10 | % of all missing | floor 10 | exp. false | max $q$ |
|---|---|---|---|---|---|---|---|---|---|
| Subclone | **6,597** | 8,220 | **0.8×** | 280,319 | 266,688 | **20.07%** | 19.09% | 10.5 | 2e-4 |
| Pre-TX | **2,019** | 1,783 | 1.1× | 32,146 | 22,256 | 2.46% | 1.70% | 8.7 | 1e-3 |
| Mouse 2 | **452** | 388 | 1.2× | 23,233 | 21,277 | 7.02% | 6.43% | 10.6 | 1e-3 |
| Mouse 1 | **424** | 320 | 1.3× | 17,620 | 14,770 | 4.34% | 3.64% | 10.5 | 1e-3 |
| Mouse 3 | **132** | 73 | 1.8× | 4,033 | 2,916 | 3.15% | 2.28% | 8.9 | 2e-3 |

⇒ **1.1–1.8× more events than floor 10, and *fewer* for Subclone.** Not 3–14×. Every arm carries
$\le11$ expected false candidates and every event clears $q\le2\times10^{-3}$. This is the catalogue
to quote.

#### ⚠⚠ Correction 1 — "the floor is conservatism, not rigour" was wrong in magnitude

Written 2026-09-04: *"At FDR 0.002–0.019% the events at $\Lambda=10$–12 are almost certainly real;
we are hundreds of times inside the cliff, so the reported counts are loose lower bounds."*

The direction was right, the magnitude badly wrong. **I measured our distance from the
false-positive cliff in FDR units (262–2,511× below target) when the relevant distance is in
$\Lambda$ units.** Under a defensible criterion the per-stratum thresholds land at **5.9–12.4 nats**
— i.e. essentially where the arbitrary floor already sat. The cliff is 0–4 nats below 10, not far
below. So the floor-10 catalogue was *accidentally near-right*, and what this exercise actually
bought is **a justified threshold instead of an inherited one**, plus ~20–80% more events, plus the
knowledge of where the cliff is. That is worth having, but it is not the discovery of a large hidden
population of events.

#### ⚠⚠ Correction 2 — I called the small-clade strata "dead" and they are not

I computed each stratum's maximum achievable $\Lambda$ as $m\log((1-\varepsilon)/\bar p)$ using the
**arm-median** expected rate, concluded the 4–5 cell stratum could not clear its own threshold in
any arm, and said so. **The budget catalogue calls 18–635 events in that stratum in all five arms.**

The error: $\tilde p$ is per (cell, tape), not per arm, and the events that *do* clear are exactly
the ones where recovery was expected:

| arm | arm median $\tilde p$ | $\tilde p$ of called 4–5 events | implied cap | threshold |
|---|---|---|---|---|
| Mouse 1 | 0.405 | **0.098** | 11.6 | 9.6 |
| Pre-TX | 0.347 | **0.093** | 11.8 | 10.9 |
| Subclone | 0.116 | **0.052** | 14.7 | 12.4 |

⇒ **capture-independence is doing the work.** A 4-cell clade is certifiable only when those four
cells were well captured *and* that tape is normally reliable — which lifts its $\Lambda$ ceiling
above the threshold. The right statement is not "small clades are undetectable" but **"a small clade
is detectable only on a tape that should have been there"**, which is the same property that makes
$\Lambda$ a test of inherited loss rather than of dropout rate. Small clades are still heavily
suppressed — 34% of Mouse 1's candidates sit in 4–5 cells against 12% of its called events — but
the stratum is live.

#### ⚑ And the `MAX_D=6` conclusion survives — my worry about it was itself the 5% artefact

I flagged that the soft partial fraction rose at the lower floor (Mouse 3 d6 8.1% → 18.0%) and that
"the graded appearance is very largely clade coarseness" would need re-reading. At the **budget**
threshold it does not:

| arm | partial, budget | partial, floor 10 |
|---|---|---|
| Mouse 1 | **15.9%** | 20.0% |
| Mouse 2 | **22.0%** | 29.6% |
| Mouse 3 | **8.9%** | 8.1% |
| Pre-TX | **15.6%** | 13.7% |
| Subclone | **46.3%** | 50.3% |

Essentially unchanged, and *lower* in three arms. The rise at 5% was the mechanical consequence of
admitting weaker events — a partial loss scores lower $\Lambda$ by construction, so any threshold
drop inflates the partial fraction without any change in biology. **⚠ The partial fraction is
therefore threshold-dependent in a predictable direction and must always be quoted with its
threshold.** The depth gradient also survives (Subclone 82.1% → 35.7% over depths 1–6).

#### Still open

- **`42` (clone-wide) remains unstratified.** Clone sizes vary far more than clade sizes, so the
  criticism applies with more force; its Mouse 2 FDR of 4.36% is the one marginal number left.
  Re-run it with size strata and `--budget`.
- **An event-level FDR** would need the dedup applied to permuted candidate sets, which the count
  vectors cannot support. The absolute budget sidesteps it rather than solving it.

## ⭐⭐ Per-combo nulls: scoring each combo against its own 1,000 permutations, in nats

Justin, 2026-09-07: *"For a given cell × tape combo, why can't we just check where the real cell
labeling falls among the 1000 permutations in terms of nats?"* **We can, it is better than the
stratified global threshold it replaces, and the reason it was not done first is a storage choice
of mine.** `src/53_percombo.py` + `src/54_percombo_merge.py`.

### ⚠ This defeats my own earlier objection, and I answered the wrong question

On 2026-09-03 I argued per-combo scoring is "the wrong tool at this floor" because the p-value
censors at $1/(B+1)$: with a per-combo null tail below $10^{-6}$, every real event returns exactly
$1/1001$ and they cannot be ordered. **That objection is about the *unit*, not the idea.** In nats
it evaporates — "this combo is 25 nats above the maximum of its own 997 permutations" is fully
informative precisely where the rank saturates.

### What it dissolves

| problem in the global-threshold analysis | why it disappears |
|---|---|
| six hand-drawn **clade-size strata** | each combo's own null conditions on its clade size, its own $\tilde p$ values *and* its clone's composition — at full resolution |
| **the $\Lambda$ cap** ("a 9-cell Pre-TX clade cannot exceed 12.4 nats") | the cap applies to that clade's *null* too: topping out at 6 nats and scoring 11 is decisive whatever a global cut says |
| **candidate-vs-event FDR mismatch** | each combo is scored on its own terms, not by a global cut followed by dedup |
| **the arbitrary scan floor** | there is none; every scorable combo is scored |

### Why it was not done first

The permutation parts stored count vectors **aggregated over combos**, discarding combo identity.
That followed the original spec and was right for an aggregate FDR — and wrong for this. A combo
*does* have a stable identity: `(depth, anchor, subclade code, tape)`, all derived from `codes`,
which the permutation never touches. So per-combo accumulators are well defined; we just never kept
them. Measured: 5.6 M–103 M combo slots per arm, four accumulators in float32 = **0.02–1.7 GB**.
⚠ First version allocated $G\times K$ per block and wasted 79% of the array on clades below
`MIN_CLADE` (26.2 M slots for Mouse 3's 5.5 M combos, and ~490 M / 24 GB for Pre-TX). Compacting to
kept clades only gives **99.4% scorable** and 4.87 GB peak on Pre-TX.

### ⚑ Multiplicity, calibrated exactly — hold out permutations

The margin needs its own null *across* combos, and accumulators cannot give one. So `53` holds out
the first `HOLD = 3` permutations as **pseudo-observed** and accumulates over the remaining 997:

$$m_{\rm obs} = \Lambda_{\rm obs} - \max_{b\ge3}\Lambda_b, \qquad
  m_{\rm null} = \Lambda_{b<3} - \max_{b\ge3}\Lambda_b$$

Both are the margin of a candidate over the max of the **same** permutation set, so under $H_0$ they
are *exactly* exchangeable — no parametric null, no approximation, nothing censored. Then
$\mathrm{FDR}(t) = \bigl(\#\{m_{\rm null}\ge t\}/3\bigr)\big/\#\{m_{\rm obs}\ge t\}$, monotonised
BH-style up the grid.

⚠ **Report the margin, not a z-score.** For a small clone the null is genuinely coarse — a 4-cell
clade in a 6-cell clone has only $\binom{6}{4}=15$ distinct compositions, so its 997 permutations
sample 15 values and $\mathrm{sd}$ can be exactly 0. `sd` is written out as context only.

### Mouse 3, validated (2026-09-07) — and the exchangeable null earns its keep immediately

| margin $\ge$ | observed combos | expected false | FDR |
|---|---|---|---|
| 0 nats | 1,243,461 | **1,211,286** | **97.4%** |
| 0.5 | 27,638 | 1,117 | 4.04% |
| 5 | 16,385 | 77.0 | 0.47% |
| 11.37 | 10,063 | **2.0** | **0.020%** |
| 24.8 | 5,611 | 0.0 | 0% |

⚑ **Look at the first row.** 22% of combos have $m_{\rm obs}\ge0$ — vastly more than the $1/998$ a
continuous null would give — because for small clones the same clade composition recurs across
permutations, so $\Lambda_{\rm obs}$ *ties* the maximum constantly. The held-out null reproduces
that exactly (1,211,286 against 1,243,461), so the calibration reports FDR 97% and the ties are
correctly declared meaningless. **A naive per-combo analysis without this null would have called a
million events.**

**At expected false $\le2$: margin $\ge11.37$ nats, 10,063 combos → 102 events, 2.82% of all
missing.** For comparison: 73 events / 2.28% (floor 10, unstratified) and 132 / 3.15% (floor 2,
stratified, budget). So the per-combo route lands **between the two, with no strata, no floor, and
an exactly exchangeable null**.

Two things worth saying out loud:
- **All 102 events beat every one of their 997 permutations** ($n_{\ge}=0$). That is the sentence the
  question was asking for, and it needs no threshold argument.
- **Minimum clade size among called events is 5 cells** — independent confirmation of Correction 2:
  small clades are certifiable, on tapes that should have been there.

⇒ **This should become the primary analysis**, with the stratified global threshold kept as the
cross-check it now is. Remaining arms in flight (~70 min); `pc_merge` collects them.

---

# ⭐ Direction 2, opened 2026-09-07: does the co-integrated symbol vanish?

The first **orthogonal** test of heritable silencing. Everything up to here infers silencing from
*missingness*, which is also what technical dropout looks like. The cassette
`PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC` carries a pegRNA and a tape on **one integration**,
and pegRNAs act in ***trans***, so integration $z$'s symbol $s(z)$ is written into *every* tape in
the cell. If silencing kills the whole locus, then losing tape $z$ from the readout should also
remove symbol $s(z)$ from **every other tape** in those cells — a channel transcript capture cannot
reach.

**Scripts.** `56_symbol_cache.py` (symbol array + first premise check, *superseded* in part) ·
`57_writes.py` (the write unit, shared) · `58_tape_cis.py` (premise check, corrected) ·
`59_depletion.py` (hand-inspected clade examples) · `60_deconvolve.py` (tape → symbol).

## Design decisions, taken with Justin 2026-09-07

| decision | chosen | why |
|---|---|---|
| unit | **independent post-MRCA writes** | deduplicates identity-by-descent exactly; raw entries are pseudoreplicated |
| layer | **sub-clone clades in the large clones**, + clone-wide for Pre-TX | the rest of the clone is an internal control no clone-level effect can reach |
| contrast | **rest of clone, site-matched** | site 6's TV vs pooled is ~2× sites 1–5, and clade writes sit deeper than clone writes |
| premise check | **run first**, on the same write unit | a *cis* component would manufacture the very signal we want to claim |
| alphabet | **design-conforming `NNNNGGA` only** | junk is a per-amplicon parsing artefact, hence intrinsically tape-specific |
| tape → symbol | **regression deconvolution + assignment** | a clade that loses $k$ tapes gives ONE depletion vector |

## Step 0 — the power gate, and it opens

The README's worry was that "tapes are ~4.5–5 of 6 saturated, so most content is ancestral". The
right unit settles it. For cell group $G$, tape $y$, site $j$: take the depth-$(j{-}1)$ prefix as
parent; if its children within $G$ take ≥2 distinct values, the MRCA had site $j$ blank, so every
distinct depth-$j$ code under it is **one write that happened inside $G$**, hence after any loss on
$G$'s stem. Identical prefixes collapse to one write, so identity-by-descent is deduplicated exactly.

| group | cells | $W$ | expected copies of a median-$\xi$ symbol |
|---|---|---|---|
| Mouse2 clone 76 (whole clone) | 3,387 | 40,389 | 300 |
| Mouse2 clone 76, anchor 112 d3 clade | 2,588 | 28,933 | 216 |
| …rest of that clone | 799 | 14,266 | 106 |
| **Subclone clone 6, anchor 90 d1 clade** | 4,967 | **214,076** | **1,238** |
| …rest of that clone | 932 | 78,854 | 456 |
| Mouse2 clone 171 | 102 | 8,750 | 65 |

⇒ a complete depletion inside Subclone's clade is a ~1,200 → 0 contrast. **This contradicted the
README's expectation that the clone-wide layer would be needed for power** — but see Correction 4:
that was right for the *large* clones and wrong as a general statement.

⚠ Cost of the dedup, stated: two independent writes of the same symbol under the same parent
collapse to one code, deflating frequent symbols. It biases inside and outside identically, so a
*contrast* survives it; an absolute $\xi$ would not.

## ⚠⚠ Correction 1 — my first premise check was pseudoreplicated, and I read it as a result

`56` counted raw (cell, tape, site) entries and took its noise floor from a random split of those
entries. Both halves are wrong for one reason: a symbol written at tape $y$ in a clone founder
reappears in every descendant. It reported median per-tape TV 0.2745 against a floor of 0.0828 with
144/158 tapes "above 2× floor" — which reads as a large *cis* effect and is not one.

**The diagnosis is in the cross-arm ordering**, which tracks clonal concentration and nothing about
tapes: Pre-TX 0.0398 (many clones, median 7 cells) · Mouse1 0.2494 · Mouse3 0.2745 · Subclone 0.4210
· **Mouse2 0.5632** (one clone holding 98.9% of the pooled weight). The per-*site* statistic moved
the same way (Mouse2 s1 = 0.192 vs Pre-TX s1 = 0.020), which settles it — **a site cannot be *cis*
to a tape.** The analytic floor confirms the arithmetic: $\tfrac12\sqrt{2/\pi n}\sum_s\sqrt{p_s}=0.073$
at $n=2{,}860$, matching the observed 0.083 — the right formula for the wrong $n$.

## The corrected premise check (`58`) — small, real, and *predicted by the mechanism*

Same contrast on deduplicated writes, clone-stratified; null permutes the **tape label** within
(clone, site), preserving each clone's composition, each site's composition and every (clone, tape)
write count.

| arm | TV obs | TV null | excess | tapes > null p95 | $G$ vs null mean | $p$ |
|---|---|---|---|---|---|---|
| Mouse3 | 0.1513 | 0.1468 | +0.0035 | 15/160 | +2.6% | 0.005 |
| Mouse2 | 0.1352 | 0.1250 | +0.0064 | 32/164 | +6.3% | 0.005 |
| Mouse1 | 0.0927 | 0.0804 | +0.0097 | 73/165 | +12.1% | 0.005 |
| Subclone | 0.0473 | 0.0413 | +0.0066 | 101/166 | +16.2% | 0.005 |

⇒ pseudoreplication was ~98% of the apparent effect, but a residual survives at $p=1/201$ in every
arm and its **relative size orders by statistical power** (2.6% → 16.2%), which is what a real
effect looks like.

⭐ **The residual is what the mechanism predicts.** Tape $t$'s writes come only from cells that
*recovered* $t$ — cells where integration $t$ is expressed — and in exactly those cells pegRNA $t$ is
expressed too, so $s(t)$ is over-supplied. So the premise check is **not separable from the
hypothesis**: it is the same measurement, run clade-free on every write in the arm. Its argmax is a
candidate $z\mapsto s(z)$ map at the highest power available.

## ⚠⚠ Correction 2 — junk symbols dominated the per-tape signal

Top per-tape hits on the powered arms were `TGGGA` (Mouse2, $z=67$), `GAATGGATGAT` (Subclone,
$z=99$), `ACAAGGG`, `GACTCGGTGCC`, `GCCGGGGTGAG`. **None is a pegRNA insert.** `xi_vectors.json`
carries 6 non-conforming "promoted" symbols above the count threshold, 2.2% of edits — and junk is a
*parsing* artefact of a particular amplicon, therefore intrinsically tape-specific. A third reading
of the residual that was not on my list.

⇒ composition is now measured on **design-conforming `NNNNGGA` only** (94–98 symbols per arm).
⚠ This does **not** contradict the earlier correction that junk values *are* heritable characters.
Junk is a fine lineage character — it is inherited. It is not evidence about the writing pool. **The
two uses need different alphabets**, and conflating them is what produced the first ranking.

## ⚠⚠ Correction 3 — I declared collapse bias absent on the basis of the weakest arm

The write dedup merges two independent writes of the same symbol under one parent, so bushier tries
collapse more and collapse hits *frequent* symbols hardest. I checked this on Mouse3, found
$\mathrm{corr}(\text{bushiness},\text{excess})=+0.125$ and $\mathrm{corr}(\xi_s,\overline{|z|})=-0.150$,
and said collapse was not driving it. On the powered arms:

| arm | corr(bushiness, excess) | corr($\xi_s$, mean\|z\|) | top $z$ median |
|---|---|---|---|
| Mouse3 | +0.125 | −0.150 | 2.99 |
| Mouse2 | −0.006 | +0.156 | 3.25 |
| Subclone | −0.196 | **+0.356** | 4.11 |
| Mouse1 | **−0.767** | **+0.363** | 6.74 |

The positive correlation with symbol frequency is exactly the collapse signature. **Mouse3 simply
lacked the power to show it**, and a null result on the weakest arm is not a clearance.

## The clade examples (`59`) — the signal is real where power exists

Statistic: per symbol, a site-stratified common log-odds ratio $\theta_s$ fitted by profile
likelihood, reported as $\Lambda_s$ in **nats**. Null: whole cell rows permuted within clone,
**writes re-extracted** (which writes are post-MRCA depends on who is in the clade).

| arm | clade | $n$ | $W_{\rm in}$ | top symbol | $\theta$ | $\Lambda$ | null max | runner-up |
|---|---|---|---|---|---|---|---|---|
| Mouse2 | clone 76 / anchor 125 d3 | 1,304 | 19,304 | AAGCGGA | −0.85 | **56.2** | 6.1 | 7.1 |
| Subclone | clone 6 / anchor 90 d1 | 4,967 | 197,013 | TGGCGGA | −2.80 | **85.4** | 7.7 | 36.2 |
| Mouse1 | clone 36 / anchor 26 d2 | 1,131 | 43,050 | GAAAGGA | −0.75 | **26.1** | 5.5 | 2.3 |
| Mouse3 | clone 111 / anchor 48 d1 | 157 | 1,984 | AATGGGA | −1.35 | 5.2 | 3.0 | 3.0 |
| Pre-TX | clone 676 / anchor 69 d2 | 71 | 5,452 | GTAGGGA | −1.25 | 7.3 | 3.9 | 2.1 |

All design-conforming, all far past the permutation null, and Mouse1's is a textbook single outlier
(26.1 against a runner-up of 2.3). **Different clades within Mouse2 return different symbols**, so it
is not a clone-wide artefact.

## ⚠⚠ The structural limit — the statistic identifies CLADES, not TAPES

$\Lambda$ depends on the clade; the lost tape enters only by being excluded from the write pool. So
**a clade that loses $k$ integrations yields one depletion vector with all $k$ superposed.**

- Subclone's three catalogue "events" (tapes 65, 157, 142) are **one clade scored three times** —
  identical $\theta=-2.80$, identical $\Lambda\approx86$. They are one measurement.
- That clade shows a **ladder** of six symbols above its null max: 85, 36, 24, 21, 13, 9 — consistent
  with several integrations silenced together.
- Mouse1's two depth-2 rows are likewise one clade.

⇒ 19 catalogue rows are ~15 distinct clades, and **reading a map off single clades is impossible in
principle, not merely underpowered.** Hence `60`.

## ⚠⚠ Correction 4 — Pre-TX needs the clone-wide layer after all

Its clades run 71–95 cells and **4 of 5 returned nothing above the null** despite $W_{\rm in}$ of
5,000–7,000: the limit is clade size, not writes. The README's original recommendation to screen the
clone-wide layer first was right for this arm; my Step 0 power table generalised from Subclone and
Mouse2, where it holds, to arms where it does not.

## `60` — the deconvolution

Unit $g$ = a scored group; $X[g,z]=1$ if tape $z$ is called lost in $g$; for each symbol,
$y[g,s] \approx \sum_z X[g,z]\,B[z,s]$ by weighted least squares, $y$ = site-stratified
Mantel–Haenszel log-OR, weight = its RBG inverse variance. Under the hypothesis $B[z,\cdot]$ is
nonzero at one symbol, so $\arg\min_s Z[z,s]$ **is** the map.

⚠ **Losses cluster** ($\rho_{\rm tape}=0.25$), so tapes always lost together are not separable. The
script reports the conditioning rather than hiding it behind the ridge: effective rank of $X$, and a
per-tape VIF; tapes with VIF ≥ 5 or coverage < 5 are marked uncallable and excluded from concordance.

⚠⚠ **A correctness fix made while building it.** For the clone-wide layer the comparison group must
be **pooled from per-clone extractions**, never extracted as one group: handing `extract_writes`
every cell outside clone $c$ makes it assess polymorphism *across* clones, where almost every parent
has several children, so the dedup collapses and nearly every entry counts as a write — the exact
pseudoreplication this direction was rebuilt to avoid. The contrast is assembled by arithmetic from
$T$, $T_c$, $T_z$ and the $(c,z)$ cross terms.

**Mouse3 (smallest arm) is underdetermined and says so:** 155 usable units for 166 tapes, 89 fitted
at coverage ≥ 3, 42 callable, median $z=-2.21$ with a median gap to the runner-up of **0.25** — i.e.
noise. Powered arms in flight.

## Where this stands, and what is still owed

**Settled:** the write unit and its power; that the premise check's first version was
pseudoreplicated; that a small tape-dependence is real, is predicted by the mechanism, and is partly
collapse bias and partly junk; that the clade-level depletion signal is real and large in the three
big clones; that clades, not tapes, are what a single contrast identifies.

**Not yet done, and needed before any claim:**
1. the recovered map from `60` on the powered arms, with its conditioning;
2. **five-way cross-arm concordance** — the killer test, since all 166 TargetBCs are identical
   across arms, so agreement by chance among ~98 conforming symbols is ~1%;
3. injectivity and the collision rate against the $E=122$ distinct predicted for 166 draws from 256;
4. dose ($j$ integrations naming a symbol ⇒ ~$j\times$ base $\xi$) and $\mathrm{corr}(\beta_z,\xi_{s(z)})>0$;
5. a permutation null on the concordance itself (permute loss sets across units, refit, re-measure).

## The recovered map (`60`) — and ⚠⚠ Correction 5, the arms are not independent

| arm | units | tapes fitted | callable | distinct symbols | median $z$ | median gap |
|---|---|---|---|---|---|---|
| Pre-TX | 1,449 | 166 | **160** | 71 | **−4.05** | **1.32** |
| Mouse1 | 585 | 145 | 105 | 48 | −2.89 | 0.42 |
| Mouse2 | 348 | 118 | 70 | 36 | −2.66 | 0.38 |
| Mouse3 | 155 | 89 | 42 | 29 | −2.21 | 0.25 |

Injectivity: distinct/expected-if-injective = 0.87 (Pre-TX), 0.85 (Mouse3), 0.74 (Mouse1), 0.72
(Mouse2) — below the collision prediction, i.e. some pile-up. ⚠ Mouse2 sends **15 of 70 tapes to
`AAGCGGA`**, which is what a near rank-1 design matrix produces: one clone holds 98.9% of that arm,
so its clades are nested and share their loss sets, and the deconvolution cannot separate them. VIF
< 5 did not catch it — **the distinct-symbol count is the diagnostic that does.**

Cross-arm concordance, permutation null ($B=2{,}000$, permuting tape labels within one arm):

| pair | $n$ | agree | rate | null mean | null max | × | $p$ | **shared clones** |
|---|---|---|---|---|---|---|---|---|
| Pre-TX × Mouse1 | 105 | 17 | 0.162 | 2.08 | 8 | **8.18** | 0.0005 | **96.3%** |
| Pre-TX × Mouse3 | 42 | 5 | 0.119 | 0.73 | 5 | 6.82 | 0.0015 | 95.3% |
| Mouse1 × Mouse3 | 40 | 3 | 0.075 | 0.55 | 4 | 5.43 | 0.0125 | 42.0% |
| Mouse1 × Mouse2 | 51 | 2 | 0.039 | 0.50 | 3 | 4.00 | 0.0905 | 42.7% |
| Pre-TX × Mouse2 | 70 | 3 | 0.043 | 0.96 | 6 | 3.12 | 0.0755 | 96.3% |
| Mouse2 × Mouse3 | 25 | 1 | 0.040 | 0.58 | 4 | 1.74 | 0.4553 | 31.3% |

restricted to confident calls ($z\le-3$, gap $\ge0.5$): **Pre-TX × Mouse1 13/27 = 48.1%, 17.0× null,
null max 4, $p=0.0005$.**

⚠⚠ **Correction 5, and it retracts the headline I was about to write.** The README has treated
five-way cross-arm concordance as "essentially unfakeable" because all 166 TargetBCs are identical
across arms. They are — but **the arms are not independent samples.** Measured on `ClonalBC`:

| pair | shared clones | % of the smaller arm |
|---|---|---|
| Pre-TX × Mouse1 | 285 | **96.3** |
| Pre-TX × Mouse2 | 210 | 96.3 |
| Pre-TX × Mouse3 | 143 | 95.3 |
| Mouse1 × Mouse2 | 93 | 42.7 |
| Mouse2 × Mouse3 | 47 | 31.3 |
| Subclone × mice | 1 / 4 / 1 | 6.7 / 26.7 / 6.7 |

The mice were transplanted from the pre-TX pool. Silencing is heritable and predates
transplantation, so **a clone measured in two arms carries the same losses** — cross-arm agreement is
the same events re-measured, not a replication. My strongest number, Pre-TX × Mouse1, is also the
most confounded pair.

⇒ **The clean test is a within-arm split on disjoint clones** (`60 --split h/N`): two halves share no
cell, no clone and no lineage, so agreement between them tests whether $z\mapsto s(z)$ is a property
of the *line* rather than of the clones it was fitted on. Launched on Pre-TX and Mouse1.
Subclone × mouse is the one honest cross-arm pair (≤ 4 shared clones) — Subclone's `60` run is still
in flight.

## ⚠⚠ Correction 6 — the cell-level route fails its own validation

I proposed promoting `58`'s cell-level scan to a primary route, on the argument that the residual
tape-dependence is what the mechanism predicts. It does not reproduce:

| test | deconvolution (`60`) | cell-level (`58`) |
|---|---|---|
| best cross-arm enrichment | 8.18× ($p=0.0005$) | 3.59× ($p=0.011$) |
| median cross-arm enrichment | 5.4× | 1.5× |
| pairs at $p<0.05$ | 4 of 6 | 3 of 10 |
| cross-method agreement (within arm) | Pre-TX 3.38× ($p=0.006$); Mouse1, Mouse2, Mouse3 **0×** | |

The two methods **disagree within the same arm** in three of four cases. So the cell-level residual
is real (p = 1/201 everywhere) but is not the map — it is dominated by collapse and capture
structure, exactly the two nuisances the shape diagnostic was meant to separate and did not. The
argument that it *should* carry the signal was sound; the measurement says it does not, and the
measurement wins.

## ⭐⭐ The clean test — a within-arm split on disjoint clones — PASSES on Pre-TX

`60 --split h/2` fits the map on a random half of the arm's **clones**. The two halves share no
cell, no clone and no lineage; the only thing they have in common is the engineered line — the
integrations and their co-integrated pegRNA barcodes. Verified disjoint: Pre-TX 1,473 + 1,473 of
2,946 clones (1,431 + 1,343 of 2,774 units); Mouse1 148 + 147 of 295 (235 + 223 of 458).

| split test | $n$ | agree | rate | null mean | null max | × | $p$ |
|---|---|---|---|---|---|---|---|
| **Pre-TX half0 × half1, confident calls** | 55 | **48** | **0.873** | 1.42 | 7 | **33.9** | **0.0005** |
| Pre-TX half0 × half1, all calls | 151 | 58 | 0.384 | 2.33 | 8 | 24.9 | 0.0005 |
| Mouse1 half0 × half1, all calls | 23 | 1 | 0.043 | 0.49 | 4 | 2.04 | 0.397 |

⭐ **87.3% of confidently-called tapes name the same symbol when estimated from disjoint halves of
the clone pool**, against a permutation null whose best of 2,000 draws reached 7 of 55. The five
strongest calls are identical in both halves — `TGTAGCCGGC→AATCGGA` ($z=-17.7/-14.7$),
`GAACGACTTC→ATACGGA`, `TGTTCGATGC→GACCGGA`, `GGAAGTGGTA→GTGCGGA`, `TATCTTCGGT→ACGTGGA`.

⚠ The permutation null permutes tape labels **within** a half, so it preserves each half's marginal
symbol distribution exactly — pile-ups on common symbols cannot generate this.

**Mouse1 fails its own split, and that is a power statement, not a contradiction.** Its halves have
220–234 units against Pre-TX's ~1,400, the design goes rank-deficient (rank 94/95, cond $7\times10^{12}$),
and half1 sends 8 of its top 10 tapes to `ATGTGGA` — the same near-rank-1 pile-up as Mouse2's
`AAGCGGA`. ⇒ **the deconvolution needs many clones with *different* loss sets, which is exactly what
Pre-TX has and the mice do not** — 2,946 clones of median 7 cells, where the mice are dominated by
one or two huge clones. That inverts the usual ordering in this project, where Pre-TX has been the
underpowered arm for everything lineage-related.

⇒ **Quote the split, not the cross-arm numbers.** Pre-TX × Mouse1 at 8.18× shares 96.3% of its
clones; the split at 33.9× shares none.

## ⚑ What the deconvolution actually needs: many INDEPENDENT clones, not many clades

Subclone finished last and makes the point sharply. It has the **most units of any arm** (1,673),
full rank (166/166), excellent conditioning (cond 10.2, VIF max 2.2) and the best median $z$
($-4.56$) — and its map is the **most degenerate of all**: 8 of its top 10 tapes call `ATGTGGA`, and
distinct/expected-if-injective is **0.58**, the worst score in the table.

The reason is that its 1,673 units are clades drawn from **15 clones**. Clades nested inside one
clone share almost all of their loss sets, so they are near-duplicate rows of $X$ — which inflates
the apparent unit count and the rank without adding information about *which* tape carries a
depletion. The conditioning diagnostics do not see it, because $X$ is full rank; only the
distinct-symbol ratio does.

| arm | clones | units | distinct/expected | split reproducibility |
|---|---|---|---|---|
| **Pre-TX** | **2,946** | 2,774 | **0.87** | **87.3% (33.9×, $p=0.0005$)** |
| Mouse3 | 150 | 156 | 0.85 | — (underdetermined) |
| Mouse1 | 295 | 458 | 0.74 | 4.3% (2.0×, $p=0.40$) |
| Mouse2 | 218 | 348 | 0.72 | — |
| Subclone | **15** | **1,673** | **0.58** | — |

⇒ **the distinct-symbol ratio, not the rank or the VIF, is the diagnostic for whether a
deconvolution has separated the tapes**, and the ordering it gives is the ordering of *clone counts*,
not unit counts. Pre-TX — the arm that has been underpowered for every lineage question in this
project, because its clones are tiny — is the only arm with enough independent loss sets to identify
the map, and it is the arm where the split test passes.

## Status of Direction 2 at the end of 2026-09-07

**Established.**
- The write unit, and Step 0's power (`57`).
- Clade-level depletion is real and large in the three big clones (`59`): Mouse2 $\Lambda=56.3$ vs
  null max 5.7; Subclone 86.2 vs 8.9; Mouse1 26.1 vs 5.5 — all design-conforming symbols.
- ⭐ **The recovered $z\mapsto s(z)$ map reproduces across disjoint clone halves of Pre-TX at 87.3%
  on confident calls, 33.9× a permutation null ($p=1/2001$).** Two halves sharing no cell, no clone
  and no lineage agree on which symbol a tape's silencing removes. **This is the orthogonal
  confirmation the direction was opened to get.**

**Six corrections to my own reasoning, all recorded above:** the first premise check was
pseudoreplicated (1); junk symbols dominated the per-tape signal and needed a separate alphabet from
the lineage-character alphabet (2); I cleared collapse bias from the weakest arm (3); Pre-TX needs
the clone-wide layer after all (4); the arms share up to 96.3% of their clones, so cross-arm
concordance is *not* the killer test the notes call it (5); and the cell-level route I argued for
fails its own validation, disagreeing with the deconvolution in 4 of 5 arms (6).

**Still owed.**
1. Dose — a symbol named by $j$ integrations should carry ~$j\times$ base $\xi$.
2. $\mathrm{corr}(\beta_z,\xi_{s(z)})>0$, the fig-3b link.
3. Collision structure of the Pre-TX map against $E=122$ for 166 draws from 256.
4. Whether the map recovered from Pre-TX **predicts** the depleted symbol in the three big mouse
   clades of `59` — a genuine out-of-sample test, and the one that would tie the two layers together.
5. More split replicates (different seeds, $N=3,4$) to put an interval on the 87.3%.

### Per-combo results, all five arms (2026-09-07) — and a division of labour

**⚠⚠ Correction to the first version of this analysis.** It used ONE GLOBAL margin threshold and I
claimed the per-combo null "dissolves" the clade-size problem. **Half right.** The *null* conditions
on clade size; a *global threshold on the margin* does not. Measured on Mouse 1: clades of 4–9 cells
are **55% of scorable combos and 0.0% of called events**, while 200+ cell clades are 1.3% of combos
and **36.7%** of events. The margin is calibrated **per size stratum** now — the synthesis rather
than a retreat, since the per-combo null conditions on the individual clade in a way no six-bin
scheme can, and per-stratum calibration makes the threshold comparable across sizes in a way a
single global cut cannot. After the fix the **minimum called clade is 4 cells in all five arms.**

| arm | scorable combos | above threshold | **events** | cell × tape | % of all missing | exp. false | $n_{\ge}=0$ |
|---|---|---|---|---|---|---|---|
| Subclone | 37.8 M | 411,727 | **4,405** | 686,803 | 49.17% | 9.7 | 4,405 |
| Pre-TX | 102.6 M | 65,302 | **1,915** | 30,938 | 2.36% | 9.0 | 1,915 |
| Mouse 1 | 14.7 M | 25,580 | **347** | 20,739 | 5.11% | 11.7 | 347 |
| Mouse 2 | 7.1 M | 18,446 | **316** | 33,542 | 10.13% | 10.7 | 316 |
| Mouse 3 | 5.5 M | 11,820 | **128** | 3,952 | 3.09% | 6.7 | 128 |

⭐ **Every called event, in every arm, beats every one of its own 997 permutations.** That is the
answer to the question as posed, and it needs no threshold argument at all.

#### ⚠⚠ But the raw cell × tape share is inflated by COARSE ATTRIBUTION, and Subclone shows it

Subclone reads 49.17% here against 20.07% from the $\Lambda$ route. The two routes are **not
disagreeing about which tape was lost in which clone** — 1,269 of ~1,400 (clone, tape) pairs are
shared, 109 per-combo-only, 53 $\Lambda$-only. They disagree about **which clade to blame**:
per-combo picks a clade **1.67× larger** (median; larger in 76% of shared pairs). And that has a
testable consequence — the inside rate, $n_{\rm missing}/|{\rm clade}|$, which must be 1.00 for a
genuine Dollo loss on the clade's *stem*:

| arm | per-combo (10% / median / mean) | $\Lambda$ route (10% / median / mean) |
|---|---|---|
| Subclone | 0.264 / **0.778** / 0.694 | 0.896 / **1.000** / 0.969 |
| Mouse 2 | 0.361 / 0.888 / 0.773 | 0.939 / 1.000 / 0.984 |
| Mouse 1 | 0.313 / 0.960 / 0.809 | 0.977 / 1.000 / 0.993 |

⇒ the per-combo route is calling **large, partly-missing clades**. Their present cells *disprove* a
stem loss on that clade, so counting all their missing entries as "inside one inherited loss"
over-attributes. **Why it happens:** for a large clade the null $\Lambda$ is strongly *negative* —
a random relabelling is heavily penalised by the present cells it picks up — so even a modestly
positive $\Lambda$ yields a large margin. **The margin is an excellent test of non-exchangeability
and a poor estimator of which clade the loss sits on.** $\Lambda_{\rm hard}$'s high absolute
threshold implicitly demands a high inside rate, and that is what localises the event.

#### ⇒ The division of labour, and it resolves the confusion

| question | statistic |
|---|---|
| **is dropout exchangeable within clones?** (significance) | **per-combo margin** — exact null, no strata, no floor, nothing censored |
| **how many cell × tape entries sit inside an inherited loss?** (attribution) | **$\Lambda$ with a completeness bar** — or the per-combo margin *plus* inside rate $\ge0.9$ |

With the completeness bar applied to the per-combo events, all four routes agree to within ~1.4× in
every arm:

| arm | per-combo, all | **per-combo, inside rate $\ge0.9$** | floor 2 stratified | floor 10 |
|---|---|---|---|---|
| Subclone | 4,405 / 49.17% | **1,498 / 27.57%** | 6,597 / 20.07% | 8,220 / 19.09% |
| Mouse 2 | 316 / 10.13% | **152 / 4.73%** | 452 / 7.02% | 388 / 6.43% |
| Mouse 1 | 347 / 5.11% | **204 / 3.03%** | 424 / 4.34% | 320 / 3.64% |
| Mouse 3 | 128 / 3.09% | **70 / 2.28%** | 132 / 3.15% | 73 / 2.28% |
| Pre-TX | 1,915 / 2.36% | **1,202 / 1.63%** | 2,019 / 2.46% | 1,783 / 1.70% |

⚠ **Correction to what I said when launching this:** I wrote that the per-combo route "should become
the primary analysis" with the $\Lambda$ threshold "kept as the cross-check it now is." Wrong as
stated. They are not competing estimates of one quantity — **they answer different questions**, and
the per-combo route needs a completeness criterion bolted on before it can count anything.

#### The tie diagnostic, which validates the held-out null

Share of combos with margin $\ge0$, against the $1/998$ a continuous null would give:

| arm | observed share | FDR at margin $\ge0$ |
|---|---|---|
| Mouse 3 | 22.4% | 97.4% |
| Pre-TX | 23.9% | 92.6% |
| Mouse 1 | 16.0% | 94.2% |
| **Subclone** | **3.4%** | **6.4%** |

Small clones admit few distinct clade compositions, so $\Lambda_{\rm obs}$ *ties* the null maximum
constantly; Subclone's large clones give a near-continuous null and few ties. The held-out
permutations reproduce each arm's tie rate and correctly declare margin $\ge0$ meaningless.
**Without them a naive per-combo analysis would have reported 24.5 million events in Pre-TX.**

---

# ⭐⭐ The detector, rebuilt: nested families, the score test, and $\hat\pi$ made non-circular (2026-09-08/09)

**Scripts.** `62_percombo_soft.py` (scan) · `63_percombo_soft_merge.py` (calibrate, call, characterise) ·
`64_submit_percombo_soft.sh` (driver) · `65_attrib_diag.py` (attribution diagnostic) ·
`66_fig5_plane.py` (figure).

## Why this was opened

Justin: *the figure showing dropout is heritable needs settling, and what we have never settled is
how to define a silencing event.* Two claims had been running together:

- **Claim A, heritability.** Dropout is not exchangeable among cells within a clone. Established three
  independent ways already, none needing an event definition at all.
- **Claim B, silencing events.** The losses are discrete, complete, heritable — Dollo characters.

⚑ **$\hat\pi$ and `inside_rate` are the same quantity** ($k/m$), so the "completeness bar" bolted onto
the per-combo route on 2026-09-07 already *was* a $\hat\pi$ threshold. The hard and soft routes differ
in exactly one respect: whether completeness is fused into the test statistic or reported separately.

## The identity that dissolves hard-vs-soft

With $U=\sum_{c\in S}[X\log\tilde p+(1-X)\log(1-\tilde p)]$ common to both,

$$\Lambda_{\rm hard} = \Lambda_{\rm soft} - m\,\mathrm{KL}(\hat\pi\,\|\,1-\varepsilon)$$

Verified to $<5\times10^{-13}$ over seven cases spanning $m=9$ to $m=4{,}967$, reproducing `43`'s
docstring examples exactly ($m{=}20,k{=}16$: 2.109 and 10.682, difference 8.5734 $=m\mathrm{KL}$).

$\Lambda_{\rm hard}$ is $\Lambda_{\rm soft}$ minus a penalty for the loss not being total, scaling with
clade size. This explains: why $\Lambda_{\rm soft}\ge\Lambda_{\rm hard}$ always; why the partial fraction
moves with threshold in a predictable direction; why $\Lambda_{\rm hard}$ "implicitly demands a high
inside rate" (it literally pays $m\mathrm{KL}$ to admit anything else); and why $\Lambda_{\rm hard}$ is a
good *attributor* — for $\hat\pi\approx1$ it grows $\approx m(\log(1-\varepsilon)-\log\tilde p)$,
linearly in $m$, while stepping up to a parent with present cells costs ~4 nats each, so it is
maximised at the **largest clade that is still complete**, which is the Dollo rule.

## ⚠⚠ CORRECTION — $\Lambda_{\rm soft}$ is NOT a completeness-agnostic $\Lambda_{\rm hard}$, and the Mouse3 run proved it

I proposed detecting on $\Lambda_{\rm soft}$ and reporting $\hat\pi$. **Wrong.** The first Mouse3 run
returned 340 events / 6.33% of missing, median $\hat\pi$ **0.203** against expected 0.182, 67.4% partial.
Decomposed:

| group | n | $m$ | $k$ | $E$ | $k-E$ | $z$ | $\Lambda_{\rm soft}$ | $\Lambda_{\rm hard}$ | margin |
|---|---|---|---|---|---|---|---|---|---|
| clades 4–49 | 114 | 16 | 12.0 | 5.6 | **+6.47** | **4.00** | +10.1 | +9.5 | 5.88 |
| clades 50–199 | 226 | 134 | 17.0 | 16.9 | **+0.27** | **0.08** | **−9.6** | −478.4 | 4.62 |

**203 of the 226 large-clade calls have $z<0.5$ — no rate elevation at all.**

**The mechanism.** $H_0$ (a product of *different* Bernoullis, $\tilde p_c=\sigma(\alpha_c+\beta_z+\gamma_{C,z})$)
is **not nested** in $H_1^{\rm soft}$ (a product of *identical* ones) unless every $\tilde p_c$ in the
clade is equal — which it never is, since $\alpha_c$ varying is the whole point of the null model. The
two therefore differ in **two** respects at once, the level *and* the homogeneity, and

$$\Lambda_{\rm soft}=\underbrace{m\,\mathrm{KL}(\hat\pi\|\bar e)}_{T_{\rm elev}}+\underbrace{\ell(\bar e\mathbf 1)-\ell(\tilde p)}_{T_{\rm misfit}},\qquad \mathbb{E}_{H_0}[T_{\rm misfit}]=-\sum_c\mathrm{KL}(\tilde p_c\|\bar e)\le0$$

Decomposition verified exact to 2e−14; $\mathbb{E}_{H_0}[T_{\rm misfit}]$ verified against
$-\sum_c\mathrm{KL}$ by 200,000-draw Monte Carlo (−8.2716 vs −8.2527 at $m$=60; −39.4279 vs −39.4182 at
$m$=400).

$T_{\rm misfit}$ is the **dispersion of $\tilde p$ inside the clade**, negative on average, growing with
$m$ and with the spread. Real clades are *more homogeneous in capture quality* than a random subset of
their clone (lineage correlates with capture, $\rho_{\rm cell}=0.13$, and clade membership requires
reaching depth $d$ on the anchor, which itself selects on capture). So
$T_{\rm misfit}^{\rm obs}>\max_b T_{\rm misfit}^{(b)}$ and a positive margin appears with
$T_{\rm elev}=0$. On the 226 artefacts $T_{\rm elev}$ carried **0.1%** of the median margin.

⚠ **This was NOT a calibration failure.** The permutation null of $\Lambda_{\rm soft}$ was exactly right
(54.5 expected false against 26,356 observed). The test correctly rejected exchangeability. What was
wrong is the *direction of the alternative* — a mixture, so rejection is uninformative about which
respect was violated. $\Lambda_{\rm hard}$ escaped only by brute force: its $-m\mathrm{KL}(\hat\pi\|1-\varepsilon)\approx-3.3m$
term put those clades at −490 nats. **Right by accident, and the same term is why it is blind to
graded losses.**

⚠⚠ **A second correction, one hour earlier.** I first blamed the *overlap collapse* — 63 ordered it by
the detection margin rather than by $\Lambda_{\rm hard}$. `65_attrib_diag.py` measured it: collapsing by
$\Lambda_{\rm hard}$ gives 341 events against 340, median $\hat\pi$ 0.210 against 0.203, and only 51 of
338 (clone,tape) groups even contain both a partial and a complete combo. **Coarse attribution was ~15%
of the problem, not the cause.** The final `63` collapses by $\Lambda_{\rm hard}$ anyway — it is the
right rule — but it changes nothing (127 events either way).

## The fix: a NESTED one-parameter family, and the score test

$$H_1(\delta):\ X_c\sim\mathrm{Bern}\big(\sigma(\eta_c+\delta)\big),\qquad \eta_c=\mathrm{logit}(\tilde p_c)=\alpha_c+\beta_z+\gamma_{C,z}$$

$\delta$ is **one more additive term in the same logistic regression** — a clade-specific intercept
(no "logit of a logit": $\sigma$ inverts logit, and $\eta$ *is* the linear predictor we already had).
$\delta=0$ **is** $H_0$, in the interior; $\delta\to+\infty$ is total loss; the per-cell heterogeneity is
kept in **both** models, so $T_{\rm misfit}$ cannot arise. The log-odds scale is forced, not chosen:
it keeps every $q_c\in(0,1)$ for every $\delta\in\mathbb{R}$ (a probability-scale shift escapes the
interval, so $\delta=0$ would not be interior); it is the Bernoulli's natural parameter, which is what
makes $\partial\ell_c/\partial\delta$ collapse to $X_c-q_c$; and $e^\delta$ is an odds ratio, the
natural "the same thing happened to every cell" model.

$$S=\left.\tfrac{\partial\ell}{\partial\delta}\right|_0=\sum_c(X_c-\tilde p_c)=k-E,\qquad
I=-\left.\tfrac{\partial^2\ell}{\partial\delta^2}\right|_0=\sum_c\tilde p_c(1-\tilde p_c)=V$$

$$z=S/\sqrt I=(k-E)/\sqrt V \quad=\quad \frac{\text{observed}-\text{expected}}{\text{SD of the count}}$$

$\mathbb{E}_{H_0}[S]=0$ and $\mathrm{Var}_{H_0}[S]=V$ hold **exactly**, any $m$, any heterogeneity —
verified by Monte Carlo ($\mathbb{E}[z]$ 0.0011/0.0002/−0.0017, $\mathrm{Var}[z]$ 1.0034/1.0004/1.0010 at
$m=8/60/400$). The heterogeneity enters the **scale**, never the **location**. And since $p(1-p)$ is
concave, a homogeneous clade has the *largest* $V$ hence the *smallest* $|z|$ — the artefact direction
is suppressed, not amplified.

⚑ **Three tests, and why the score one.** Nesting makes LRT / score / Wald all available (height climbed
/ slope at the null / horizontal distance). Worked example, built to match the two Mouse3 populations:

| clade | $m$ | $k$ | $E$ | $V$ | $\hat\pi$ | $z$ | $\hat\delta$ | LRT | **Wald** | $\Lambda_{\rm soft}$ | $\Lambda_{\rm hard}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A small, complete loss | 16 | 16 | 5.94 | 3.57 | 1.000 | **5.32** | 32.05 | **33.20** | **0.000** | 16.60 | 16.44 |
| B large, count as predicted | 134 | 17 | 17.69 | 13.83 | 0.127 | **−0.19** | −0.05 | **0.04** | 0.034 | −2.13 | −490.14 |

⚠ **Wald = 0.000 on a complete loss** — complete separation sends $\hat\delta\to\infty$ and
$I(\hat\delta)=\sum q(1-q)\to0$. Never use Wald here. ⚠ **$\hat\delta$ diverges too**, so it is the
effect size for *graded* shifts and $\hat\pi$ is the readout for complete ones — the divergence is the
signature of a Dollo loss, not a numerical failure.

⚠ **Cost, not theory, is why the score screens and the LRT characterises.** The score is three sums, no
iteration. The LRT needs $\hat\delta$ per combo, ~10 bincounts instead of 1 — Mouse3 9→90 min, Pre-TX
~2h→~20h *per part*. So: score for the $10^8$-combo scan, exact LRT + $\hat\delta$ on the called
shortlist ($10^2$–$10^3$ combos, seconds).

⚠ **$\chi^2_1$ does not replace the permutation.** A nominal 5% test using $\chi^2_1$ actually rejects
**9.07%** at $m=8$ (only 9 possible values of $k$; the quantiles plateau on a lattice), 4.60% at $m=16$,
5.33% at $m=134$. Small clades are the bulk (4–9 cells are 55% of Mouse1's scorable combos). Two further
reasons $\chi^2$ cannot help: the clade search is $10^8$ tests, and $\gamma_{C,z}$ forces
$\sum_{c\in C}(X_c-\tilde p_c)=0$ within each clone — making the clade score exactly minus its
complement's (a genuine within-clone contrast, which is what we want) but imposing a finite-population
correction of roughly $1-m/n_C$ that $V$ ignores, conservatively.

## ⚑ Calibration draws are APPENDED, not held out — the number stops being a launch-time decision

`53` carved `HOLD=3` draws out of $B$ in advance. Unnecessary. Let $M=\max_{b<B}z_b$; compare the
observed to $M$, and later compare **fresh** permutations $B,B+1,\dots$ to the same stored $M$. Swapping
"observed" for "draw $B{+}j$" leaves $M$ untouched, so both are "a value outside the max-set minus the
max of that set" — the identical exchangeability argument. **Each extra draw costs one scan (~0.1% of
the run), so `62 --calib N` can be re-run alone at any time.**

⚠ **NCAL is the ONLY source of precision on the null count** — the $B$ accumulated permutations build
$M$, they are not null draws of the margin, so SE $\approx\mathrm{sd}/\sqrt{\rm NCAL}$. At the budget
operating point (expected false $\approx2$, near-Poisson) HOLD=3 gave $\pm40\%$. **NCAL=12** here.

**Validated on data generated under $H_0$ exactly** (540 cells, 30 tapes, 40 clones of 6–100, independent
Bernoulli, clade structure random and independent of the data, full 3-margin fit, within-clone
permutation, 99,058 scorable combos):
- ties reproduced — $P(m_{\rm obs}\ge0)=0.0117$ observed, 0.0113 calibration, against $1/189=0.0053$ for
  a continuous null;
- quantiles of $m_{\rm obs}$ vs pooled $m_{\rm null}$ agree to 3 decimals through the 99.99th percentile;
- leave-one-out among the held-out draws: mean ratio 1.000, 1.000, 1.001, 1.002, 1.015 → unbiased by
  construction;
- **12 independent replicates**: ratio 0.998 / 1.027 / 1.056 / 1.131 / 1.218 at $t=0,0.5,1,2,3$, $t$-stats
  −0.19 to 1.75, none significant. ⚠ A single replicate had shown 1.55 at $t=2$; it did not reproduce.
  **Bound: up to ~20% optimistic at the smallest counts, not distinguishable from zero at $n=12$.**

## Implementation checks that had to pass

- **`62` reproduces `53` bitwise** — depth-4, observed-only: identical `valid` mask, identical `size`
  vector, max $|\text{diff}|$ on observed $\Lambda_{\rm hard}$ of **exactly 0.000** over 5,574,944 slots.
  Re-checked after the score-test edit: still exact.
- **Depth-6 contains depth-4.** The depth-1–4 partitions are identical between `prefix_codes_{arm}.npz`
  and `prefix_codes6_{arm}.npz` in all five arms, so this run is a superset of the committed one, not a
  substitute. Slot counts at $d\le6$: Mouse3 6,849,492 · Mouse2 9,250,516 · Mouse1 19,059,290 · Pre-TX
  105,877,788 · Subclone 61,727,266 (1.03–1.62× the $d\le4$ counts).
- **Same seed, same 1,000 relabellings as the committed run** — the two are paired permutation by
  permutation, so hard-vs-soft differences are properties of the statistic with Monte Carlo differenced out.
- ⚑ **$V\to0$ cannot blow up $z$.** Since $k\le m$, $z\le m(1-\bar e)/\sqrt V$, so as $\bar e\to1$ the
  numerator vanishes at least as fast as the denominator: measured max $z$ is **0.06** across the 638,355
  slots with $\bar e>0.999$ (11% of all slots — tapes essentially absent from their whole clone). That is
  the $\gamma_{C,z}$ margin working: a clone-wide loss is fitted away, so no sub-clade of it can look
  elevated. ⚠ A float32 check of the Jensen bound $V\le m\bar e(1-\bar e)$ "fails" on 4.5% of slots purely
  from catastrophic cancellation in computing $1-E/m$ when $E\approx m$; $V\le E$ holds on all
  5,541,360 slots and the bound holds in exact arithmetic.

## Mouse 3 results (2026-09-09) — the three detectors on the SAME 1,000 permutations

| detector | combos above threshold | exp. false | median $\hat\pi$ | $\hat\pi\ge0.99$ | clades $\ge50$ | **of those $z<0.5$** |
|---|---|---|---|---|---|---|
| **$z$ (score)** | 19,844 | 9.2 | **0.993** | 52.7% | 8,168 | **0** |
| $\Lambda_{\rm soft}$ | 63,920 | 8.8 | 0.085 | 12.9% | 48,718 | **45,916** |
| $\Lambda_{\rm hard}$ | 13,674 | 8.8 | 0.989 | 49.8% | 4,588 | 0 |

$\Lambda_{\rm soft}$'s calls are **72% artefact**; the score detector admits **zero** of them while
admitting 78% more combos than $\Lambda_{\rm hard}$ at the same budget.

**Per-stratum margin thresholds** (expected false $\le2$ each): 4–5 **1.720** · 6–9 **1.260** ·
10–19 **1.060** · 20–49 **1.080** · 50–199 **0.760** · 200+ none (no combos).
Tie diagnostic: $P(\text{margin}\ge0)$ 0.2103 observed vs 0.1828 calibration, against 0.00100 continuous.

**⭐ 127 events, 3.19% of all missing entries**, and the completeness readout is now **non-circular**:

$$\hat\pi\ \text{median}=\mathbf{1.000}\quad\text{against}\quad \bar e=0.376\ \text{predicted}$$

$\hat\pi\ge0.99$ in 60.6%; $z$ median 5.23, min 2.26; exact LRT median 22.99; $\hat\delta$ median 30.87
with **60% diverged** (complete separation); smallest called clade **4 cells**; every event beats every
one of its 1,000 permutations.

⚑ **Convergence from three directions**: 127 events / 3.19% (score, $d6$) · 128 / 3.09% (per-combo
$\Lambda_{\rm hard}$ margin, $d4$) · 132 / 3.15% ($\Lambda$ floor-2 stratified budget). Membership overlaps
without being identical — 98 of 124 (clone,tape) pairs shared with the floor-2 catalogue (Jaccard 0.64),
82 common to all three.

⚠⚠ **Two numbers that must now replace older ones.**
1. **The partial fraction is 29.9%, not 8.9%.** The hard route cannot see partial losses by construction,
   so its 8.9% was circular. This is the first honest estimate. ⚠ Still quote it with its threshold.
2. **The $\hat\pi$-by-depth gradient is GONE.** Median $\hat\pi$ by clade depth 1→6 is 1.000, 0.992, 1.000,
   1.000, 0.992, 0.963 — no monotone trend. Under $\Lambda_{\rm soft}$ it ran 0.155→1.000. **So the
   "graded losses are largely clade coarseness" gradient was itself an artefact of the detector.** With an
   honest detector events are complete at *every* depth, and the ~30% partial tail sits in the 6–9 cell
   band (median $\hat\pi$ 0.826) — small clades, where one present cell moves $\hat\pi$ a long way.
   ⚠ This supersedes the `MAX_D=6` reading recorded on 2026-09-03 and re-confirmed on 2026-09-07.

## Fig 5 (`66_fig5_plane.py`) — three standalone panels, and ⚠ panel a has a defect

| panel | file | claim |
|---|---|---|
| a | `fig5a_plane_{arm}.png` | the plane: log-density of all 6.8M combos in $(\hat\pi,\text{margin})$ |
| b | `fig5b_margin_{arm}.png` | survival curves, observed vs null: identical below 0.5, then $10^2$–$10^4$ apart; **no permutation exceeds 2.00**, observed runs to 7.14 |
| c | `fig5c_completeness_{arm}.png` | $\hat\pi$ of called events (84 of 127 in the top bin) against $\bar e$ predicted for the same cells |

⚠⚠ **Panel a's first version was misleading; REPLACED 2026-09-09 by a facet per clade-size stratum**,
where the single orange line in each facet **is** the decision rule and every called event sits above it.
The record of what was wrong, because the mistake is easy to remake: the orange line was a *pooled* maximum —
per $\hat\pi$ bin, the highest margin any calibration draw of any combo reached — so it mixes clade sizes,
which is exactly what the per-combo design avoids. **58 of 127 called events sit below it** (21 in 10–19,
15 in 20–49, 9 in 50–199, 7 in 6–9, 6 in 4–5), every one of them clearing its *own* stratum's threshold.
The root cause is structural: **the decision variable (clade size) was on neither axis**, so no curve in
that plane could be the decision boundary. ⇒ faceted, which puts clade size on the facet. The faceted
panel also shows something the pooled one hid: **the threshold falls monotonically with clade size**
(1.72 / 1.26 / 1.06 / 1.08 / 0.76) while the observed cloud extends further right — a bigger clade needs
a smaller margin because its null is near-continuous, a 4-cell clade needs a larger one because its null
is a coarse lattice with a low ceiling.
Two defects already fixed on inspection: a filled contour also drew a meaningless *lower* null boundary
at margin $\approx-7$ (replaced by the upper envelope), and panel b was a per-bin histogram mislabelled
as cumulative (now a true survival curve, which is also literally the FDR's numerator and denominator).

## Status and what is owed

**Settled.** The detector ($z$), the attributor ($\Lambda_{\rm hard}$), the readout ($\hat\pi$, plus
$\hat\delta$/LRT on the shortlist); the appended-calibration design; Mouse 3 end to end.
**Owed.** (1) the four other arms — same code, Pre-TX 24 G / Subclone 16 G, ~17 min per accumulation part
at Mouse3 scale; (2) panel a refaceted; (3) `42` (clone-wide) is still unstratified — the last marginal
number in the project.

---

# ⭐⭐ Exact permutation moments and analytic p-values (2026-09-10) — `67`, `68`

Justin: *the 12 draws are convoluted and hard to justify; why not LRT-style p-values, accept that
small clades have little power, and drop the bespoke machinery?* Substantially right, and pursuing
it turned up an error in the previous design.

## ⚠⚠ CORRECTION — $V$ is the wrong scale, and my verification of it was circular

I asserted $\mathbb{E}[T]=0$ and $\mathrm{Var}[T]=V=\sum_{c\in S}\tilde p(1-\tilde p)$ hold *exactly*,
and "verified" it by Monte Carlo — **simulating under the model null with $\tilde p$ FIXED**, which is
the assumption being tested. Under the **permutation** null with $\tilde p$ **fitted**, measured from
the stored $s_1,s_2$ of the $B=1000$ run:

| clade | combos | perm mean of $z$ | **perm SD of $z$** |
|---|---|---|---|
| 4–5 | 2,314,455 | −0.0020 | **0.667** |
| 6–9 | 1,679,700 | −0.0024 | 0.599 |
| 10–19 | 1,604,130 | −0.0001 | 0.506 |
| 20–49 | 753,060 | 0.0000 | 0.614 |
| 50–199 | 456,885 | 0.0007 | 0.511 |

Mean right, **scale wrong by ~2×**. Event *calling* was unaffected (the margin compared $z_{\rm obs}$
to permuted $z$ on the same scale, so it cancels), but every interpretation of $z$ as "standard
deviations" was wrong. ⚠ **Notation:** the score is now $T$; earlier notes used $S$ for both the clade
and the score.

## The exact moments (closed form, no permutations)

The $\gamma_{C,z}$ fit forces $\sum_{c\in C}(X_{cz}-\tilde p_{cz})=0$ per (clone, tape) — verified,
$\max|\sum_c r_c| = 2.1\times10^{-4}$. So with $r_c=X_{cz}-\tilde p_{cz}$, a clade is a simple random
sample **without replacement** of size $m$ from the clone's residuals. With $I_c$ the membership
indicator, $\mathbb{E}[I_c]=m/n$, $\mathrm{Cov}[I_c,I_{c'}]=-m(n-m)/[n^2(n-1)]$ (negative: fixed
sample size makes membership competitive), and $\sum_{c\neq c'}r_cr_{c'}=-\sum_c r_c^2$ by the
constraint:

$$\mathbb{E}_\pi[T] = m\bar r$$

$$\mathrm{Var}_\pi[T] = m\frac{n-m}{n-1}\sigma^2$$

The factor $(n-m)/(n-1)$ is what $V$ lacked. Sanity: at $m=n$ the clade is the whole population,
$T\equiv0$, variance vanishes; at $m=1$ the factor is 1 and the variance is $\sigma^2$.

⚠ The population is **block-conditional** — a clade is drawn from its clone's cells *that reached that
depth on that anchor*, not from the whole clone. Using the whole clone overstates the SD by 3–6%
(conservative). **Verified against the stored permutations**: median predicted/empirical SD ratio
1.031–1.064 by stratum, overall 1.046, correlation 0.93, empirical perm mean −0.0027 vs a predicted 0.

## ⚠⚠ CORRECTION 2 — Edgeworth was the wrong tool, and is not used

I proposed an Edgeworth skewness correction. It is a *central* approximation: at $z=6$ with skewness 1
the normal tail is $9.9\times10^{-10}$ and the correction term is $3.6\times10^{-8}$ — **37× the
leading term**, i.e. the expansion has diverged exactly where decisions are made. Replaced by the
**hypergeometric**, which is the *exact* permutation distribution in the constant-$\alpha_c$ limit
(Fisher's exact test on clade membership × missingness): discrete, skew-exact, no expansion.

Two p-values per combo; the **conservative (larger)** one is used.

## Multiplicity: one global threshold, and the strata disappear

$p \le \mathcal{E}/N_{\rm test}$ (Bonferroni at level $\mathcal{E}$), so expected false *combos*
$\le\mathcal{E}$; events are unions of combos and the collapse only merges or drops, so expected false
*events* $\le\mathcal{E}$ too. ⇒ **the candidate-vs-event mismatch that forced the absolute budget
simply does not arise.**
⚑ **The six clade-size strata are gone.** $\mathrm{Var}_\pi$ depends explicitly on $m$ and $n$, so $p$
is already conditioned on clade size; a global cut on $p$ is legitimate where a global cut on the
margin was not. The strata were a patch for an unstandardised statistic.

## Testability is structural, not a choice

A combo is testable iff $\sigma^2>0$ **and** $m<n$. Measured on Mouse3: **53.3%** of valid combos.
The other 47% are cases where the clade *is* its whole block-conditional population, so $T\equiv0$ and
no test exists — the clone-level layer was already fitted away by $\gamma$. ⇒ **This is the principled
version of "don't call tiny clades": for many of them there is literally no test**, and for the rest
the p-value decides.

## Mouse 3 results (2026-09-10)

| criterion | combos passing | events | cell × tape | % missing | min clade | $\hat\pi$ median | $\ge0.99$ | partial |
|---|---|---|---|---|---|---|---|---|
| $\mathcal{E}=1$ | 12,093 | **92** | 3,806 | **2.97%** | 6 | **1.000** (exp. 0.402) | 66.3% | 23.9% |
| $\mathcal{E}=12$ | 15,492 | 128 | 4,852 | 3.79% | 5 | 1.000 (exp. 0.392) | 64.1% | 28.1% |

⭐ **Stringency-matched, the three routes converge**: 128 events (exact, $\mathcal{E}=12$) · 127
(per-combo margin) · 132 (committed $\Lambda$ floor-2 budget). Pairwise overlap 96/154 (clone,tape)
pairs exact-vs-margin, 89 exact-vs-committed, **85 common to all three**. Three different statistics,
the same population.

**Validation against the stored 1,000 permutations** (124,490 combos where they resolve):
median $p_{\rm norm}/p_{\rm emp}=2.31$, median $p_{\rm hyper}/p_{\rm emp}=13.0$ — **both conservative**.
⚠ The hypergeometric's 13× conservatism is a real power cost: at $\mathcal{E}=1$ the normal alone
clears 26,973 combos, the hypergeometric 13,009, the max 12,093. Taking the max is a deliberate
conservative choice, not a free one.

## What this improves

| | margin route (`63`) | exact route (`67`) |
|---|---|---|
| scale | model $V$, wrong by ~2× | exact permutation variance, verified to 3–6% |
| null reference | max over 1,000 permutations | closed-form moments + normal/hypergeometric |
| calibration draws | 12 appended | none |
| multiplicity | 6 hand-drawn strata × budget 2 | one global threshold |
| ties | 18–32% of combos tie the max | no max taken, so none |
| censoring | $p$ floored at $1/1001$ | continuous |
| error control | expected false *candidates* | expected false **events** $\le\mathcal{E}$ |
| permutations needed | 1,012 per arm | **0** (validation only) |
| runtime | 10 min merge | **9 s**, 1.76 GB |

⚠ **What is given up.** The normal tail is validated only to $p\approx10^{-3}$ (1,000 permutations
cannot resolve further) while events sit near $10^{-9}$; beyond that it is extrapolation, mitigated
by the hypergeometric which does not degrade in the tail. Bonferroni over ~$10^6$ *dependent* tests
(each loss detected ~160×) over-corrects, costing power. And this is the **third** detector design in
three days — the first two looked right until measured, so it is held to the same standard: it
reproduces the previous catalogues at matched stringency, and its two approximations are checked
against the stored permutations above.

## Five arms, exact route (2026-09-10)

| arm | testable | events | (clone,tape) pairs | ev/pair | **cell × tape** | **% of missing** | $\hat\pi$ med | min clade | $p_{\rm norm}/p_{\rm emp}$ |
|---|---|---|---|---|---|---|---|---|---|
| Subclone | 99.9% | 5,863 | 1,382 | **4.2** | 355,849 | **25.48%** | 0.969 | 4 | **0.46** |
| Pre-TX | 76.3% | 821 | 817 | 1.0 | 17,477 | 1.33% | 1.000 | 6 | 1.32 |
| Mouse1 | 76.3% | 343 | 317 | 1.1 | 27,366 | 6.74% | 1.000 | 5 | 1.66 |
| Mouse2 | 63.4% | 326 | 163 | 2.0 | 25,493 | 7.70% | 1.000 | 6 | 1.85 |
| Mouse3 | 53.3% | 92 | 91 | 1.0 | 3,806 | 2.97% | 1.000 | 6 | 2.31 |

Runtime 8 s – 12 min per arm on stored parts. Against the committed $\Lambda$-budget catalogue
(20.07 / 2.46 / 4.34 / 7.02 / 3.15%) the ordering is the same apart from Mouse3 and Pre-TX swapping.
**(clone, tape) membership agrees strongly**: shared with the $\Lambda$-budget catalogue 1,240/1,382
(Subclone), 774/817 (Pre-TX), 154/163 (Mouse2), 217/317 (Mouse1), 72/91 (Mouse3).

⚠⚠ **Subclone's normal p-value is ANTI-conservative** ($p_{\rm norm}/p_{\rm emp}=0.46$; every other arm
is 1.32–2.31) and **I have no clean explanation.** The obvious candidate — that its clades are a tiny
fraction of their clones (median 0.006) — fails, because Mouse2 is nearly the same (0.014) and comes
out conservative at 1.85. ⇒ **the $\max(p_{\rm norm},p_{\rm hyper})$ guard is load-bearing, not
belt-and-braces**: the hypergeometric is the binding (larger) p in **89.4%** of Subclone's called
events. I had suggested dropping it for power; that suggestion is **withdrawn**.

⚠⚠ **Subclone fragments: 4.2 events per (clone, tape), up to 26 on one**, median called clade 0.6% of
its clone. Most likely the larger clade carrying the loss is only partially missing (Subclone's
$\hat\pi$ median 0.969 is the lowest of the five) so it fails while its complete sub-pieces pass.
Pre-TX and Mouse3 sit at 1.0. ⇒ **quote cell × tape, never event counts** — the standing guidance,
and this is its sharpest demonstration. Cell × tape is safe: the collapse guarantees kept events for
one (clone, tape) share no cells, so nothing is double-counted.
