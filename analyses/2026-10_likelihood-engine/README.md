# likelihood engine — findings

## Session 26 (2026-10-12/13) — sampler discussion; HMC deep dive (Betancourt §1–§3.1); teaching figures

No scientific run. HMC deep dive under way (Betancourt 2017, `refs/hmc.pdf`, his notation); records in
`notes/sciphy_notes.md` §G.4 (sampler roles; HMC-within-Gibbs §G.4.9; reporting §G.4.10; Betancourt
§1–§3.1 §G.4.11–G.4.16). **Scaling probe** (`src/11_probe_rwm_hmc_scaling.py`, Gaussian targets, seconds):
random-walk Metropolis ≈2.5×D evaluations per independent draw vs HMC 3–18 gradients for D = 10–1,000
(≈117× fewer at D=400, ≈50× wall-clock at our gradient/likelihood cost ratio); along an SD ratio k, RWM ∝ k²,
HMC ∝ k. *Reading rule: evaluations per independent draw — fewer is better; the ratio is the speed-up on an
ideal target, not against BEAST's one-node operators.*
**Teaching figures** (Justin's request; `src/07_fig_coordinates.py`, seconds, login node):
`figures/fig_coord_A_exp1.png` — $\lambda\sim$Exp(1) summarised in $\lambda$ and $\log\lambda$;
`figures/fig_coord_B_node.png` — an **illustrative** late node, $h\sim$Beta(1.5, 6), not fitted,
summarised in calendar / ratio-coordinate / recorder time. *Reading rule (both): where the coloured
marks coincide the summary is a property of the posterior; where they split, of the coordinates.*
Mode, HPD and "average the coordinate then convert" split; median, equal-tailed interval and the mean
of the reported quantity coincide. Fig B: MAP date 36.0 / 40.9 / 42.6 across the three rulers vs a
ruler-free median of 37.2. Grid values cross-checked by 2×10⁶-draw Monte Carlo. ⚠ The palette
validator could not run (cluster Node v10 too old); colours are the reference palette's slots 1–3,
documented as passing all-pairs.
**Teaching figures for Betancourt §1.2** (`src/08_fig_density_volume.py`): `figures/fig_volume_1_rings2d.png`
— 2-D Gaussian in rings of width 0.5 SD: density per ring falls, ring area grows, probability per ring
peaks at 0.5–1.5 SD (28% / 28%) while the innermost ring holds 12%; 2×10⁶ draws agree to 3 decimals.
`figures/fig_volume_2_shellsD.png` — density × shell size in $D$ = 1/2/10/100; peak at $\sqrt{D-1}$ =
0/1/3/9.95 SD. *Reading rule: the shaded curve is where draws land; its peak moves away from the mode as
$D$ grows.*
`figures/fig_volume_3_allcoords.png` (`src/09_fig_all_coordinates.py`) — "near the mode means near the centre
in every coordinate at once": $P$(all within ±$c$ SD) $=p_c^D$ against $D$ (below 1% at $D$=12 / 99 / ~1,700
for $c$=1/2/3; our $D$≈200–800 marked), one $D$=400 draw coordinate by coordinate, and its distance
accumulating as $\sqrt k$ to 19.5 SD. *Reading rule: in panel a, a curve near the bottom means
"essentially never".*
`figures/fig_funnel_pathology.png` (`src/10_fig_funnel.py`) — Betancourt §2.3 on Neal's funnel: random-walk
Metropolis at 37% acceptance, 4 chains from exact draws: $P(v<-3)$ estimates 0.05–0.34 against truth 0.159,
split $\hat R$ 1.08 / 1.16; non-centred, same sampler and step: 0.14–0.16, $\hat R$ 1.00. *Reading rule: a healthy
run converges onto the dashed truth with all chains together.*

**Validation figures for the PI update (2026-10-13; Justin approved fig 1 of a proposed set of four).**
`figures/fig_val1_sciphy.png` (`src/12_fig_val_sciphy.py`, replots `results/validate_sciphy.json`, nothing
recomputed). Revised to one panel at Justin's request: total log-likelihood per dataset, SciPhy vs ours, 50
datasets, V6a coding only (blanks as observed blanks), on the identity line; annotated with the worst per-tape
relative difference, 1.4e-15 over 1,500 tapes. *Reading rule: points on the line = same number.* Framing:
SciPhy's likelihood = ours at π = 0, Γ_T = 0, no censoring, λ₀ constant in time, one ξ, blank founder (per-tape
r_z kept — SciPhy has per-tape clocks). The simulated data contain no dropout or silencing.
`figures/fig_val2_exact.png` (`src/13_fig_val_exact.py`, replots `validate_engine.json` + `validate_closeout.json`):
one row per exact check — revised at Justin's request to V1 enumeration, V1 full state space, V2, V3 only (V4, V8a dropped from the figure, still PASS), mark = worst error over the
check (V1/V4 saved only their worst, so all rows show the worst), against each row's pre-fixed tolerance and the
float64 rounding band. Worst 1.3e-15 to 1.8e-14; margin to tolerance 90× (V4) to 7.5×10⁶× (V8a). V1/V2 at random
parameters with every layer on; V3 dropout + censoring in 90/90 clones, silencing 67, time-varying rate/comp 22.

`figures/fig_val4_gradients.png` (`src/15_fig_val_gradients.py`; approved "you can make this figure"): (a) all 360
V5 relative errors by quantity vs the 1e-5 line — 359 at ≤2.5e-7, the floor ~1e-8 set by the h = 1e-5 reference itself;
(b) the outlier's step ladder, **re-run** by `src/14_probe_grad_ladder.py` (job 19801344, seconds) → `results/v5_ladder.json`.
The probe reproduces the recorded autodiff (2679.07139773103) and h = 1e-5 value (2679.7352913490613) exactly; central-difference
error vs autodiff 2.5e-2 / 2.5e-4 / 2.5e-6 / 2.3e-8 / 5.2e-8 at h = 1e-4 … 1e-8 (h = 1e-3 undefined: evaluates Γ_T < 0);
the rise at 1e-8 is roundoff. ⚠ Minor: Session 22's "2.6e-2 at h 1e-4" divided by |autodiff|; V5's definition divides by
max(|fd|, 1), giving 2.5e-2. Same diagnosis, no conclusion moves. V5′ (formal ladder over all entries) still open.

## Session 25 (2026-10-12) — validation closeout: V8 and V6 PASS; V5′ re-proposed before running

Approved by Justin 2026-10-12 ("go ahead with the engine validation"). Proposal as approved: V5′
(gradients against a Richardson-extrapolated finite difference at a fixed step $h=10^{-4}\max(|\theta|,1)$,
one-sided at a boundary, tolerance $10^{-5}$ unchanged), V8a (new hard check), V8⁺ (2,000 fresh trees at
the failing cell + 4 neighbours), V6 (SciPhy's Java likelihood, base model). `03_validate_engine.py` is
untouched as the record of the 2026-10-06 run. New: `src/05_validate_closeout.py` (V8a, V8⁺),
`src/06_validate_sciphy.py` + `src/06_sciphy_harness/V6Harness.java` (V6). Results
`results/validate_closeout.json`, `results/validate_sciphy.json`.

**Reading rule (pre-fixed).** Hard checks: worst error below tolerance. Statistical: $z$ = mean/se across
trees, no effect = 0; $|z|\le3$ PASS, 3–5 marginal, $>5$ defect (a real 20-tree deficit would return near
$z\approx-30$ at 2,000 trees).

| check | compared against | configs | worst | tol | verdict |
|---|---|---|---|---|---|
| V8 regression | the V8 statistic recomputed on the library trees of the failing cell vs the recorded $z$ | 10 deciles | 0.0 | 1e-9 | reproduces −3.00 / −3.03 exactly |
| V8a | engine's prior density on one branching height, integrated $0\to x$, vs the simulator's CDF `bds_G` | 360 (18 cells × 20 $x$) | 1.3e-15 | 1e-8 | **PASS** |
| V8⁺ | fresh trees (2,000/cell, seed 61012026, disjoint from the library): decile occupancy + root PIT, $|z|$ | 55 (5 cells × 11) | \|z\| 2.55 | 3.0 | **PASS** |
| V6a | engine vs SciPhy's `SciPhyTreeLikelihood`, per-tape ℓ, blank tape coded as an observed blank | 1,500 (50 trees × 30 tapes) | 1.4e-15 | 1e-8 | **PASS** |
| V6b | the same, blank tape coded MISSING (Park's coding; exercises the all-missing-subtree path) | 1,500 | 1.4e-15 | 1e-8 | **PASS** |
| V6 teeth | SciPhy given every clock rate × (1 + 10⁻⁴) | 30 | 1.6e-5 | must exceed 1e-8 | detected, 1,600× tol |

**V8⁺ in words.** At the cell that tripped (64 cells, ρ = 0.25, turnover δ/b = 0.3) the two low deciles
are **+0.53 and +0.35** at 2,000 trees, against −3.00 / −3.03 at 20. Worst |z| over the 55 statistics is
2.55 (n = 25, decile 9), about the expected maximum of 55 chance draws; every decile holds 0.1 ± 0.0035
of the branchings. ⇒ the 2026-10-06 trip was chance.
⚠⚠ **And the 2026-10-06 rule was miscalibrated — a defect in my check design, third time in this
project.** With 20 trees per cell, $z$ is $t_{19}$, not normal: $P(|t_{19}|>3)=0.0074$, so **1.46 trips
were expected by chance over 198 cells, not the 0.53 the rule assumed**. The observed SD of $z$ over all
198 cells, 1.058, matches $t_{19}$'s 1.057. Same lesson as sessions 18 (check E) and 21 (rung-2 smoke):
*a t-statistic over ~20 trees needs a t threshold, or far more trees.*

**V6 in words.** Our engine reproduces SciPhy's own code (azwaans/SciPhy commit `2a49194`, source
unmodified, compiled against BEAST 2.7.7; one `SciPhyTreeLikelihood` per tape) to $1.4\times10^{-15}$ —
float rounding — with dropout, censoring and silencing switched off. Coverage: library trees $n$ = 25/64/210
× 3 ρ × 2 turnovers, $\Lambda_T$ ∈ {1, 2.5, 4, 5.5, 8} edits/tape, $k$=30, $M$=114, Park $\xi^L$; full tapes
0–79% of entries (the saturation tail), blank tapes 0–47% (in V6b these all enter as missing). ⇒ **our
likelihood contains SciPhy's as an exact special case.** This is the only check written by someone else;
it rules out a misconception shared by engine and simulator *for the base model*. The emission and
silencing layers still rest on first-principles enumeration (V1/V2). Infrastructure (not committed):
BEAST 2.7.7 and the compiled classes under `/data1/choij10/justin/tools/`.

**⚠ V5′ was NOT run — the approved step size was wrong, found by a probe before writing it.** The
2026-10-06 diagnosis implies the failing entry ($\partial J/\partial\Gamma_T$ at $\Gamma_T=0$, tree 16,
autodiff 2,679.0714) behaves like $\log(a+b\Gamma_T)$ with $b/a\approx2.7\times10^3$, i.e. it curves on
a scale $\approx3.7\times10^{-4}$ in $\Gamma_T$ — comparable to the proposed $h=10^{-4}$. Probe (one-sided
Richardson, relative error vs autodiff): **2.6e-3 at h = 1e-4** (would FAIL), 4.7e-6 at 1e-5, 4.6e-9 at
1e-6, 3.9e-11 at 1e-7 — converging on autodiff at the $O(h^3)$ rate. ⇒ the engine gradient is correct
(again); my proposed reference would have manufactured a failure. A fixed smaller step is not safe either:
roundoff grows as $\epsilon|J|/h$ and $|J|$ reaches $10^4$ nats on the larger trees. Re-proposed to Justin
as a step ladder selected by agreement between adjacent steps (never by agreement with autodiff).

## Session 22 (2026-10-06) — engine built; validation FAILED BY RULE on two check-set issues

Approved by Justin 2026-10-06 (CLAUDE.md, PROPOSAL parts 1–2, "proceed with building using all proposed
decisions"). Write-up: `writeup/likelihood_engine.tex` (10 pp). Code: `src/01_engine_ref.py` (NumPy
reference twin + full-state-space and literal-enumeration references), `src/02_engine.py` (JAX
engine, float64), `src/03_validate_engine.py`, `src/04_time_engine.py`. Environment
`/data1/choij10/justin/envs/pando-jax` (python 3.12, jax 0.11.2 with CUDA 12, numpy 2.5.3, scipy
1.18.1), built **from conda-forge only** — the first attempt was pulling from Anaconda's commercial
`repo.anaconda.com` channel (inherited from the user's conda config) and was stopped. `scripts/submit.sh`
gained optional `--py` and `--gres` flags (default behaviour unchanged).

**Full validation (seed 6102026, job 17612197, 6 min, 2.5 GB). Pre-fixed rule: hard checks below
tolerance; statistical worst |z| ≤ 3. Verdict: FAIL — V5 over tolerance; V8 two trips.**

| check | compared against | configs | worst | tol | verdict |
|---|---|---|---|---|---|
| V1 | twin vs dense full state space (41 states, no slot logic) | 340 | 3.7e-15 | 1e-10 | PASS |
| V1 | twin vs literal enumeration of every hidden history | 144 | 5.8e-15 | 1e-10 | PASS |
| V2 | sum of L over every possible set of reads = 1 | 30 | 1.8e-14 | 1e-10 | PASS |
| V3 | JAX vs twin, per-tape ℓ, simulator data on library trees n 25–210 | 2,700 | 1.9e-15 | 1e-10 | PASS |
| V4 | closed-form kernel vs matrix exponential | 1,000 | 1.1e-14 | 1e-12 | PASS |
| V5 | autodiff vs central finite differences (step 1e-5) | 360 | **2.5e-4** | 1e-5 | **FAIL** |
| V7 | mean score at the truth (12 cells: 2 settings × 6 scores) | 12 | \|z\| 1.96 | 3.0 | PASS |
| V8 | tree prior vs library branching heights (18 cells × 11) | 198 | \|z\| **3.03** | 3.0 | **2 trips** |

**⚠ V5 — a defect in my CHECK, not the engine (diagnosed, job 17615735).** One comparison of 360
failed: d(Σℓ)/dΓ_T on a tree with the "no silencing" setting, i.e. at Γ_T = 0, where the gradient is
large (2,679.07) and strongly curved. All other 359 agree to ≤ 2.5e-7. Shrinking the finite-difference
step: central-difference error 2.6e-2 (h 1e-4) → 2.5e-4 (1e-5) → 2.5e-6 (1e-6) → 2.3e-8 (1e-7) —
falling 100× per 10× step, the signature of step-size (truncation) error in the reference; it converges
to the autodiff value. (At h 1e-3 the difference is undefined: it evaluates Γ_T < 0.) ⇒ the gradient
is correct; the check's fixed step was too coarse at a boundary point of high curvature. **Not patched —
an amended V5 is proposed to Justin.**

**⚠ V8 — two adjacent deciles in ONE library cell (n = 64, ρ = 0.25, θ = 0.3): z −3.00 and −3.03**
— fewer branchings in the lowest two deciles of the prior's CDF (i.e. near harvest) than the iid
conditioned density predicts. Same cell, same numbers as the quick run, because V8 reads the FIXED
library trees (20 per cell): the "re-run at a new seed" rule cannot apply. No global bias (198 cells:
mean z ≈ −0.1, SD ≈ 1.1). Within a cell the deciles share trees and sum to zero, so the two trips are
closer to one deviation than two. Unresolved: chance in 20 trees, or a real mismatch between the prior
formula and the library's conditioning. **Proposed follow-up: fresh trees for that cell.**

**Timing (descriptive, k = 30 tapes; median of 5; float64).**

| cells n | likelihood, CPU 1 core | gradient, CPU | likelihood, A100 | gradient, A100 |
|---|---|---|---|---|
| 25 | 0.006 s | 0.008 s | 0.025 s | 0.034 s |
| 210 | 0.039 s | 0.094 s | 0.19 s | 0.25 s |
| 1,976 | 0.33 s | 2.1 s | 1.7 s | 2.2 s |
| 11,081 | 2.0 s | **121 s** | 9.5 s | 12.9 s |

Compile 5–8 s (130 s CPU / 38 s GPU at n = 11,081, gradient included). ⚠ **The forward likelihood
scales linearly (CPU), but the CPU gradient is super-linear at large n** (2.1 s → 121 s for 5.6× more
cells, ≈ n^2.3): reverse-mode through the per-node `lax.scan` saves the whole per-node array at every
step. The GPU is slower than one CPU core on the forward pass at every size (thousands of small
sequential steps). Both point to the same fix — process nodes **level by level** (sequential steps =
tree depth, not node count) — an implementation change that would need V3/V5 re-run. Not done.

**Not yet done:** V6 (SciPhy's Java likelihood cross-check; separate setup: BEAST 2.7 + SciPhy package
+ the `java` module).
