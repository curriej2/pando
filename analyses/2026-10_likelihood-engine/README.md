# likelihood engine — findings

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
