#!/usr/bin/env python3
r"""
What exactly is rho a correlation between?  -- checked by brute force.

rho_cell is the correlation between TWO 0/1 ENTRIES OF THE MATRIX drawn from the
same cell.  Build the dataset literally: draw a cell c and two distinct tapes
t, t'; record U = Y_ct and V = Y_ct'.  Millions of such rows, two columns of 0s
and 1s, Pearson correlation of the two columns.  rho_tape is the same with a tape
and two distinct cells.

We never compute it that way in practice -- the variance identity gives the same
number in O(nk) instead of O(nk^2):

    Var(R_c) = Var(sum_t Y_ct) = sum_t Var(Y_ct) + sum_{t != t'} Cov(Y_ct, Y_ct')
             = k p(1-p) + k(k-1) sigma^2
    => VIF = Var(R_c)/[k p(1-p)] = 1 + (k-1) rho

so rho is literally the average covariance over all within-unit pairs, rescaled.
This script confirms the two routes agree.

Also checked here:
  * Spearman == Pearson for two binary variables (ranking a 0/1 variable is an
    affine map, and Pearson is invariant to affine maps), so the question
    "Pearson or Spearman?" has no content for these data.
  * the ANOVA reading: rho = between-unit variance / total variance, with the
    naive sample version and the noise-corrected population version reported
    separately -- they differ for cells (k = 166 items) and not for tapes.
  * the conditional-probability reading P(miss | miss) = (1-p) + rho p.
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ARM, NPAIR, SEED = "Mouse1", 3_000_000, 20260901

z = np.load(RES / f"dropout_matrix_{ARM}.npz", allow_pickle=False)
Y = z["recovered"]
n, k = Y.shape
p = float(Y.mean())
rng = np.random.default_rng(SEED)
out = {"arm": ARM, "n_cells": n, "n_tapes": k, "p": p, "n_pairs_sampled": NPAIR}
print(f"{ARM}: {n:,} cells x {k} tapes, {n*k:,} entries, p = {p:.4f}\n")

for name, m in (("cell", k), ("tape", n)):
    # --- route 1: the variance identity
    R = (Y.sum(1) if name == "cell" else Y.sum(0)).astype(float)
    vif = R.var(ddof=1) / (m * p * (1 - p))
    rho_var = (vif - 1) / (m - 1)

    # --- route 2: brute force over pairs of entries sharing a unit
    if name == "cell":
        u = rng.integers(0, n, NPAIR)                       # the cell
        a = rng.integers(0, k, NPAIR); b = rng.integers(0, k - 1, NPAIR)
        b += (b >= a)                                       # b != a
        U, V = Y[u, a], Y[u, b]
    else:
        u = rng.integers(0, k, NPAIR)                       # the tape
        a = rng.integers(0, n, NPAIR); b = rng.integers(0, n - 1, NPAIR)
        b += (b >= a)
        U, V = Y[a, u], Y[b, u]
    U = U.astype(np.float64); V = V.astype(np.float64)
    r_pear = float(pearsonr(U, V).statistic)
    r_spear = float(spearmanr(U, V).statistic)
    se = 1 / np.sqrt(NPAIR - 3)                             # rough Fisher-z se

    # --- route 3: the ANOVA reading
    means = Y.mean(1) if name == "cell" else Y.mean(0)
    between_naive = float(means.var(ddof=1))
    total = p * (1 - p)
    noise = total / m                                       # sampling var of a unit mean
    between_corr = max(between_naive - noise, 0.0)

    # --- the conditional-probability reading
    miss = U == 0
    pmm_emp = float((V[miss] == 0).mean())
    pmm_pred = (1 - p) + rho_var * p

    out[name] = {"items_summed": m, "vif": float(vif), "rho_variance_identity": float(rho_var),
                 "rho_pearson_pairs": r_pear, "rho_spearman_pairs": r_spear,
                 "pair_se": float(se), "sigma": float(np.sqrt(rho_var * total)),
                 "anova_total_var": total, "anova_between_naive": between_naive,
                 "anova_between_noise_corrected": between_corr,
                 "R2_naive": between_naive / total, "R2_corrected": between_corr / total,
                 "P_miss": 1 - p, "P_miss_given_miss_empirical": pmm_emp,
                 "P_miss_given_miss_predicted": float(pmm_pred)}
    print(f"--- rho_{name}   (two entries sharing a {name}; {m:,} items summed per unit)")
    print(f"  variance identity   (VIF-1)/(m-1)      {rho_var:.5f}")
    print(f"  Pearson  over {NPAIR:,} sampled pairs  {r_pear:.5f}   (se ~{se:.5f})")
    print(f"  Spearman over the same pairs           {r_spear:.5f}   <- identical: binary")
    print(f"  ANOVA  between-unit var / total var    {between_naive/total:.5f} naive, "
          f"{between_corr/total:.5f} noise-corrected")
    print(f"     total var p(1-p) = {total:.5f};  between-unit var = {between_naive:.5f} "
          f"(naive) / {between_corr:.5f} (corrected);  sampling noise = {noise:.6f}")
    print(f"  P(miss) = {1-p:.4f}  ->  P(miss | miss in same {name}) "
          f"= {pmm_emp:.4f} measured, {pmm_pred:.4f} predicted by (1-p)+rho*p\n")

(RES / "rho_check.json").write_text(json.dumps(out, indent=1))
print("wrote results/rho_check.json")
