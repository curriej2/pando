#!/usr/bin/env python3
r"""
Is the pooled rho one clone in disguise?  -- per-clone decomposition.

rho_hat = SUM_C num_C / SUM_C den_C is a RATIO OF SUMS, so each clone enters
weighted by its pair count n_C(n_C-1).  Measured shares of the pooled weight held
by the single largest clone: Mouse2 98.9%, Mouse1 84.0%, Subclone 47.6%,
Mouse3 28.8%, Initial 1.6%.  So for two arms the pooled number is essentially one
clone, and "5/5 arms replicate" overstates the independence of the replication.

Three things computed here, all against the same within-sample permutation null:

  1. per-CLONE excess, pooling tapes inside each clone:
         rho_C = SUM_z [ (SUM_c r)^2 - SUM_c r^2 ] / SUM_z [ (SUM_c sqrt v)^2 - SUM_c v ]
     -> answers "which clones contribute", which the pooled statistic cannot.
  2. pooled excess with the LARGEST clone dropped -- does the signal survive?
  3. equal-weight-per-clone excess (mean over clones of rho_C), which weights a
     4-cell clone the same as a 1,607-cell one.

Output: results/pool_robustness_{arm}.json
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
JUNK, MIN_CLONE = 2, 5

arm = sys.argv[1]
B = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 500
rng = np.random.default_rng(20260902)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, dep, trm, clone, sample = (z0["recovered"], z0["depth"].astype(np.int16),
                              z0["term"], z0["clone"], z0["sample"])
s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, dep, trm, clone, sample = Y[s], dep[s], trm[s], clone[s], sample[s]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
_, g_samp = np.unique(sample, return_inverse=True)
Gc, Gs = g_clone.max() + 1, g_samp.max() + 1
sizes = np.bincount(g_clone).astype(float)
big = int(np.argmax(sizes))
print(f"{arm}: {n:,} cells, {Gc:,} clones, largest {int(sizes[big]):,} cells", flush=True)


def fit_twoway(M, W, iters=40):
    cm = np.clip((M * W).sum(0) / np.maximum(W.sum(0), 1), 1e-3, 1 - 1e-3)
    beta = np.log(cm / (1 - cm)); alpha = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))
            if ax == 0:
                g = (W * (M - p)).sum(1); h = (W * p * (1 - p)).sum(1)
                alpha = np.clip(alpha + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
                alpha -= alpha.mean()
            else:
                g = (W * (M - p)).sum(0); h = (W * p * (1 - p)).sum(0)
                beta = np.clip(beta + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :]))), 1e-6, 1 - 1e-6)


def per_clone(r, v, sv):
    """num_C and den_C, pooling tapes within each clone."""
    num = np.zeros(Gc); den = np.zeros(Gc)
    for k in range(K):
        S1 = np.bincount(g_clone, weights=r[:, k], minlength=Gc)
        S2 = np.bincount(g_clone, weights=r[:, k] ** 2, minlength=Gc)
        T1 = np.bincount(g_clone, weights=sv[:, k], minlength=Gc)
        T2 = np.bincount(g_clone, weights=v[:, k], minlength=Gc)
        num += S1 * S1 - S2
        den += T1 * T1 - T2
    return num, den


# the full T0 sweep, so the calibration curve can be rebuilt on the
# equal-clone-weighted estimator rather than the one-clone-dominated pooled one
families = {"missing": ((~Y).astype(float), np.ones_like(Y, dtype=float))}
for L in range(2, 7):
    families[f"depth>={L}"] = ((dep >= L).astype(float),
                               (Y & ((dep >= L) | (trm != JUNK))).astype(float))
samp_pos = [np.flatnonzero(g_samp == s_) for s_ in range(Gs)]
out = {"arm": arm, "n_cells": int(n), "n_clones": int(Gc), "n_perm": B,
       "largest_clone_cells": int(sizes[big]),
       "largest_clone_weight_share": float((sizes[big] * (sizes[big] - 1)) /
                                           (sizes * (sizes - 1)).sum()), "families": {}}

for name, (M, W) in families.items():
    P = fit_twoway(M, W)
    r = np.where(W > 0, M - P, 0.0); v = np.where(W > 0, P * (1 - P), 0.0)
    sv = np.sqrt(v)
    num, den = per_clone(r, v, sv)
    keep = (sizes >= MIN_CLONE) & (den > 0)
    nulls = np.empty((B, Gc)); nulls_den = np.empty((B, Gc))
    for b in range(B):
        perm = np.arange(n)
        for pos in samp_pos:
            perm[pos] = rng.permutation(pos)
        nb, db = per_clone(r[perm], v[perm], sv[perm])
        nulls[b], nulls_den[b] = nb, db
    with np.errstate(invalid="ignore", divide="ignore"):
        rho_C = np.where(den > 0, num / den, np.nan)
        rho_C_null = np.where(nulls_den > 0, nulls / nulls_den, np.nan)
    exc_C = rho_C - np.nanmean(rho_C_null, 0)
    # pooled, all clones and with the dominant clone removed
    msk = np.ones(Gc, bool); msk[big] = False
    pooled = num.sum() / den.sum()
    pooled_null = (nulls.sum(1) / nulls_den.sum(1)).mean()
    drop = num[msk].sum() / den[msk].sum()
    drop_null = (nulls[:, msk].sum(1) / nulls_den[:, msk].sum(1)).mean()
    eq = float(np.nanmean(exc_C[keep]))
    out["families"][name] = {
        "marginal": float((M * W).sum() / max(W.sum(), 1)),
        "pooled_excess": float(pooled - pooled_null),
        "pooled_excess_drop_largest": float(drop - drop_null),
        "equal_weight_excess": eq,
        "n_clones_tested": int(keep.sum()),
        "frac_clones_positive": float(np.nanmean(exc_C[keep] > 0)),
        "clone_excess_quantiles": {str(q): float(np.nanquantile(exc_C[keep], q))
                                   for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        "clone_sizes": sizes[keep].tolist(),
        "clone_excess": np.nan_to_num(exc_C[keep]).tolist()}
    d = out["families"][name]
    print(f"  {name:9s} pooled {d['pooled_excess']:+.4f} | "
          f"drop largest {d['pooled_excess_drop_largest']:+.4f} | "
          f"equal-weight {d['equal_weight_excess']:+.4f} | "
          f"{d['n_clones_tested']} clones >= {MIN_CLONE} cells, "
          f"{100*d['frac_clones_positive']:.0f}% positive | "
          f"median clone excess {d['clone_excess_quantiles']['0.5']:+.4f}", flush=True)

(RES / f"pool_robustness_{arm}.json").write_text(json.dumps(out, indent=1))
print(f"wrote results/pool_robustness_{arm}.json")
