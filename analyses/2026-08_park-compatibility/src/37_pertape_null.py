#!/usr/bin/env python3
r"""
================================================================================
 Per-tape permutation nulls -> a p-value and a z-score for every tape
================================================================================

34_lineage_rho.py permutes cells within harvest sample and keeps only the POOLED
rho.  But rho is computed per tape, and the per-tape values are wildly
heterogeneous -- in Mouse1 the median tape sits at +0.049 while the top decile
carries 59% of the excess.  A single pooled number hides exactly the structure a
Dollo mechanism predicts: discrete losses at PARTICULAR loci, not a uniform shift.

So: keep the whole per-feature vector on every permutation draw.  For tape z,

    rho_obs(z)                          on the true clone labels
    rho_perm(z, b),  b = 1..B           cells reshuffled within sample, clone
                                        sizes and sample composition preserved

    z-score  = [rho_obs - mean_b] / sd_b
    p        = (1 + #{b : rho_perm >= rho_obs}) / (B + 1)     one-sided, upper
               (the +1 makes it a valid permutation p-value; it is bounded below
                by 1/(B+1), which is why B must exceed ~m/alpha for m tapes)

One-sided upper because the hypothesis is EXCESS concordance; nothing predicts
sisters agreeing LESS than strangers, so the count of significantly negative
tapes is a null-calibration check -- it should be near zero, and if it is not the
null is wrong.

BH-FDR across the 166 tapes within each arm, since we expect a subset to be real
and want to count them rather than test one pre-chosen tape.

⚠ POWER IS NOT UNIFORM ACROSS TAPES.  A tape recovered in 0.6% of cells has
almost no variance to detect (Var <= p(1-p)), so it will not reach significance
however heritable it is.  The denominator SUM_C [(SUM sqrt(v))^2 - SUM v] is the
natural "effective information" per tape and is reported alongside, so the count
of significant tapes can be read against the count that was ever detectable.

A T0 family (depth>=6, known heritable, marginal matched to missingness) is run
identically as a positive control: it says how many tapes light up when the
feature definitely IS heritable.

Usage: 37_pertape_null.py <arm> [--screen] [--nperm B] [--fdr 0.05]
Output: results/pertape_{arm}{_screen}.json  (+ .npz of the raw null draws)
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
JUNK = 2

arm = sys.argv[1]
SCREEN = "--screen" in sys.argv
B = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 5000
QFDR = float(sys.argv[sys.argv.index("--fdr") + 1]) if "--fdr" in sys.argv else 0.05
rng = np.random.default_rng(20260902)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, dep, trm, clone, sample = (z0["recovered"], z0["depth"].astype(np.int16),
                              z0["term"], z0["clone"], z0["sample"])
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, dep, trm, clone, sample = Y[sel], dep[sel], trm[sel], clone[sel], sample[sel]
if SCREEN:
    good = ~np.load(RES / f"collision_flags_{arm}.npz", allow_pickle=False)["flag"]
    Y, dep, trm, clone, sample = Y[good], dep[good], trm[good], clone[good], sample[good]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
_, g_samp = np.unique(sample, return_inverse=True)
Gc, Gs = g_clone.max() + 1, g_samp.max() + 1
print(f"{arm}{' +screen' if SCREEN else ''}: {n:,} cells, {K} tapes, {Gc:,} clones, "
      f"{Gs} samples, B={B:,}", flush=True)


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


def pair_sums(x2, x1, lab, G):
    out = np.empty(x1.shape[1])
    for k in range(x1.shape[1]):
        S1 = np.bincount(lab, weights=x1[:, k], minlength=G)
        S2 = np.bincount(lab, weights=x2[:, k], minlength=G)
        out[k] = (S1 * S1 - S2).sum()
    return out


def bh(p, q):
    m = p.size
    o = np.argsort(p)
    thresh = q * np.arange(1, m + 1) / m
    ok = p[o] <= thresh
    k = np.flatnonzero(ok).max() + 1 if ok.any() else 0
    rej = np.zeros(m, bool)
    if k:
        rej[o[:k]] = True
    return rej, (p[o[k - 1]] if k else 0.0)


families = {"missing": ((~Y).astype(float), np.ones_like(Y, dtype=float)),
            "depth>=6": ((dep >= 6).astype(float),
                         (Y & ((dep >= 6) | (trm != JUNK))).astype(float))}
samp_pos = [np.flatnonzero(g_samp == s) for s in range(Gs)]
out = {"arm": arm, "screen": SCREEN, "n_cells": int(n), "n_perm": B,
       "fdr_q": QFDR, "families": {}}
draws = {}

for name, (M, W) in families.items():
    P = fit_twoway(M, W)
    r = np.where(W > 0, M - P, 0.0)
    v = np.where(W > 0, P * (1 - P), 0.0)
    sv = np.sqrt(v)
    nw = pair_sums(r * r, r, g_clone, Gc)
    dw = pair_sums(v, sv, g_clone, Gc)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_obs = np.where(dw > 0, nw / np.maximum(dw, 1e-300), np.nan)
    null = np.empty((B, K))
    for b in range(B):
        perm = np.arange(n)
        for pos in samp_pos:
            perm[pos] = rng.permutation(pos)
        rp, vp = r[perm], v[perm]
        nb = pair_sums(rp * rp, rp, g_clone, Gc)
        db = pair_sums(vp, np.sqrt(vp), g_clone, Gc)
        null[b] = np.where(db > 0, nb / np.maximum(db, 1e-300), np.nan)
        if (b + 1) % max(B // 5, 1) == 0:
            print(f"    {name}: {b+1:,}/{B:,}", flush=True)
    mu, sd = np.nanmean(null, 0), np.nanstd(null, 0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        zs = (rho_obs - mu) / np.where(sd > 0, sd, np.nan)
        skew = np.nanmean(((null - mu) / np.where(sd > 0, sd, np.nan)) ** 3, 0)
    p_up = (1 + (null >= rho_obs[None, :]).sum(0)) / (B + 1)
    p_dn = (1 + (null <= rho_obs[None, :]).sum(0)) / (B + 1)
    ok = np.isfinite(rho_obs) & (sd > 0)
    rej = np.zeros(K, bool); cut = 0.0
    if ok.any():
        rj, cut = bh(p_up[ok], QFDR)
        rej[np.flatnonzero(ok)] = rj
    draws[f"{name}_null"] = null.astype(np.float32)
    out["families"][name] = {
        "marginal": ((M * W).sum(0) / np.maximum(W.sum(0), 1)).tolist(),
        "info_denom": dw.tolist(), "rho_obs": np.nan_to_num(rho_obs).tolist(),
        "null_mean": np.nan_to_num(mu).tolist(), "null_sd": np.nan_to_num(sd).tolist(),
        "z": np.nan_to_num(zs).tolist(), "null_skew": np.nan_to_num(skew).tolist(),
        "p_upper": p_up.tolist(), "excess": np.nan_to_num(rho_obs - mu).tolist(),
        "testable": int(ok.sum()), "n_sig_fdr": int(rej.sum()),
        "n_sig_low": int((p_dn[ok] <= QFDR / max(ok.sum(), 1)).sum()),
        "bh_cut": float(cut), "sig": rej.tolist(),
        "pooled_obs": float(np.nansum(nw) / np.nansum(dw))}
    d = out["families"][name]
    print(f"\n  {name}: testable {d['testable']}/{K} | "
          f"significant at FDR {QFDR}: {d['n_sig_fdr']} "
          f"({100*d['n_sig_fdr']/max(d['testable'],1):.0f}% of testable) | "
          f"significantly NEGATIVE: {d['n_sig_low']} (should be ~0)")
    zf = np.array(d["z"])[ok]
    print(f"     z: median {np.median(zf):+.1f}, 90th {np.quantile(zf,0.9):+.1f}, "
          f"max {zf.max():+.1f} | null skew median {np.median(np.array(d['null_skew'])[ok]):+.2f}",
          flush=True)

tag = "_screen" if SCREEN else ""
(RES / f"pertape_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
np.savez_compressed(RES / f"pertape_null_draws_{arm}{tag}.npz", **draws)
print(f"\nwrote results/pertape_{arm}{tag}.json")
