#!/usr/bin/env python3
r"""
================================================================================
 T0 / T1 -- do related cells lose the SAME tapes?      (row A9)
================================================================================

THE QUESTION.  Two cells of one clone, one tape: are they missing it together
more often than two unrelated cells of the same harvest sample?

THE STATISTIC.  For cell c and feature f (a tape, for missingness) write

    r_cf = Y_cf - p_hat_cf          the SURPRISE: what happened minus what we
                                    expected from the cell and the feature alone
    v_cf = p_hat_cf (1 - p_hat_cf)  the variance that surprise should have had

Multiply the surprises of two cells sharing a group and average: positive means
they went the same way relative to expectation.  Normalised by the typical size of
such products, that is a correlation -- the intraclass correlation of the feature
within the group.  Pairs are never enumerated; for any group g,

    sum_{c != c' in g} r_c r_c'  =  (sum_c r_c)^2 - sum_c r_c^2

so two running totals per (group, feature) deliver every pair.  Hence

    rho_within = SUM_g [ (S1)^2 - S2 ]  /  SUM_g [ (T1)^2 - T2 ]
    S1 = sum_g r, S2 = sum_g r^2, T1 = sum_g sqrt(v), T2 = sum_g v

Computed at TWO nesting levels -- clone and harvest sample -- so the between-clone
term is free:  same-sample pairs minus same-clone pairs = different clone, same
sample.  rho_within - rho_between is then the lineage signal, with batch,
cell quality and feature quality all held fixed.

THREE CENTRINGS, because the nuisances are themselves of interest
(Justin, 2026-09-02): a clone whose cells merely share a low CAPTURE RATE is a
real and separate phenomenon -- it is the birth-death uniform-sampling question --
not only a confound to be removed.  So p_hat is built three ways:

    grand    p_hat = p_bar                total within-clone agreement
    tape     p_hat = beta_z               after feature quality; still contains
                                          whole-cell capture clustering
    twoway   p_hat = sigma(alpha_c+beta_z) tape-SPECIFIC agreement only

(tape - twoway) is the capture-clustering contribution; twoway is heritable loss
of particular tapes.  rho_clone of R_c itself is reported alongside as the direct,
readable form of the same thing.

T0 -- POSITIVE CONTROL, AND IT SETS THE SCALE.  rho depends on a feature's marginal
frequency, so a single number is not a control.  We run the identical statistic on
features we KNOW are heritable -- "has tape z reached site L", L = 2..6, whose
marginals sweep ~0.95 down to ~0.4 and so overlap the missingness marginals -- and
report rho excess against marginal frequency.  Missingness is then read against
that curve.  The within-sample permutation supplies the zero line.

GROUPING IS A PARAMETER (--group).  Clone today; a prefix clade tomorrow, to test
heritability at a finer level than clonal, with no rework.

Usage: 34_lineage_rho.py <arm> [--screen] [--nperm N] [--group clone]
Output: results/lineage_rho_{arm}{_screen}.json
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
ABSENT, UNEDITED, JUNK, COMPLETE = 0, 1, 2, 3

arm = sys.argv[1]
SCREEN = "--screen" in sys.argv
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 200
GROUP = sys.argv[sys.argv.index("--group") + 1] if "--group" in sys.argv else "clone"
rng = np.random.default_rng(20260902)

# ---------------------------------------------------------------- data
z = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, dep, trm, clone, sample = (z["recovered"], z["depth"].astype(np.int16),
                              z["term"], z["clone"], z["sample"])
keep = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, dep, trm, clone, sample = Y[keep], dep[keep], trm[keep], clone[keep], sample[keep]

if SCREEN:
    fl = np.load(RES / f"collision_flags_{arm}.npz", allow_pickle=False)
    names = [str(x) for x in fl["clone_names"]]
    # flags were built on the same filter, in the same row order
    assert fl["flag"].size == Y.shape[0], "collision flags do not match the cell set"
    good = ~fl["flag"]
    Y, dep, trm, clone, sample = Y[good], dep[good], trm[good], clone[good], sample[good]
    print(f"collision screen: dropped {int((~good).sum()):,} cells")

n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
_, g_samp = np.unique(sample, return_inverse=True)
grp = g_clone if GROUP == "clone" else None
print(f"{arm}{' +screen' if SCREEN else ''}: {n:,} cells, {K} tapes, "
      f"{g_clone.max()+1:,} clones, {g_samp.max()+1} samples, grouping='{GROUP}'",
      flush=True)


# ---------------------------------------------------------------- two-way fit
def fit_twoway(M, W, iters=40):
    """p_hat = sigma(alpha_c + beta_z) by alternating Newton, on weighted entries."""
    with np.errstate(divide="ignore", invalid="ignore"):
        cm = (M * W).sum(0) / np.maximum(W.sum(0), 1)
    beta = np.log(np.clip(cm, 1e-3, 1 - 1e-3) / (1 - np.clip(cm, 1e-3, 1 - 1e-3)))
    alpha = np.zeros(M.shape[0])
    for _ in range(iters):
        for axis in (0, 1):                       # 0: update alpha, 1: update beta
            p = 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))
            g = (W * (M - p)).sum(1 - axis + 0) if axis == 0 else (W * (M - p)).sum(0)
            h = (W * p * (1 - p)).sum(1) if axis == 0 else (W * p * (1 - p)).sum(0)
            step = g / np.maximum(h, 1e-9)
            if axis == 0:
                alpha = np.clip(alpha + np.clip(step, -2, 2), -12, 12)
                alpha -= alpha.mean()
            else:
                beta = np.clip(beta + np.clip(step, -2, 2), -12, 12)
    return 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))


def group_pairs(x2, x1, g, G):
    """SUM_g [ (sum_g x1)^2 - sum_g x2 ], per feature. Returns (K,).

    One bincount per feature column -- an order of magnitude faster than
    np.add.at on the (n, K) block, and this runs inside the permutation loop."""
    K = x1.shape[1]
    out = np.empty(K)
    for k in range(K):
        S1 = np.bincount(g, weights=x1[:, k], minlength=G)
        S2 = np.bincount(g, weights=x2[:, k], minlength=G)
        out[k] = (S1 * S1 - S2).sum()
    return out


def rho_pair(M, W, mode):
    """rho_within, rho_between per feature and pooled, for one centring."""
    if mode == "grand":
        pb = (M * W).sum() / max(W.sum(), 1)
        P = np.full(M.shape, pb)
    elif mode == "tape":
        cm = (M * W).sum(0) / np.maximum(W.sum(0), 1)
        P = np.repeat(cm[None, :], M.shape[0], axis=0)
    else:
        P = fit_twoway(M, W)
    P = np.clip(P, 1e-6, 1 - 1e-6)
    r = np.where(W > 0, M - P, 0.0)
    v = np.where(W > 0, P * (1 - P), 0.0)
    sv = np.sqrt(v)
    Gc, Gs = g_clone.max() + 1, g_samp.max() + 1
    nw = group_pairs(r**2, r, g_clone, Gc); dw = group_pairs(v, sv, g_clone, Gc)
    ns = group_pairs(r**2, r, g_samp, Gs); ds = group_pairs(v, sv, g_samp, Gs)
    nb, db = ns - nw, ds - dw
    with np.errstate(divide="ignore", invalid="ignore"):
        rw = np.where(dw > 0, nw / np.maximum(dw, 1e-12), np.nan)
        rb = np.where(db > 0, nb / np.maximum(db, 1e-12), np.nan)
    return {"per_feature_within": rw, "per_feature_between": rb,
            "pooled_within": float(nw.sum() / dw.sum()),
            "pooled_between": float(nb.sum() / db.sum()),
            "marginal": (M * W).sum(0) / np.maximum(W.sum(0), 1),
            "_r": r, "_v": v}


# ---------------------------------------------------------------- features
W_all = np.ones_like(Y, dtype=float)
families = {"missing": (~Y).astype(float)}
Wt = {"missing": W_all}
for L in range(2, 7):
    # "has tape z reached site L" -- known-heritable, and its marginal (~0.95 down
    # to ~0.4) overlaps the missingness marginals, which is what makes it a control
    M = (dep >= L).astype(float)
    det = Y & ((dep >= L) | (trm != JUNK))     # junk before L: we cannot know
    families[f"depth>={L}"] = M
    Wt[f"depth>={L}"] = det.astype(float)

out = {"arm": arm, "screen": SCREEN, "n_cells": int(n), "grouping": GROUP,
       "n_clones": int(g_clone.max() + 1), "families": {}}
store = {}
for name, M in families.items():
    W = Wt[name]
    fam = {}
    for mode in ("grand", "tape", "twoway"):
        res = rho_pair(M, W, mode)
        fam[mode] = {"pooled_within": res["pooled_within"],
                     "pooled_between": res["pooled_between"],
                     "excess": res["pooled_within"] - res["pooled_between"]}
        if mode == "twoway":
            fam["marginal"] = res["marginal"].tolist()
            fam["per_feature_within"] = np.nan_to_num(res["per_feature_within"], nan=0).tolist()
            fam["per_feature_between"] = np.nan_to_num(res["per_feature_between"], nan=0).tolist()
            store[name] = (res["_r"], res["_v"])
    out["families"][name] = fam
    print(f"  {name:12s} marg {np.nanmean(res['marginal']):.3f} | " +
          "  ".join(f"{m}: w {fam[m]['pooled_within']:+.4f} b "
                    f"{fam[m]['pooled_between']:+.4f} exc {fam[m]['excess']:+.4f}"
                    for m in ("grand", "tape", "twoway")), flush=True)

# ---------------------------------------------------------------- rho_clone of R_c
Rc = Y.sum(1).astype(float)
res_c = Rc - Rc.mean()
Gc, Gs = g_clone.max() + 1, g_samp.max() + 1
def icc_cont(x, g, G):
    S1 = np.bincount(g, weights=x, minlength=G)
    S2 = np.bincount(g, weights=x**2, minlength=G)
    cnt = np.bincount(g, minlength=G).astype(float)
    return float((S1**2 - S2).sum()), float((cnt * (cnt - 1)).sum())
nw_c, dw_c = icc_cont(res_c, g_clone, Gc)
ns_c, ds_c = icc_cont(res_c, g_samp, Gs)
var = float((res_c**2).mean())
out["capture_icc"] = {"within_clone": nw_c / max(dw_c, 1) / var,
                      "between_clone_same_sample": (ns_c - nw_c) / max(ds_c - dw_c, 1) / var}
print(f"  R_c ICC: within-clone {out['capture_icc']['within_clone']:+.4f}  "
      f"between-clone/same-sample {out['capture_icc']['between_clone_same_sample']:+.4f}",
      flush=True)

# ---------------------------------------------------------------- permutation null
print(f"  permuting cells within sample, {NPERM} draws ...", flush=True)
samp_pos = [np.flatnonzero(g_samp == s) for s in range(Gs)]
null = {name: [] for name in store}
for b in range(NPERM):
    perm = np.arange(n)
    for pos in samp_pos:
        perm[pos] = rng.permutation(pos)
    for name, (r, v) in store.items():
        rp, vp = r[perm], v[perm]
        nw = group_pairs(rp**2, rp, g_clone, Gc)
        dw = group_pairs(vp, np.sqrt(vp), g_clone, Gc)
        null[name].append(float(nw.sum() / dw.sum()))
for name in null:
    a = np.array(null[name])
    obs = out["families"][name]["twoway"]["pooled_within"]
    out["families"][name]["null"] = {
        "excess_over_null": float(obs - a.mean()),
        "mean": float(a.mean()), "sd": float(a.std(ddof=1)),
        "q025": float(np.quantile(a, 0.025)), "q975": float(np.quantile(a, 0.975)),
        "z": float((obs - a.mean()) / max(a.std(ddof=1), 1e-12)),
        "p_one_sided": float((a >= obs).mean())}
    d = out["families"][name]["null"]
    print(f"  {name:12s} marg {np.mean(out['families'][name]['marginal']):.3f}  "
          f"observed {obs:+.4f}  null {d['mean']:+.4f} +/- {d['sd']:.4f}  "
          f"EXCESS {d['excess_over_null']:+.4f}  z = {d['z']:+.1f}", flush=True)

tag = "_screen" if SCREEN else ""
# ⚠ rho_between is NOT the comparator.  The two-way fit has no per-sample term, so
# the mean residual within a harvest sample is not exactly zero; squared over group
# sums of thousands of cells that small offset dominates, and it inflates the
# between term (which pools the largest groups).  The PERMUTATION null has the same
# offset by construction -- it recomputes rho_within on random within-sample groups
# of the same sizes and sample composition -- so observed-minus-null is the sound
# contrast.  rho_between is retained in the output for the record only.
out["comparator"] = "permutation null (excess_over_null); rho_between is contaminated by per-sample residual offsets"
(RES / f"lineage_rho_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/lineage_rho_{arm}{tag}.json")
