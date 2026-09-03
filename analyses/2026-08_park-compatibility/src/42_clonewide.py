#!/usr/bin/env python3
r"""
The CLONE-WIDE layer: tapes lost across an entire clone.

40_event_catalogue.py deliberately removes this layer (the gamma_{C,z} offset
absorbs it) so that it can ask about structure BELOW the clone.  This script
measures the layer that was removed, with the same statistic one level up:

    Lambda_clone(C,z) = SUM_{c in C, missing} log[(1-eps)/p_hat_cz]
                      + SUM_{c in C, present} log[eps/(1 - p_hat_cz)]

with the clade = the whole clone and p_hat the TWO-margin fit (alpha_c + beta_z,
no gamma -- gamma is exactly what we now want to measure rather than absorb).

NULL: cells permuted ACROSS clones WITHIN sample, destroying clone identity while
preserving each cell's whole missingness profile and each tape's marginal.  Same
FDR-curve machinery as script 40.

⚠ THE LAYERS DO NOT ADD TO SHARES OF A TOTAL.  Each fit matches its own margins,
so every layer sums to zero against the layer beneath it.  Report each layer
against the one below: Lambda in nats (the likelihood a model gains) and the count
of missing entries inside confidently-called blocks.  Do NOT present these as a
partition of missingness.

Usage: 42_clonewide.py <arm> [--eps 0.01] [--nperm 200] [--lam 10]
Output: results/clonewide_{arm}.json, results/clonewide_{arm}.tsv.gz
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLONE = 5

arm = sys.argv[1]
EPS = float(sys.argv[sys.argv.index("--eps") + 1]) if "--eps" in sys.argv else 0.01
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 200
LAM0 = float(sys.argv[sys.argv.index("--lam") + 1]) if "--lam" in sys.argv else 10.0
rng = np.random.default_rng(20260903)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone, sample = z0["recovered"], z0["clone"], z0["sample"]
s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone, sample = Y[s], clone[s], sample[s]
n, K = Y.shape
cl_names, g = np.unique(clone, return_inverse=True)
_, gs = np.unique(sample, return_inverse=True)
Gc, Gs = g.max() + 1, gs.max() + 1
miss = ~Y
sizes = np.bincount(g)
print(f"{arm}: {n:,} cells, {K} tapes, {Gc:,} clones, {Gs} samples", flush=True)


def fit_twoway(M, iters=40):
    cm = np.clip(M.mean(0), 1e-3, 1 - 1e-3)
    beta = np.log(cm / (1 - cm)); alpha = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))
            if ax == 0:
                gg = (M - p).sum(1); h = (p * (1 - p)).sum(1)
                alpha = np.clip(alpha + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
                alpha -= alpha.mean()
            else:
                gg = (M - p).sum(0); h = (p * (1 - p)).sum(0)
                beta = np.clip(beta + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :]))), 1e-6, 1 - 1e-6)


P = fit_twoway(miss.astype(float))
W = np.where(miss, np.log(1 - EPS) - np.log(P), np.log(EPS) - np.log1p(-P))
big = sizes >= MIN_CLONE


def lam_all(Wm):
    L = np.empty((Gc, K))
    for k in range(K):
        L[:, k] = np.bincount(g, weights=Wm[:, k], minlength=Gc)
    L[~big] = -np.inf
    return L


L = lam_all(W)
obs = L[np.isfinite(L)]
obs = obs[obs >= LAM0]
samp_pos = [np.flatnonzero(gs == i) for i in range(Gs)]
null = []
for b in range(NPERM):
    perm = np.arange(n)
    for pos in samp_pos:
        perm[pos] = rng.permutation(pos)
    Ln = lam_all(W[perm])
    v = Ln[np.isfinite(Ln)]
    null.append(v[v >= LAM0])
    if (b + 1) % max(NPERM // 4, 1) == 0:
        print(f"    perm {b+1}/{NPERM}", flush=True)

grid = np.unique(np.round(np.geomspace(LAM0, max(obs.max(), LAM0 * 2), 60), 2))
on = np.array([(obs >= t).sum() for t in grid], float)
nn = np.array([np.mean([(v >= t).sum() for v in null]) for t in grid])
fdrc = np.where(on > 0, nn / np.maximum(on, 1), 1.0)
ok = np.flatnonzero(fdrc <= 0.05)
LAM = float(grid[ok[0]]) if ok.size else float(grid[-1])
FDR = float(fdrc[ok[0]]) if ok.size else float(fdrc[-1])
print(f"  candidates >= {LAM0}: {obs.size:,} | null {np.mean([v.size for v in null]):,.0f}")
print(f"  FDR curve: " + "  ".join(f"{t:.0f}n:{100*f:.0f}%" for t, f in zip(grid[::6], fdrc[::6])))
print(f"  => threshold {LAM:.1f} nats at FDR {100*FDR:.1f}%", flush=True)

gg, zz = np.nonzero(L >= LAM)
mc = np.array([miss[g == c, t].sum() for c, t in zip(gg, zz)], float)
ec = np.array([P[g == c, t].sum() for c, t in zip(gg, zz)], float)
szc = sizes[gg].astype(float)
tot_missing = int(miss.sum())
out = {"arm": arm, "eps": EPS, "lambda_threshold": LAM, "fdr": FDR, "n_perm": NPERM,
       "n_clone_tape_pairs_tested": int(big.sum()) * K,
       "n_clonewide_losses": int(gg.size),
       "distinct_tapes": int(len(set(zz.tolist()))),
       "distinct_clones": int(len(set(gg.tolist()))),
       "slots": float(szc.sum()), "missing_in_losses": float(mc.sum()),
       "expected_missing": float(ec.sum()),
       "share_of_all_missing": float(mc.sum() / max(tot_missing, 1)),
       "excess_over_twoway": float(mc.sum() - ec.sum()),
       "lambda_total_nats": float(L[gg, zz].sum()),
       "median_inside_rate": float(np.median(mc / szc)) if gg.size else None,
       "median_expected_rate": float(np.median(ec / szc)) if gg.size else None,
       "total_missing_entries": tot_missing,
       "fdr_curve": {"lambda": grid.tolist(), "observed": on.tolist(),
                     "null": nn.tolist(), "fdr": fdrc.tolist()}}
print(f"\n  clone-wide losses: {gg.size:,} (clone,tape) pairs on "
      f"{out['distinct_tapes']} tapes in {out['distinct_clones']} clones")
print(f"  inside rate median {out['median_inside_rate']} vs expected "
      f"{out['median_expected_rate']}")
print(f"  missing entries inside them: {mc.sum():,.0f} = "
      f"{100*out['share_of_all_missing']:.2f}% of all missing")
print(f"  Lambda total {out['lambda_total_nats']:,.0f} nats", flush=True)
(RES / f"clonewide_{arm}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"clonewide_{arm}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\ttape\tclone_cells\tn_missing\texpected\tinside_rate\tlambda_nats\n")
    for i in range(gg.size):
        fh.write(f"{gg[i]}\t{cl_names[gg[i]]}\t{zz[i]}\t{int(szc[i])}\t{int(mc[i])}\t"
                 f"{ec[i]:.3f}\t{mc[i]/szc[i]:.4f}\t{L[gg[i], zz[i]]:.3f}\n")
print(f"wrote results/clonewide_{arm}.json")
