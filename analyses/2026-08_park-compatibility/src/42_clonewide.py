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

B = 1,000 (2026-09-03): one scan here costs ~0.02-0.12 s, so this layer needs no
parallelisation -- run it with --nperm 1000 directly.  It reports the same two
quantities as the sub-clone catalogue: a global permutation p-value at the
chosen threshold and a per-(clone,tape) q-value attached to every output row.

Usage: 42_clonewide.py <arm> [--eps 0.01] [--nperm 1000] [--lam 10] [--seed S]
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
SEED = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260903

# same fixed grid as scripts 40 and 43, so the layers are directly comparable
GRID = np.unique(np.round(np.geomspace(LAM0, 1.0e7, 600), 4))


def counts_ge(v):
    """#{v >= t} for every t in GRID."""
    sv = np.sort(v)
    return (sv.size - np.searchsorted(sv, GRID, side="left")).astype(np.int64)

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
grid = GRID
on = counts_ge(obs).astype(float)
# per-permutation count vectors, not candidate lists: B = 1,000 x a few hundred
# integers, and both the mean (for q) and the exceedance count (for p) come out.
null_counts = np.empty((NPERM, grid.size))
for b in range(NPERM):
    rb = np.random.default_rng([SEED, b])       # permutation b independent of B
    perm = np.arange(n)
    for pos in samp_pos:
        perm[pos] = rb.permutation(pos)
    Ln = lam_all(W[perm])
    v = Ln[np.isfinite(Ln)]
    null_counts[b] = counts_ge(v[v >= LAM0])
    if (b + 1) % max(NPERM // 10, 1) == 0:
        print(f"    perm {b+1}/{NPERM}", flush=True)
B = NPERM
nn = null_counts.mean(0)

# BH-style monotonisation: running minimum UP the grid (see 40_event_catalogue.py)
fdr_raw = np.where(on > 0, nn / np.maximum(on, 1.0), 1.0)
fdrc = np.minimum(np.minimum.accumulate(fdr_raw), 1.0)
ok = np.flatnonzero(fdrc <= 0.05)
j0 = int(ok[0]) if ok.size else grid.size - 1
LAM = float(grid[j0]); FDR = float(fdrc[j0])
p_at = lambda j: (1 + int((null_counts[:, j] >= on[j]).sum())) / (B + 1)
p_thresh, p_floor = p_at(j0), p_at(0)
show = np.searchsorted(grid, [10, 15, 20, 30, 50, 75, 100, 200, 500, 1000])
show = np.unique(np.clip(show, 0, len(grid) - 1))
print(f"  candidates >= {LAM0}: {obs.size:,} | null mean {nn[0]:,.1f}")
print(f"  FDR curve: " + "  ".join(f"{grid[i]:.0f}n:{100*fdrc[i]:.0f}%" for i in show))
print(f"  => threshold {LAM:.1f} nats at FDR {100*FDR:.1f}%")
print(f"  global permutation p (B={B:,}): {p_thresh:.4g} at the threshold, "
      f"{p_floor:.4g} at the scan floor {LAM0:.0f} nats", flush=True)

gg, zz = np.nonzero(L >= LAM)
mc = np.array([miss[g == c, t].sum() for c, t in zip(gg, zz)], float)
ec = np.array([P[g == c, t].sum() for c, t in zip(gg, zz)], float)
szc = sizes[gg].astype(float)
tot_missing = int(miss.sum())
qrow = fdrc[np.clip(np.searchsorted(grid, L[gg, zz], side="right") - 1, 0, grid.size - 1)]
out = {"arm": arm, "eps": EPS, "lambda_threshold": LAM, "fdr": FDR, "n_perm": NPERM,
       "seed": SEED, "p_global_at_threshold": p_thresh,
       "p_global_at_scan_floor": p_floor, "p_floor_resolution": 1.0 / (B + 1),
       "null_mean_at_threshold": float(nn[j0]),
       "null_max_at_threshold": float(null_counts[:, j0].max()),
       "null_mean_at_scan_floor": float(nn[0]),
       "null_max_at_scan_floor": float(null_counts[:, 0].max()),
       "max_q_of_kept": float(qrow.max()) if qrow.size else 1.0,
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
    fh.write("clone\tclone_bc\ttape\tclone_cells\tn_missing\texpected\tinside_rate\t"
             "lambda_nats\tq_value\n")
    for i in range(gg.size):
        fh.write(f"{gg[i]}\t{cl_names[gg[i]]}\t{zz[i]}\t{int(szc[i])}\t{int(mc[i])}\t"
                 f"{ec[i]:.3f}\t{mc[i]/szc[i]:.4f}\t{L[gg[i], zz[i]]:.3f}\t"
                 f"{qrow[i]:.3g}\n")
print(f"wrote results/clonewide_{arm}.json")
