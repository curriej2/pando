#!/usr/bin/env python3
r"""
How much of Lambda_soft is just "we fitted a parameter"?  -- measured, not argued.

THE TENSION (Justin, 2026-09-03).  pi_hat = k/m is fitted from the very clade
being tested, so it matches that clade's observed count EXACTLY.  Is it then
almost guaranteed to beat H0, making Lambda_soft uninformative?

The answer has to be empirical: compute Lambda_soft for EVERY candidate (clade,
tape) -- not only those above threshold -- in the real data and under within-clone
permutation, and compare the two whole distributions.  If the fitted-parameter
advantage were driving the statistic, the two would look alike.

For reference, if H0 were nested in H1 (it is not, exactly -- H0 gives each cell
its own p~, H1 one shared pi), Wilks would put the null at 2*Lambda ~ chi^2_1,
i.e. mean Lambda = 0.5 nats.  Measuring where the null actually sits is the point.

Output: results/lrt_calibration_{arm}.json
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4
arm = sys.argv[1]
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 3
MAXD = int(sys.argv[sys.argv.index("--maxd") + 1]) if "--maxd" in sys.argv else 4
DEPTHS = list(range(1, MAXD + 1))
rng = np.random.default_rng(20260903)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s], clone[s]
codes = np.load(RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4 else
                       f"prefix_codes_{arm}.npz"), allow_pickle=False)["codes"]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
miss = ~Y
print(f"{arm}: {n:,} cells, {K} tapes, {Gc:,} clones, depths 1-{MAXD}", flush=True)


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
eta = np.log(P / (1 - P))
o0 = np.argsort(g_clone, kind="stable")
b0 = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])
for a_, b_ in zip(b0[:-1], b0[1:]):
    ix = o0[a_:b_]
    e = eta[ix]; tgt = miss[ix].sum(0).astype(float); gam = np.zeros(K)
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)

MISS = miss.astype(float)
A1 = MISS * np.log(P)
A2 = (1 - MISS) * np.log1p(-P)
EDGES = np.array([-np.inf, -1, -0.25, 0, 0.25, 0.5, 1, 2, 4, 8, 10, 16, 25, 50, 100, np.inf])


def scan_all(mi, a1, a2):
    """Every candidate's Lambda_soft, no threshold. Returns (hist, sum, sumsq, n, tail)."""
    hist = np.zeros(len(EDGES) - 1); tot = 0; s1 = 0.0; s2 = 0.0
    tails = {t: 0 for t in (2, 4, 10, 16, 25)}
    for d in DEPTHS:
        for a_ in range(K):
            cd = codes[:, a_, d - 1]
            ok = cd >= 0
            if ok.sum() < MIN_CLADE:
                continue
            idx = np.flatnonzero(ok)
            key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
            _, sub = np.unique(key, return_inverse=True)
            G = sub.max() + 1
            m = np.bincount(sub, minlength=G).astype(float)
            k = np.empty((G, K)); A = np.empty((G, K)); B = np.empty((G, K))
            for t in range(K):
                k[:, t] = np.bincount(sub, weights=mi[idx, t], minlength=G)
                A[:, t] = np.bincount(sub, weights=a1[idx, t], minlength=G)
                B[:, t] = np.bincount(sub, weights=a2[idx, t], minlength=G)
            mm = m[:, None]
            with np.errstate(divide="ignore", invalid="ignore"):
                t1 = np.where(k > 0, k * np.log(np.maximum(k, 1e-300) / mm), 0.0)
                t2 = np.where(mm - k > 0, (mm - k) * np.log(np.maximum(mm - k, 1e-300) / mm), 0.0)
            ls = t1 + t2 - (A + B)
            keep = (m >= MIN_CLADE)[:, None] & np.ones((1, K), bool)
            keep[:, a_] = False
            v = ls[keep]
            hist += np.histogram(v, bins=EDGES)[0]
            tot += v.size; s1 += v.sum(); s2 += (v * v).sum()
            for t in tails:
                tails[t] += int((v >= t).sum())
    return hist, s1, s2, tot, tails


h, s1, s2, tot, tl = scan_all(MISS, A1, A2)
mean, sd = s1 / tot, np.sqrt(s2 / tot - (s1 / tot) ** 2)
print(f"\n  OBSERVED: {tot:,} candidate (clade,tape) pairs, "
      f"mean Lambda_soft {mean:+.3f}, sd {sd:.3f}", flush=True)

nh, nm, nsd, nt, ntl = np.zeros_like(h), [], [], 0, {t: [] for t in tl}
for b in range(NPERM):
    pm = o0.copy()
    for a_, b_ in zip(b0[:-1], b0[1:]):
        pm[a_:b_] = rng.permutation(pm[a_:b_])
    inv = np.empty_like(pm); inv[o0] = pm
    hh, q1, q2, qt, qtl = scan_all(MISS[inv], A1[inv], A2[inv])
    nh += hh / NPERM; nm.append(q1 / qt)
    nsd.append(np.sqrt(q2 / qt - (q1 / qt) ** 2)); nt = qt
    for t in tl:
        ntl[t].append(qtl[t])
    print(f"    perm {b+1}/{NPERM}: mean {q1/qt:+.3f}, sd {np.sqrt(q2/qt-(q1/qt)**2):.3f}",
          flush=True)

out = {"arm": arm, "max_depth": MAXD, "n_candidates": int(tot),
       "observed_mean": float(mean), "observed_sd": float(sd),
       "null_mean": float(np.mean(nm)), "null_sd": float(np.mean(nsd)),
       "wilks_reference_mean": 0.5,
       "edges": [float(x) for x in EDGES],
       "observed_hist": h.tolist(), "null_hist": nh.tolist(),
       "tails": {str(t): {"observed": tl[t], "null": float(np.mean(ntl[t])),
                          "enrichment": tl[t] / max(np.mean(ntl[t]), 1e-9)} for t in tl}}
print(f"\n  null mean Lambda_soft {out['null_mean']:+.3f} (Wilks reference for a "
      f"1-parameter fit: +0.500), sd {out['null_sd']:.3f}")
print(f"\n  {'Lambda >=':>10} {'observed':>12} {'null':>12} {'enrichment':>12}")
for t in (2, 4, 10, 16, 25):
    d = out["tails"][str(t)]
    print(f"  {t:>10} {d['observed']:>12,} {d['null']:>12,.0f} {d['enrichment']:>11.1f}x")
(RES / f"lrt_calibration_{arm}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/lrt_calibration_{arm}.json")
