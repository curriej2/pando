#!/usr/bin/env python3
r"""
Diagnostic: is the soft route's low pi_hat a real graded effect, or is it COARSE
ATTRIBUTION -- the overlap collapse crediting a big shallow clade when the loss
actually sits on a small deep one?

Rebuilds 63's above-threshold set, then compares two collapse rules on it:
  (a) by MARGIN      (what 63 does)  -- detection statistic used for attribution
  (b) by Lambda_hard (the Dollo rule) -- by the identity Lambda_hard = Lambda_soft
      - m KL(pi||1-eps), this is maximised at the LARGEST clade that is still
      COMPLETE, which is where a Dollo loss on a stem belongs.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4
arm = sys.argv[1]; MAXD = int(sys.argv[2]) if len(sys.argv) > 2 else 6
BUDGET = 2.0
tag = f"{arm}_d{MAXD}"

parts = sorted((RES / "percombo_soft").glob(f"percombo_soft_{tag}_p*.npz"))
obs_s = obs_h = obs_k = size = valid = mx_s = None
B = 0
for f in parts:
    z = np.load(f, allow_pickle=False)
    if obs_s is None:
        obs_s, obs_h, obs_k = z["obs_s"], z["obs_h"], z["obs_k"]
        size, valid = z["size"], z["valid"]; mx_s = z["mx_s"].copy()
    else:
        np.maximum(mx_s, z["mx_s"], out=mx_s)
    B += int(z["meta"][4]); del z
zc = np.load(RES / "percombo_soft" / f"percombo_soft_{tag}_cal.npz", allow_pickle=False)
LSC = zc["lam_s"]; NCAL = LSC.shape[0]
V = valid & np.isfinite(obs_s) & (mx_s > -1e29)
m_obs = np.where(V, obs_s - mx_s, -np.inf)
m_cal = [np.where(V, LSC[j] - mx_s, -np.inf) for j in range(NCAL)]
pi = np.where(size > 0, obs_k.astype(np.float64) / np.maximum(size, 1), np.nan)

SIZE_EDGES = np.array([4, 6, 10, 20, 50, 200, 10 ** 9])
SIZE_LABEL = ["4-5", "6-9", "10-19", "20-49", "50-199", "200+"]
strat = np.clip(np.searchsorted(SIZE_EDGES, size, side="right") - 1, 0, 5)
grid = np.unique(np.round(np.concatenate(
    [np.arange(0.0, 5.0, 0.05), np.geomspace(5.0, 2000.0, 300)]), 4))
THRS = np.full(6, np.inf)
for i in range(6):
    m = V & (strat == i)
    if not m.any(): continue
    u = np.array([np.mean([(v[m] >= t).sum() for v in m_cal]) for t in grid])
    o = np.array([(m_obs[m] >= t).sum() for t in grid], float)
    ok = np.flatnonzero(u <= BUDGET)
    if ok.size and o[int(ok[0])] > 0: THRS[i] = grid[int(ok[0])]
THRM = THRS[strat]

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s_ = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s_], clone[s_]
codes = np.load(RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4
                       else f"prefix_codes_{arm}.npz"), allow_pickle=False)["codes"]
n, K = Y.shape
cl_names, g_clone = np.unique(clone, return_inverse=True)
miss = ~Y
hits, off = [], 0
for d in range(1, MAXD + 1):
    for a_ in range(K):
        cd = codes[:, a_, d - 1]; ok = cd >= 0
        if ok.sum() < MIN_CLADE: continue
        idx = np.flatnonzero(ok)
        key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
        _, sub = np.unique(key, return_inverse=True)
        G = int(sub.max() + 1); ns = np.bincount(sub, minlength=G)
        keep = np.flatnonzero(ns >= MIN_CLADE)
        if keep.size == 0: continue
        sl0 = slice(off, off + keep.size * K)
        blk = m_obs[sl0].reshape(keep.size, K)
        gg, zz = np.nonzero(blk >= THRM[sl0].reshape(keep.size, K))
        if gg.size:
            o2 = np.argsort(sub, kind="stable"); cb = np.concatenate([[0], np.cumsum(ns)])
            owner = np.zeros(G, np.int64); owner[sub] = g_clone[idx]
            for gi, z_ in zip(gg, zz):
                gk = int(keep[gi]); cells = idx[o2[cb[gk]:cb[gk + 1]]]
                sl = off + gi * K + z_
                hits.append(dict(clone=int(owner[gk]), depth=d, tape=int(z_),
                                 size=int(ns[gk]), margin=float(m_obs[sl]),
                                 lam_h=float(obs_h[sl]), pi=float(pi[sl]),
                                 nmiss=int(obs_k[sl]), cells=frozenset(cells.tolist())))
        off += keep.size * K
print(f"{tag}: {len(hits):,} combos above their stratum threshold")


def collapse(hs, keyfn):
    hs = sorted(hs, key=keyfn)
    kept, seen = [], {}
    for h in hs:
        k_ = (h["clone"], h["tape"])
        if any(h["cells"] & prev for prev in seen.get(k_, [])): continue
        seen.setdefault(k_, []).append(h["cells"]); kept.append(h)
    return kept


tot = int(miss.sum())
print(f"\n{'collapse rule':>22} {'events':>8} {'cell x tape':>12} {'% missing':>10} "
      f"{'median pi':>10} {'pi>=0.99':>9} {'partial':>8} {'med depth':>10}")
for name, kf in [("by MARGIN (63 as built)", lambda h: -h["margin"]),
                 ("by Lambda_hard (Dollo)", lambda h: -h["lam_h"])]:
    kp = collapse(hits, kf)
    p_ = np.array([h["pi"] for h in kp]); cx = sum(h["nmiss"] for h in kp)
    dp = np.array([h["depth"] for h in kp])
    print(f"{name:>22} {len(kp):>8,} {cx:>12,} {100*cx/tot:>9.2f}% {np.median(p_):>10.3f} "
          f"{100*(p_>=0.99).mean():>8.1f}% {100*(p_<0.90).mean():>7.1f}% {np.median(dp):>10.1f}")

# ---- does a shallow partial hit CONTAIN a deeper complete hit for the same (clone,tape)?
from collections import defaultdict
grp = defaultdict(list)
for h in hits: grp[(h["clone"], h["tape"])].append(h)
nested = both = 0
for k_, hs in grp.items():
    shallow = [h for h in hs if h["pi"] < 0.90]
    deep = [h for h in hs if h["pi"] >= 0.99]
    if shallow and deep:
        both += 1
        if any(d["cells"] < s["cells"] for s in shallow for d in deep): nested += 1
print(f"\n(clone,tape) groups above threshold: {len(grp):,}")
print(f"  groups holding BOTH a partial (pi<0.90) and a complete (pi>=0.99) combo: {both:,}")
print(f"  ...of which the complete one is NESTED INSIDE the partial one: {nested:,}")
