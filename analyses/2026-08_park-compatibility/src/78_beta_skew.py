#!/usr/bin/env python3
r"""
IS THE CLONE-SPECIFIC TAPE COEFFICIENT RIGHT-SKEWED, AS ONE-DIRECTIONAL LOSS PREDICTS?

THE IDEA, and it needs no event definition, no threshold and no clade search.
Inside clone C the M3 fit gives a tape offset beta_{C,z}.  Centre it across clones
for each tape to remove global tape quality:

    gamma_{C,z} = beta_{C,z} - (weighted mean over clones of beta_{.,z})

WHAT THE TWO HYPOTHESES PREDICT ABOUT ITS SHAPE
  technical dropout only -- gamma is sampling noise about 0, and SYMMETRIC: there
      is no reason a clone should be systematically worse OR better at a tape.
  heritable silencing    -- loss is ONE-DIRECTIONAL (Dollo: a silenced integration
      is not un-silenced), so a clone either never lost tape z (gamma ~ 0) or lost
      it in some subclade (gamma > 0).  Never gamma < 0.  ⇒ a spike at zero with a
      HEAVY RIGHT TAIL, i.e. positive skew.

⚠ Skewness is location-invariant, so the centring convention cannot create it.

THE NULL.  Permute cell membership ACROSS clones within the arm, preserving clone
sizes, then refit everything.  Whole cell ROWS move, so each cell's profile, every
tape's marginal and every clone size are preserved, and the only thing destroyed is
which cells share a clone -- exactly the structure under test.  Any asymmetry the
ESTIMATOR itself induces (logit nonlinearity, finite n_C, the +-12 clip) appears
identically in the null and cancels in the comparison.

⚠ No A/B tape split here: nothing in this test involves relatedness, so all 166
tapes are used and beta is estimated as well as possible.

Usage: 78_beta_skew.py <arm> [--mincl 100] [--nperm 200] [--iters 60] [--seed 11]
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
CLIP = 12.0


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
MINCL, NPERM = arg("--mincl", 100, int), arg("--nperm", 200, int)
ITERS, SEED = arg("--iters", 60, int), arg("--seed", 11, int)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
_, g0 = np.unique(clone, return_inverse=True)
sizes = np.bincount(g0)
keep = np.isin(g0, np.flatnonzero(sizes >= MINCL))
X = (~Y[keep]).astype(np.float64)
_, g = np.unique(g0[keep], return_inverse=True)
o = np.argsort(g, kind="stable"); X, g = X[o], g[o]
n, K = X.shape; nC = int(g.max()) + 1
off = np.concatenate([[0], np.cumsum(np.bincount(g))]); csz = np.diff(off)
print(f"{arm}: {n:,} cells, {nC} clones (>= {MINCL}), {K} tapes, {NPERM} permutations",
      flush=True)
sig = lambda e: 1.0 / (1.0 + np.exp(-e))


def gammas(Xa):
    """Fit alpha_c + beta_z inside every clone, then centre beta across clones."""
    B = np.zeros((nC, K))
    for b_ in range(nC):
        Xb = Xa[off[b_]:off[b_ + 1]]
        m = Xb.shape[0]
        r0 = np.clip(Xb.mean(0), 1e-3, 1 - 1e-3)
        be = np.log(r0 / (1 - r0)); al = np.zeros(m)
        for _ in range(ITERS):
            p = sig(al[:, None] + be[None, :])
            gr = (Xb - p).sum(1); he = (p * (1 - p)).sum(1)
            al = np.clip(al + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -CLIP, CLIP)
            al -= al.mean()
            p = sig(al[:, None] + be[None, :])
            gr = (Xb - p).sum(0); he = (p * (1 - p)).sum(0)
            be = np.clip(be + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -CLIP, CLIP)
        B[b_] = be
    w = csz / csz.sum()
    return B - (w[:, None] * B).sum(0)[None, :]


def stats(G):
    v = G.ravel()
    s = v.std()
    g1 = float(((v - v.mean()) ** 3).mean() / max(s, 1e-12) ** 3)
    q = np.quantile(v, [0.01, 0.25, 0.50, 0.75, 0.95, 0.99])
    # quantile skew, robust to the clip: (q99-q50) vs (q50-q01)
    qs = float((q[5] - q[2]) / max(q[2] - q[0], 1e-12))
    return dict(g1=g1, qskew=qs, sd=float(s),
                tail_pos=float(np.mean(v > 2.0)), tail_neg=float(np.mean(v < -2.0)),
                clip_hi=float(np.mean(v >= CLIP - 1e-6)), clip_lo=float(np.mean(v <= -CLIP + 1e-6)),
                q=[float(x) for x in q])


t0 = time.time()
obs = stats(gammas(X))
print(f"  observed: g1 {obs['g1']:+.3f}  qskew {obs['qskew']:.3f}  sd {obs['sd']:.3f}  "
      f"P(g>2) {100*obs['tail_pos']:.2f}%  P(g<-2) {100*obs['tail_neg']:.2f}%  "
      f"[{time.time()-t0:.0f}s]", flush=True)
nul = []
for b in range(NPERM):
    rb = np.random.default_rng([SEED, b])
    Xp = X[rb.permutation(n)]              # whole cell ROWS move; clone sizes fixed
    nul.append(stats(gammas(Xp)))
    if (b + 1) % 25 == 0:
        print(f"  perm {b+1}/{NPERM}  g1 {np.mean([x['g1'] for x in nul]):+.3f} "
              f"[{time.time()-t0:.0f}s]", flush=True)


def cmp(key, bigger=True):
    o_ = obs[key]; nn = np.array([x[key] for x in nul])
    ex = int((nn >= o_).sum() if bigger else (nn <= o_).sum())
    p = (1 + ex) / (len(nn) + 1)
    zz = (o_ - nn.mean()) / max(nn.std(ddof=1), 1e-12)
    return o_, nn.mean(), nn.std(ddof=1), nn.max() if bigger else nn.min(), zz, p


print(f"\n  {'statistic':>12}{'observed':>11}{'null mean':>11}{'null sd':>10}"
      f"{'null max':>10}{'z':>9}{'p':>9}")
for k in ("g1", "qskew", "tail_pos", "sd"):
    o_, m_, s_, x_, zz, p = cmp(k)
    print(f"  {k:>12}{o_:>+11.4f}{m_:>+11.4f}{s_:>10.4f}{x_:>+10.4f}{zz:>9.1f}{p:>9.4f}")
o_, m_, s_, x_, zz, p = cmp("tail_neg", bigger=False)
print(f"  {'tail_neg':>12}{o_:>+11.4f}{m_:>+11.4f}{s_:>10.4f}{x_:>+10.4f}{zz:>9.1f}{p:>9.4f}"
      f"   (tested for being SMALLER)")
print(f"\n  asymmetry of the observed tails: P(g>2)/P(g<-2) = "
      f"{obs['tail_pos']/max(obs['tail_neg'],1e-12):.2f}   "
      f"null {np.mean([x['tail_pos']/max(x['tail_neg'],1e-12) for x in nul]):.2f}")
print(f"  clipped at +{CLIP:.0f}: obs {100*obs['clip_hi']:.3f}%  "
      f"null {100*np.mean([x['clip_hi'] for x in nul]):.3f}%   |  at -{CLIP:.0f}: "
      f"obs {100*obs['clip_lo']:.3f}%  null {100*np.mean([x['clip_lo'] for x in nul]):.3f}%")
print(f"  quantiles (1,25,50,75,95,99): obs {['%.2f'%v for v in obs['q']]}")
print(f"                                null {['%.2f'%v for v in np.mean([x['q'] for x in nul],0)]}")
(RES / f"betaskew_{arm}_mincl{MINCL}.json").write_text(
    json.dumps(dict(arm=arm, mincl=MINCL, n_cells=int(n), n_clones=int(nC), K=int(K),
                    nperm=NPERM, observed=obs, null=nul), indent=1))
print(f"\nwrote results/betaskew_{arm}_mincl{MINCL}.json")
