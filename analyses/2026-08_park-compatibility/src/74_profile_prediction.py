#!/usr/bin/env python3
r"""
Can a cell's dropout be predicted from its RELATIVES?  A held-out prediction, and
three readings of the same fit at three levels of digestibility.

THE SETUP.  eta_cz = alpha_c + beta_z + gamma_{C,z} is the existing model's linear
predictor (cell quality, tape quality, clone-by-tape rate), already fitted and held
fixed.  r*_cz = (X_cz - p~_cz)/sqrt(p~_cz(1-p~_cz)) is the Pearson residual.

  1. split the 166 tapes into halves A and B
  2. for each cell c, take the k most related cells IN THE SAME CLONE using ONLY
     half-A tapes, excluding c itself -> N(c)
  3. two neighbour signals at each half-B tape z, neither involving c's own data:
        u_cz = mean over N(c) of r*_c'z          (standardised; used for the fit)
        f_cz = fraction of N(c) missing tape z   (interpretable; used for the table)
  4. fit ONE scalar w on training cells:  logit P(X_cz = 1) = eta_cz + w * signal
  5. evaluate on HELD-OUT cells

THREE READINGS
  anyone         observed dropout rate by neighbour signal, WITHIN strata of the
                 model's own predicted probability.  Holding p~ fixed means cell
                 and tape quality cannot explain the contrast.
  stats-literate exp(w_f) = the odds ratio for "all relatives missing it" against
                 "none of them missing it".
  us             nats per cell (scaled to 166 tapes) -- whether it is worth adding
                 to the likelihood at all.

⚠ THE EFFECT IS CONCENTRATED, so an average understates it and an extreme stratum
overstates it.  Entry counts are printed beside every rate for that reason.
⚠ The gamma fit forces sum_c r_cz = 0 within a clone, so a random clone-mate's
residual is negatively correlated with c's by about -1/(n_C - 1) -- comparable to
the signal in a 20-cell clone.  The within-clone permutation null carries the same
constraint, so only OBSERVED MINUS NULL is quoted.

Usage: 74_profile_prediction.py <arm> [--clone C | --pooled] [--nsub 1500]
       [--nsplit 5] [--kk 5,20,50] [--mincl 20] [--seed S]
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
CL = arg("--clone", 76, int)
POOLED = "--pooled" in sys.argv
NSUB, NSPLIT = arg("--nsub", 1500, int), arg("--nsplit", 5, int)
MINCL, SEED = arg("--mincl", 20, int), arg("--seed", 11, int)
KK = [int(x) for x in arg("--kk", "5,20,50").split(",")]

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
K = Y.shape[1]
_, g = np.unique(clone, return_inverse=True)
miss = ~Y
MISSF = miss.astype(np.float64)


def fit_twoway(M, iters=40):
    cm = np.clip(M.mean(0), 1e-3, 1 - 1e-3)
    be = np.log(cm / (1 - cm)); al = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :])))
            if ax == 0:
                gg = (M - p).sum(1); h = (p * (1 - p)).sum(1)
                al = np.clip(al + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
                al -= al.mean()
            else:
                gg = (M - p).sum(0); h = (p * (1 - p)).sum(0)
                be = np.clip(be + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :]))), 1e-6, 1 - 1e-6)


P = fit_twoway(MISSF); eta0 = np.log(P / (1 - P))
o0 = np.argsort(g, kind="stable"); b0 = np.concatenate([[0], np.cumsum(np.bincount(g))])
for a_, b_ in zip(b0[:-1], b0[1:]):
    ix = o0[a_:b_]
    e = eta0[ix]; tgt = MISSF[ix].sum(0); gam = np.zeros(K)
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)
ETA = np.log(P / (1 - P))
SDF = np.sqrt(0.01 * 0.99)
RP = (MISSF - P) / np.maximum(np.sqrt(P * (1 - P)), SDF)

rng = np.random.default_rng(SEED)
if POOLED:
    sizes = np.bincount(g)
    groups = []
    for c in range(sizes.size):
        if sizes[c] < MINCL:
            continue
        ix = np.flatnonzero(g == c)
        if ix.size > NSUB:
            ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
        groups.append(ix)
else:
    ix = np.flatnonzero(g == CL)
    if ix.size > NSUB:
        ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
    groups = [ix]
cells = np.concatenate(groups)
n = cells.size
off = np.concatenate([[0], np.cumsum([x.size for x in groups])])
codes = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"][cells]
ETAc, RPc, Mc = ETA[cells], RP[cells], MISSF[cells]
Pc = P[cells]
print(f"{arm} {'POOLED ' + str(len(groups)) + ' clones' if POOLED else 'clone ' + str(CL)}: "
      f"{n:,} cells, k in {KK}, {NSPLIT} tape splits", flush=True)


def fit_w(Xv, etav, sv, iters=60):
    """One-parameter logistic fit with a fixed offset."""
    w = 0.0
    for _ in range(iters):
        q = 1.0 / (1.0 + np.exp(-(etav + w * sv)))
        gr = np.dot(sv, Xv - q); he = np.dot(sv * sv, q * (1 - q))
        if he < 1e-12:
            break
        st = gr / he; w += float(np.clip(st, -2, 2))
        if abs(st) < 1e-10:
            break
    return w


def ll(Xv, etav, sv, w):
    q = np.clip(1.0 / (1.0 + np.exp(-(etav + w * sv))), 1e-12, 1 - 1e-12)
    return float(np.sum(Xv * np.log(q) + (1 - Xv) * np.log1p(-q)))


out = {"arm": arm, "pooled": POOLED, "clone": None if POOLED else CL,
       "n_cells": n, "nsplit": NSPLIT, "k_list": KK, "runs": []}
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    agree = np.zeros((n, n), np.int16); denom = np.zeros((n, n), np.int16)
    for a_ in A:
        v0 = codes[:, a_, 0] >= 0
        denom += v0[:, None] & v0[None, :]
        for d in range(6):
            cdv = codes[:, a_, d]; okd = cdv >= 0
            agree += (cdv[:, None] == cdv[None, :]) & okd[:, None] & okd[None, :]
    rel = np.where(denom > 0, agree / np.maximum(denom, 1), -1.0)
    np.fill_diagonal(rel, -np.inf)                       # never a neighbour of itself
    for b_ in range(len(groups)):                        # never across clones
        lo, hi = off[b_], off[b_ + 1]
        rel[lo:hi, :lo] = -np.inf; rel[lo:hi, hi:] = -np.inf
    tr = rs.random(n) < 0.5                              # train / held-out cells
    pmN = np.concatenate([lo + rs.permutation(hi - lo)
                          for lo, hi in zip(off[:-1], off[1:])])   # within-clone null
    for k in KK:
        nb = np.argsort(-rel, axis=1)[:, :k]
        for tag, RPu, Mu in (("obs", RPc, Mc), ("null", RPc[pmN], Mc[pmN])):
            u = RPu[nb][:, :, Bh].mean(1)                # [cells, |B|]
            f = Mu[nb][:, :, Bh].mean(1)
            X = Mc[:, Bh]; E = ETAc[:, Bh]
            wu = fit_w(X[tr].ravel(), E[tr].ravel(), u[tr].ravel())
            fc = f - f.mean()
            wf = fit_w(X[tr].ravel(), E[tr].ravel(), fc[tr].ravel())
            te = ~tr
            gain = ll(X[te].ravel(), E[te].ravel(), u[te].ravel(), wu) - \
                   ll(X[te].ravel(), E[te].ravel(), u[te].ravel(), 0.0)
            nats = gain / max(te.sum(), 1) * (K / Bh.size)
            out["runs"].append(dict(split=s, k=k, kind=tag, w_u=wu, w_f=wf,
                                    nats_per_cell=nats, n_test=int(te.sum())))
        print(f"  split {s} k={k}: nats/cell obs {out['runs'][-2]['nats_per_cell']:+.3f} "
              f"null {out['runs'][-1]['nats_per_cell']:+.3f}  "
              f"OR(all vs none) {np.exp(out['runs'][-2]['w_f']):.2f}  [{time.time()-t0:.0f}s]",
              flush=True)
    if s == 0:                                            # the digestible table, once
        k = KK[len(KK) // 2]
        nb = np.argsort(-rel, axis=1)[:, :k]
        f = Mc[nb][:, :, Bh].mean(1); X = Mc[:, Bh]; pv = Pc[:, Bh]
        FB = np.array([-.01, .001, .25, .5, .75, .999, 1.01])
        FL = ["none", "<25%", "25-50%", "50-75%", "75-<100%", "all"]
        qs = np.quantile(pv.ravel(), [0, .25, .5, .75, 1.0])
        print(f"\n  CONDITIONAL TABLE (k={k}, split 0): observed dropout rate, "
              f"holding the model's prediction fixed")
        print(f"  {'model says':>14}" + "".join(f"{l:>13}" for l in FL))
        tab = []
        for i in range(4):
            m0 = (pv >= qs[i]) & (pv <= qs[i + 1])
            row, cnt = [], []
            for j in range(len(FL)):
                m = m0 & (f > FB[j]) & (f <= FB[j + 1])
                row.append(float(X[m].mean()) if m.sum() else np.nan); cnt.append(int(m.sum()))
            tab.append(dict(p_lo=float(qs[i]), p_hi=float(qs[i+1]),
                            p_mid=float(np.median(pv[m0])), rates=row, counts=cnt))
            print(f"  {100*np.median(pv[m0]):>13.0f}%" +
                  "".join(f"{100*r:>12.0f}%" if np.isfinite(r) else f"{'-':>13}" for r in row))
            print(f"  {'n =':>14}" + "".join(f"{c:>13,}" for c in cnt))
        out["table"] = dict(k=k, f_labels=FL, rows=tab)
tag = f"{arm}_" + ("pooled" if POOLED else f"c{CL}")
(RES / f"profilepred_{tag}.json").write_text(json.dumps(out, indent=1))
R = out["runs"]
print(f"\n  {'k':>5}{'nats/cell obs':>15}{'null':>9}{'obs-null':>10}{'sd':>8}"
      f"{'OR obs':>9}{'OR null':>9}")
for k in KK:
    o = np.array([r["nats_per_cell"] for r in R if r["k"] == k and r["kind"] == "obs"])
    u = np.array([r["nats_per_cell"] for r in R if r["k"] == k and r["kind"] == "null"])
    fo = np.array([np.exp(r["w_f"]) for r in R if r["k"] == k and r["kind"] == "obs"])
    fu = np.array([np.exp(r["w_f"]) for r in R if r["k"] == k and r["kind"] == "null"])
    print(f"{k:>5}{o.mean():>+15.3f}{u.mean():>+9.3f}{(o-u).mean():>+10.3f}"
          f"{(o-u).std(ddof=1) if o.size>1 else 0:>8.3f}{fo.mean():>9.2f}{fu.mean():>9.2f}")
print(f"\nwrote results/profilepred_{tag}.json")
