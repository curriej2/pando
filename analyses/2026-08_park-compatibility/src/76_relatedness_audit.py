#!/usr/bin/env python3
r"""
Is "lineage relatedness" secretly a measure of how many tapes a pair SHARES?

JUSTIN'S WORRY (2026-09-11), and it is the right one to have.  Relatedness is
computed on half A, dropout on half B -- but if relatedness were larger for pairs
that jointly recover more A-tapes, then relatedness would be partly a *capture*
statistic.  Capture is exactly what dropout is.  The disjoint tape split stops A's
dropout from being read directly as B's dropout, but it does NOT by itself stop a
pair's overall capture quality from leaking across the split.  So: measure it.

THE METRIC, as actually computed in 73/74 (verified against the builder):
    denom_{cc'} = # A-tapes DETERMINED IN BOTH cells at depth 1
    agree_{cc'} = # (tape, depth) cells where both codes are valid AND equal
    rel_{cc'}   = agree / denom  =  MEAN shared prefix depth per jointly
                  determined A-tape,  in [0, 6]
Validity is nested in depth (44_prefix_codes6.py: `ok &= S[:,:,d] >= 0`), so
`agree` can only count tapes already counted in `denom`.  ⇒ a pair sharing just
two tapes at shared depths 4 and 5 scores 9/2 = 4.5, the same as a pair sharing
eighty tapes at a mean depth of 4.5.  The metric is a MEAN, not a total.
⚠ "determined" is stricter than "recovered": it needs a kept-alphabet symbol at
site 1, so a tape read as junk at site 1 is absent from denom though the cell
counted it towards the >=20/>=100 tape filter.

WHAT THIS SCRIPT REPORTS
 A. brute force -- recompute rel for sampled pairs with an independent loop.
 B. the distribution of denom, and corr(rel, denom).  A POSITIVE correlation here
    is expected and is not by itself a problem; it is what C and D test.
 C. mean denom per relatedness bin -- is the top bin the well-captured pairs?
 D. ⭐ THE CONTROL THAT SETTLES IT: the variogram computed WITHIN narrow strata of
    denom.  If the relatedness gradient survives at fixed joint-capture, then the
    gradient is not a capture artefact.  If it flattens, it was.

Usage: 76_relatedness_audit.py <arm> [--pooled] [--clone C] [--nsub 0]
       [--nsplit 4] [--sdfloor 0.01] [--seed S]
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
NSUB, NSPLIT = arg("--nsub", 0, int), arg("--nsplit", 4, int)
MINCL, SEED = arg("--mincl", 20, int), arg("--seed", 11, int)
PFLOOR = arg("--sdfloor", 0.01, float)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
K = Y.shape[1]
_, g = np.unique(clone, return_inverse=True)
MISSF = (~Y).astype(np.float64)


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


P = fit_twoway(MISSF); eta = np.log(P / (1 - P))
o0 = np.argsort(g, kind="stable"); b0 = np.concatenate([[0], np.cumsum(np.bincount(g))])
for a_, b_ in zip(b0[:-1], b0[1:]):
    ix = o0[a_:b_]
    e = eta[ix]; tgt = MISSF[ix].sum(0); gam = np.zeros(K)
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)
SDF = np.sqrt(PFLOOR * (1 - PFLOOR))
RP = (MISSF - P) / np.maximum(np.sqrt(P * (1 - P)), SDF)

rng = np.random.default_rng(SEED)
sizes = np.bincount(g)
if POOLED:
    groups = []
    for c in range(sizes.size):
        if sizes[c] < MINCL:
            continue
        ix = np.flatnonzero(g == c)
        if NSUB and ix.size > NSUB:
            ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
        groups.append(ix)
else:
    ix = np.flatnonzero(g == CL)
    if NSUB and ix.size > NSUB:
        ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
    groups = [ix]
cells = np.concatenate(groups); n = cells.size
off = np.concatenate([[0], np.cumsum([x.size for x in groups])])
codes = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"][cells]
RPc = RP[cells]
print(f"{arm} {'POOLED ' + str(len(groups)) + ' clones' if POOLED else 'clone ' + str(CL)}: "
      f"{n:,} cells, {NSPLIT} tape splits", flush=True)

RR, DD, SS, NN = [], [], [], []
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    RPB = RPc[:, Bh]
    pm = rs.permutation(n); RPN = RPB[pm]
    for b_ in range(len(groups)):
        lo, hi = off[b_], off[b_ + 1]; m = hi - lo
        if m < 2:
            continue
        cb = codes[lo:hi]
        agree = np.zeros((m, m), np.int16); denom = np.zeros((m, m), np.int16)
        for a_ in A:
            v0 = cb[:, a_, 0] >= 0
            denom += v0[:, None] & v0[None, :]
            for d in range(6):
                cdv = cb[:, a_, d]; ok = cdv >= 0
                agree += (cdv[:, None] == cdv[None, :]) & ok[:, None] & ok[None, :]
        tu = np.triu_indices(m, 1)
        dn = denom[tu].astype(np.int16); ag = agree[tu]
        if s == 0 and b_ == 0 and m >= 4:
            # ---- A. BRUTE FORCE, independent of the vectorised path -----------
            print("\n  A. brute-force check on 5 pairs of the first clone "
                  "(independent reimplementation):", flush=True)
            rs2 = np.random.default_rng(7)
            for _ in range(5):
                i, j = rs2.choice(m, 2, replace=False)
                shared, ntap = [], 0
                for a_ in A:
                    if cb[i, a_, 0] < 0 or cb[j, a_, 0] < 0:
                        continue
                    ntap += 1
                    d_ = 0
                    for dd_ in range(6):
                        if cb[i, a_, dd_] >= 0 and cb[i, a_, dd_] == cb[j, a_, dd_]:
                            d_ += 1
                        else:
                            break
                    shared.append(d_)
                bf = sum(shared) / ntap if ntap else np.nan
                vec = agree[i, j] / denom[i, j] if denom[i, j] else np.nan
                print(f"     jointly determined tapes = {ntap:3d}  "
                      f"depths {sorted(shared, reverse=True)[:6]}{'...' if ntap > 6 else ''}  "
                      f"brute force {bf:.6f}  vectorised {vec:.6f}  "
                      f"{'OK' if abs(bf - vec) < 1e-9 else '*** MISMATCH ***'}", flush=True)
        keep = dn > 0
        RR.append((ag[keep] / dn[keep]).astype(np.float32))
        DD.append(dn[keep])
        Xb, Nb = RPB[lo:hi], RPN[lo:hi]
        SS.append(((Xb @ Xb.T) / Bh.size)[tu][keep].astype(np.float32))
        NN.append(((Nb @ Nb.T) / Bh.size)[tu][keep].astype(np.float32))
    print(f"  split {s+1}/{NSPLIT} done [{time.time()-t0:.0f}s]", flush=True)

rel = np.concatenate(RR); dnm = np.concatenate(DD).astype(np.float64)
sim = np.concatenate(SS); nul = np.concatenate(NN)
del RR, DD, SS, NN
print(f"\n  {rel.size:,} pair-observations ({NSPLIT} splits x within-clone pairs)")

# ---- B. does relatedness track joint capture? --------------------------------
r_pear = np.corrcoef(rel, dnm)[0, 1]
rr = np.argsort(np.argsort(rel)).astype(np.float64)
dr = np.argsort(np.argsort(dnm)).astype(np.float64)
r_spear = np.corrcoef(rr, dr)[0, 1]
del rr, dr
print(f"\n  B. joint capture per pair (# A-tapes determined in BOTH cells), |A|={len(A)}")
print(f"     denom: min {dnm.min():.0f}  q1 {np.quantile(dnm,.25):.0f}  median "
      f"{np.median(dnm):.0f}  q3 {np.quantile(dnm,.75):.0f}  max {dnm.max():.0f}")
print(f"     corr(rel, denom)  Pearson {r_pear:+.4f}   Spearman {r_spear:+.4f}")

# ---- C. is the top relatedness bin just the well-captured pairs? -------------
QS = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .98, .995, 1.0])
E = np.unique(np.quantile(rel, QS)); E[0] -= 1e-6; E[-1] += 1e-6
bi = np.clip(np.searchsorted(E, rel, side="right") - 1, 0, len(E) - 2)
nb = len(E) - 1
cnt = np.bincount(bi, minlength=nb).astype(float)
mean = lambda v: np.bincount(bi, weights=v, minlength=nb) / np.maximum(cnt, 1)
print(f"\n  C. {'rel bin':>16}{'mean rel':>10}{'mean denom':>12}{'pairs':>14}{'obs-null':>11}")
mr, md, ms, mn = mean(rel.astype(np.float64)), mean(dnm), mean(sim.astype(np.float64)), mean(nul.astype(np.float64))
for i in range(nb):
    print(f"     {E[i]:>7.3f}-{E[i+1]:<8.3f}{mr[i]:>10.3f}{md[i]:>12.1f}{cnt[i]:>14,.0f}"
          f"{ms[i]-mn[i]:>+11.4f}")

# ---- D. ⭐ THE CONTROL: the gradient WITHIN strata of joint capture ----------
print(f"\n  D. ⭐ variogram WITHIN strata of joint capture (denom).  If the "
      f"relatedness\n     gradient survives at fixed denom, it is not a capture artefact.")
dq = np.quantile(dnm, [0, .25, .5, .75, 1.0])
NRB = 6
for q in range(4):
    m0 = (dnm >= dq[q]) & (dnm <= dq[q + 1]) if q == 3 else (dnm >= dq[q]) & (dnm < dq[q + 1])
    if m0.sum() < 1000:
        continue
    rq = rel[m0]
    eb = np.unique(np.quantile(rq, np.linspace(0, 1, NRB + 1)))
    eb[0] -= 1e-6; eb[-1] += 1e-6
    b2 = np.clip(np.searchsorted(eb, rq, side="right") - 1, 0, len(eb) - 2)
    c2 = np.bincount(b2, minlength=len(eb) - 1).astype(float)
    f = lambda v: np.bincount(b2, weights=v[m0].astype(np.float64), minlength=len(eb) - 1) / np.maximum(c2, 1)
    d2 = f(sim) - f(nul); r2 = np.bincount(b2, weights=rq.astype(np.float64),
                                           minlength=len(eb) - 1) / np.maximum(c2, 1)
    print(f"\n     denom stratum {q+1}/4: {dq[q]:.0f}-{dq[q+1]:.0f} tapes, "
          f"{int(m0.sum()):,} pairs, mean denom {dnm[m0].mean():.1f}")
    print(f"       {'rel':>8}" + "".join(f"{v:>9.2f}" for v in r2))
    print(f"       {'obs-null':>8}" + "".join(f"{v:>+9.4f}" for v in d2))
    print(f"       {'pairs':>8}" + "".join(f"{int(v):>9,}" for v in c2))
    print(f"       ⇒ lowest {d2[0]:+.4f}  highest {d2[-1]:+.4f}  "
          f"spread {d2[-1]-d2[0]:+.4f}")
