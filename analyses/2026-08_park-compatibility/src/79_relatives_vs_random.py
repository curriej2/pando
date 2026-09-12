#!/usr/bin/env python3
r"""
DATA FOR FIG C -- "relatives against strangers", the panel that needs no statistics.

THE QUESTION, stripped of vocabulary: if I tell you what fraction of a cell's
NEAREST RELATIVES are missing a tape, does that predict whether the cell is missing
it?  And what if I use ARBITRARY cells from the same clone instead?

The second series IS the null the ladder uses (k random clone-mates, self excluded),
so the flat line in the figure is not a claim -- it is the null drawn from data.

CONSTRUCTION, identical to 77 so the figure and the ladder are the same analysis:
  * tape split A / B: relatedness from A only, dropout from B only
  * entry mask: fit on training entries, PLOT ONLY HELD-OUT ENTRIES
  * M3 = alpha_c + beta_z fitted inside each clone on training entries
  * f_rel = fraction of the k nearest relatives missing that tape
    f_rnd = fraction of k random clone-mates missing it, self excluded
    both averaged over the neighbours' TRAINING entries only
  * bin held-out entries by f, WITHIN quartiles of the model's own prediction p~,
    so cell quality and tape quality are held fixed by construction and cannot
    explain any contrast the figure shows

Usage: 79_relatives_vs_random.py <arm> [--mincl 100] [--k 20] [--nsplit 5]
       [--nmask 2] [--frac 0.5] [--seed 11]
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
FB = np.array([-.01, .001, .25, .5, .75, .999, 1.01])
FL = ["none", "<25%", "25-50%", "50-75%", "75-<100%", "all"]


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
MINCL, KNN = arg("--mincl", 100, int), arg("--k", 20, int)
NSPLIT, NMASK = arg("--nsplit", 5, int), arg("--nmask", 2, int)
FRAC, SEED = arg("--frac", 0.5, float), arg("--seed", 11, int)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
K = Y.shape[1]
_, g0 = np.unique(clone, return_inverse=True)
sizes = np.bincount(g0)
keep = np.isin(g0, np.flatnonzero(sizes >= MINCL))
Y, g0 = Y[keep], g0[keep]
_, g = np.unique(g0, return_inverse=True)
codes_all = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"]
assert codes_all.shape[0] == keep.size
codes_all = codes_all[keep]
X = (~Y).astype(np.float64)
n, nC = X.shape[0], int(g.max()) + 1
o = np.argsort(g, kind="stable"); X, g, codes_all = X[o], g[o], codes_all[o]
off = np.concatenate([[0], np.cumsum(np.bincount(g))])
print(f"{arm}: {n:,} cells, {nC} clones (>= {MINCL}), k = {KNN}, "
      f"{NSPLIT} splits x {NMASK} masks", flush=True)
sig = lambda e: 1.0 / (1.0 + np.exp(-e))


def fit_m3(Xb, Mb, iters=80):
    m, q = Xb.shape
    tw = Mb.sum(0)
    r0 = np.clip((Mb * Xb).sum(0) / np.maximum(tw, 1), 1e-3, 1 - 1e-3)
    be = np.log(r0 / (1 - r0)); al = np.zeros(m)
    for _ in range(iters):
        p = sig(al[:, None] + be[None, :])
        gr = (Mb * (Xb - p)).sum(1); he = (Mb * p * (1 - p)).sum(1)
        al = np.clip(al + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -12, 12)
        al -= al.mean()
        p = sig(al[:, None] + be[None, :])
        gr = (Mb * (Xb - p)).sum(0); he = (Mb * p * (1 - p)).sum(0)
        be = np.clip(be + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -12, 12)
    return sig(al[:, None] + be[None, :])


NQ = 4
cnt = {t: np.zeros((NQ, len(FL))) for t in ("rel", "rnd")}
hit = {t: np.zeros((NQ, len(FL))) for t in ("rel", "rnd")}
pmid = np.zeros(NQ); pn = 0
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    XB = X[:, Bh]; nb = Bh.size
    NBR = np.zeros((n, KNN), np.int64)
    for b_ in range(nC):
        lo, hi = off[b_], off[b_ + 1]; m = hi - lo
        cb = codes_all[lo:hi]
        agree = np.zeros((m, m), np.int16); denom = np.zeros((m, m), np.int16)
        for a_ in A:
            v0 = cb[:, a_, 0] >= 0
            denom += v0[:, None] & v0[None, :]
            for d in range(6):
                cdv = cb[:, a_, d]; ok = cdv >= 0
                agree += (cdv[:, None] == cdv[None, :]) & ok[:, None] & ok[None, :]
        rel = np.where(denom > 0, agree / np.maximum(denom, 1), -1.0)
        np.fill_diagonal(rel, -np.inf)
        keff = min(KNN, m - 1)
        nbb = np.argpartition(-rel, keff - 1, axis=1)[:, :keff]
        if keff < KNN:
            nbb = nbb[:, np.arange(KNN) % keff]
        NBR[lo:hi] = nbb + lo
        del rel, agree, denom
    for r in range(NMASK):
        rm = np.random.default_rng([SEED, s, r, 99])
        M = rm.random((n, nb)) < FRAC
        Tr = M.astype(np.float64); Va = ~M
        P = np.zeros((n, nb))
        for b_ in range(nC):
            lo, hi = off[b_], off[b_ + 1]
            P[lo:hi] = np.clip(fit_m3(XB[lo:hi], Tr[lo:hi]), 1e-6, 1 - 1e-6)
        rn = np.random.default_rng([SEED, s, r, 1234])
        RND = np.zeros_like(NBR)
        for b_ in range(nC):
            lo, hi = off[b_], off[b_ + 1]; m = hi - lo
            loc = np.arange(m)
            dr = rn.integers(0, m - 1, size=(m, KNN))
            dr += (dr >= loc[:, None])
            RND[lo:hi] = dr + lo
        qs = np.quantile(P[Va], np.linspace(0, 1, NQ + 1))
        for i in range(NQ):
            m0 = Va & (P >= qs[i]) & (P <= qs[i + 1])
            pmid[i] += np.median(P[m0]); 
            for tag, IDX in (("rel", NBR), ("rnd", RND)):
                mz = Tr[IDX].sum(1)
                f = (Tr[IDX] * XB[IDX]).sum(1) / np.maximum(mz, 1.0)
                f[mz == 0] = np.nan
                for j in range(len(FL)):
                    mm = m0 & (f > FB[j]) & (f <= FB[j + 1])
                    c = int(mm.sum())
                    cnt[tag][i, j] += c
                    hit[tag][i, j] += float(XB[mm].sum())
        pn += 1
        print(f"  split {s} mask {r} done [{time.time()-t0:.0f}s]", flush=True)
pmid /= pn
out = dict(arm=arm, mincl=MINCL, k=KNN, n_cells=int(n), n_clones=int(nC),
           nsplit=NSPLIT, nmask=NMASK, f_labels=FL, p_mid=[float(v) for v in pmid],
           counts={t: cnt[t].tolist() for t in cnt}, hits={t: hit[t].tolist() for t in hit})
(RES / f"relvsrnd_{arm}_mincl{MINCL}_k{KNN}.json").write_text(json.dumps(out, indent=1))
print(f"\n  observed dropout rate by neighbour signal, held-out entries only")
for i in range(NQ):
    print(f"\n  model says {100*pmid[i]:.0f}%")
    for tag, nm in (("rel", "nearest relatives"), ("rnd", "random clone-mates")):
        rr = np.where(cnt[tag][i] > 0, hit[tag][i] / np.maximum(cnt[tag][i], 1), np.nan)
        print(f"    {nm:>19}" + "".join(f"{(f'{100*v:.0f}%' if np.isfinite(v) else '-'):>11}" for v in rr))
        print(f"    {'n =':>19}" + "".join(f"{int(c):>11,}" for c in cnt[tag][i]))
print(f"\nwrote results/relvsrnd_{arm}_mincl{MINCL}_k{KNN}.json")
