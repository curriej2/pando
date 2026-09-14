#!/usr/bin/env python3
r"""
PER-TAPE UNANIMITY -- Justin's design, 2026-09-12.

THE STATEMENT THE FIGURE MAKES, for ONE tape at a time:
  "Take k cells that are close relatives.  Take k cells at random.  How often is
   this tape missing in ALL k of them?"

WHY PER TAPE IS NECESSARY, not merely nicer.  Subclone's overall missing rate is
21.8%, so under independence a random set of 20 would all lack a tape with
probability 0.22^20 = 5e-14.  The pooled figure showed 4.3%.  That gap of twelve
orders of magnitude is entirely tapes that are ALREADY NEARLY DEAD inside a clone.
⇒ the pooled unanimity statistic measures how many tapes are nearly dead, not how
much relatedness clusters dropout.  Splitting by tape is the only way to see the
effect itself.

THREE SET TYPES, and the middle one is the real control:
  rel  the anchor plus its k-1 nearest relatives (relatedness from tape half A)
  mat  ⭐ CAPTURE-MATCHED random: every cell of the related set is replaced by a
       random cell FROM THE SAME CAPTURE DECILE of the same clone.  The matched set
       therefore has the same capture composition BY CONSTRUCTION -- not merely the
       same average -- and lineage is the only thing that differs.
  rnd  plain random from the clone, kept so we can see whether matching mattered

⚠⚠ WHY MARGINAL MATCHING WOULD NOT BE ENOUGH.  What drives all-k-missing is
within-set HOMOGENEITY of capture, not the mean.  Twenty uniformly badly-captured
cells reach unanimity far more easily than a mixed set with the same average, so
two groups can have identical capture marginals and still differ for a purely
technical reason.  Matching cell-by-cell on decile removes that.

⚠ Capture is measured on the A tapes only.  Using all 166 would correlate the
matching variable with the outcome being measured on B.
⚠⚠ MATCHING IS BY RANK CALIPER, not deciles.  Deciles left a residual gap -- on
Mouse2 the related sets averaged 37.8 A-tapes recovered against 40.9 for the
decile-matched control, because within a decile the related cells sit at the low
end.  A better-captured control produces FEWER all-missing sets, so the fold change
came out INFLATED.  Each cell is now replaced by one drawn within +-WIN of its own
capture rank, WIN = max(3, n_C/40).  The printed capture check is the validation.

⚠⚠ ALL-PRESENT IS COMPUTED TOO, AND IT IS A DISCRIMINATOR.  Shared capture quality
should cluster BOTH unanimity types -- uniformly good cells agree on presence,
uniformly bad ones agree on absence.  One-directional silencing should cluster
MISSING much more than PRESENT, because presence is just the default state.  So the
asymmetry between the two argues against the technical explanation from inside this
same figure.

⚠ P(all k missing) ~ q^k, so k = 20 needs q > 0.70 before the random-set
probability even reaches 1e-3, and most tapes would have a ZERO denominator.  Hence
the sweep over k.

⚠ Sets anchored on every cell would overlap heavily and be strongly correlated.
Anchors are SAMPLED, and no confidence interval is put on the result.

Usage: 83_unanimity.py <arm> [--mincl 100] [--kk 3,5,10,20] [--nsets 2000]
       [--nsplit 5] [--ndec 10] [--seed 11]
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
MINCL = arg("--mincl", 100, int)
KK = [int(v) for v in arg("--kk", "3,5,10,20").split(",")]
NSETS, NSPLIT = arg("--nsets", 2000, int), arg("--nsplit", 5, int)
NDEC, SEED = arg("--ndec", 10, int), arg("--seed", 11, int)

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
X = (~Y).astype(bool)                            # True = missing
n, nC = X.shape[0], int(g.max()) + 1
o = np.argsort(g, kind="stable"); X, g, codes_all = X[o], g[o], codes_all[o]
off = np.concatenate([[0], np.cumsum(np.bincount(g))]); csz = np.diff(off)
print(f"{arm}: {n:,} cells, {nC} clones (>= {MINCL}), {K} tapes; k in {KK}; "
      f"{NSETS} anchors/clone/split; {NSPLIT} splits", flush=True)

SER = ("rel", "mat", "rnd")
miss = {k: {s: np.zeros(K, np.int64) for s in SER} for k in KK}
pres = {k: {s: np.zeros(K, np.int64) for s in SER} for k in KK}
nset = {k: np.zeros(K, np.int64) for k in KK}
capsum = {k: {s: 0.0 for s in SER} for k in KK}
capsd = {k: {s: 0.0 for s in SER} for k in KK}
capn = {k: 0 for k in KK}
t0 = time.time()
for s_ in range(NSPLIT):
    rs = np.random.default_rng([SEED, s_])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    RA = (~X[:, A]).sum(1)                       # capture measured on A tapes only
    for b_ in range(nC):
        lo, hi = off[b_], off[b_ + 1]; m = hi - lo
        if m < max(KK) + 1:
            continue
        cb = codes_all[lo:hi]
        agree = np.zeros((m, m), np.int16); denom = np.zeros((m, m), np.int16)
        for a_ in A:
            v0 = cb[:, a_, 0] >= 0
            denom += v0[:, None] & v0[None, :]
            for d in range(6):
                cdv = cb[:, a_, d]; ok = cdv >= 0
                agree += (cdv[:, None] == cdv[None, :]) & ok[:, None] & ok[None, :]
        rel = np.where(denom > 0, agree / np.maximum(denom, 1), -1.0)
        del agree, denom
        # ⭐ CAPTURE MATCHING BY RANK CALIPER, from A tapes only.
        # Each cell is replaced by a random cell whose capture RANK inside the clone
        # is within +-WIN of its own, so the matched set reproduces the related set's
        # capture profile cell by cell rather than only on average.
        ra = RA[lo:hi]
        ordr = np.argsort(ra, kind="stable")
        rank = np.empty(m, np.int64); rank[ordr] = np.arange(m)
        WIN = max(3, m // 40)
        XBc = X[lo:hi][:, Bh]                    # dropout on B tapes only
        S = min(NSETS, m)
        anch = rs.choice(m, size=S, replace=False)
        for k in KK:
            # related set: the anchor together with its k-1 nearest
            r2 = rel[anch].copy()
            r2[np.arange(S), anch] = np.inf      # anchor is a member of its own set
            idx_rel = np.argpartition(-r2, k - 1, axis=1)[:, :k]
            # capture-matched random: same decile, cell by cell
            rk = rank[idx_rel]                   # (S, k) capture ranks
            jit = rs.integers(-WIN, WIN + 1, size=rk.shape)
            idx_mat = ordr[np.clip(rk + jit, 0, m - 1)]
            idx_rnd = rs.integers(0, m, size=(S, k))
            nset[k] += S * np.isin(np.arange(K), Bh)
            for tag, IDX in (("rel", idx_rel), ("mat", idx_mat), ("rnd", idx_rnd)):
                blk = XBc[IDX]                   # (S, k, |B|) bool
                miss[k][tag][Bh] += blk.all(1).sum(0)
                pres[k][tag][Bh] += (~blk).all(1).sum(0)
                cap = RA[lo:hi][IDX]             # capture profile of the set
                capsum[k][tag] += float(cap.mean()) * S
                capsd[k][tag] += float(cap.std(1).mean()) * S
            capn[k] += S
        del rel
    print(f"  split {s_+1}/{NSPLIT} [{time.time()-t0:.0f}s]", flush=True)

out = dict(arm=arm, mincl=MINCL, kk=KK, nsets=NSETS, nsplit=NSPLIT, ndec=NDEC,
           n_cells=int(n), n_clones=int(nC), K=int(K),
           nset={str(k): nset[k].tolist() for k in KK},
           miss={str(k): {s: miss[k][s].tolist() for s in SER} for k in KK},
           pres={str(k): {s: pres[k][s].tolist() for s in SER} for k in KK},
           cap_mean={str(k): {s: capsum[k][s] / capn[k] for s in SER} for k in KK},
           cap_within_sd={str(k): {s: capsd[k][s] / capn[k] for s in SER} for k in KK})
(RES / f"unanimity_{arm}_mincl{MINCL}.json").write_text(json.dumps(out))
print(f"\n  CAPTURE CHECK (mean A-tapes recovered per set, and mean WITHIN-set sd)")
for k in KK:
    print(f"   k={k:>3}: " + "  ".join(
        f"{s} {capsum[k][s]/capn[k]:.1f}/{capsd[k][s]/capn[k]:.1f}" for s in SER))
print(f"\n  TAPES ON SCALE (matched-random all-k-missing > 0) and median fold change")
for k in KK:
    ok = (nset[k] > 0) & (miss[k]["mat"] > 0)
    pr = miss[k]["rel"][ok] / nset[k][ok]; pm = miss[k]["mat"][ok] / nset[k][ok]
    fold = pr / pm
    up = int((pr > pm).sum())
    print(f"   k={k:>3}: {int(ok.sum()):>3}/{K} tapes on scale | "
          f"median fold {np.median(fold):>6.2f} | above the line {up}/{int(ok.sum())} | "
          f"pooled P_rel {miss[k]['rel'].sum()/max(nset[k].sum(),1):.4f} "
          f"P_mat {miss[k]['mat'].sum()/max(nset[k].sum(),1):.4f}")
print(f"\n  ALL-PRESENT (the discriminator: shared capture should lift BOTH)")
for k in KK:
    ok = (nset[k] > 0) & (pres[k]["mat"] > 0)
    fp = (pres[k]["rel"][ok] / nset[k][ok]) / (pres[k]["mat"][ok] / nset[k][ok])
    okm = (nset[k] > 0) & (miss[k]["mat"] > 0)
    fm = (miss[k]["rel"][okm] / nset[k][okm]) / (miss[k]["mat"][okm] / nset[k][okm])
    print(f"   k={k:>3}: median fold  MISSING {np.median(fm):>6.2f}   "
          f"PRESENT {np.median(fp):>6.2f}   ratio {np.median(fm)/np.median(fp):>5.2f}")
print(f"\nwrote results/unanimity_{arm}_mincl{MINCL}.json")
