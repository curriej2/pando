#!/usr/bin/env python3
r"""
Is dropout structured by LINEAGE, continuously?  A variogram of dropout residuals
against lineage relatedness -- no clades, no thresholds, no event definition.

⚠⚠ THE CONFOUND THAT WOULD OTHERWISE KILL IT.  Lineage relatedness is read from the
edit data, and dropout decides which edits are readable.  Two cells that both lack
tapes 1-50 "agree" at those positions, so a naive distance would call them related
BECAUSE their dropout matches -- predicting the answer from itself.
⇒ DISJOINT TAPE SPLIT.  Partition the 166 tapes into halves A and B.  Lineage
  relatedness uses ONLY half A; dropout similarity uses ONLY half B.  No shared
  data, so no circularity.  Repeated over R random splits for power and a spread.

TWO FURTHER CONTROLS
 * residuals, not raw missingness.  r = X - p~ with p~ = sigma(alpha_c + beta_z +
   gamma_{C,z}), so the per-cell capture margin (rho_cell = 0.13), the per-tape
   margin (rho_tape = 0.25) and the clone's own per-tape rate are already removed.
   What is left is sub-clonal.
 * within one clone.  Clone identity cannot do the work.

DEFINITIONS, per cell pair (c, c') inside the clone
  relatedness = (sum over a in A, d in 1..6 of 1[prefix codes agree at depth d])
                / (number of a in A both cells recovered)
              = mean shared prefix depth per jointly recovered A-tape, in [0, 6].
    Prefixes are nested, so counting matching depths IS the shared depth.
  similarity   = (1/|B|) * sum over z in B of r_cz * r_c'z        (a covariance)
    This is the numerator of the rho_within statistic from the row-A9 work, so the
    curve is a continuous generalisation of the rho-by-subclade-depth gradient.

NULL.  Permute the cell order of the half-B residual matrix only.  That destroys
the pairing between lineage (from A) and dropout (from B) while preserving every
marginal, and is binned against the SAME relatedness values.

⚠ E[similarity] over all within-clone pairs is slightly NEGATIVE, not zero: the
gamma fit forces sum_c r_cz = 0, so sum_{c != c'} r_c r_c' = -sum_c r_c^2.  The
reference level is therefore the all-pairs mean, which is printed.

Usage: 73_dropout_variogram.py <arm> --clone C [--nsub 2000] [--nsplit 10] [--seed S]
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
NSUB = arg("--nsub", 2000, int)
NSPLIT = arg("--nsplit", 10, int)
SEED = arg("--seed", 11, int)

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
R = MISSF - P
# ⚠ PEARSON residuals: each residual divided by its own null SD, sqrt(p~(1-p~)),
# so every term has mean 0 and variance 1 under H0 regardless of where the cell and
# tape sit on the p~ range.  Without this the mean product is scaled by p~(1-p~),
# which varies ~250-fold across (cell, tape), and any covariation between
# relatedness and p~ level would move the curve for non-lineage reasons.
# ⚠ The SD is floored so a residual cannot exceed ~10: p~ is clipped to [1e-6,
# 1-1e-6], and an unfloored Pearson residual could reach 1000 for a cell missing a
# tape the model says is essentially always present.  Those terms are the most
# informative ones but would single-handedly set the mean.
SDFLOOR = np.sqrt(0.01 * 0.99)
SD = np.maximum(np.sqrt(P * (1 - P)), SDFLOOR)
RP = R / SD
print(f"  Pearson residuals: SD floored at {SDFLOOR:.4f} (max |r*| = "
      f"{1/SDFLOOR:.1f}); floor binds on "
      f"{100*np.mean(np.sqrt(P*(1-P)) < SDFLOOR):.2f}% of (cell, tape) entries")

POOLED = "--pooled" in sys.argv
MINCL = arg("--mincl", 20, int)
rng = np.random.default_rng(SEED)
if POOLED:
    sizes = np.bincount(g)
    use = [c for c in range(sizes.size) if sizes[c] >= MINCL]
    groups = []
    for c in use:
        ix = np.flatnonzero(g == c)
        if ix.size > NSUB:
            ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
        groups.append(ix)
    tot_pairs = sum(x.size * (x.size - 1) // 2 for x in groups)
    allp = int((sizes * (sizes - 1) // 2).sum())
    print(f"{arm} POOLED: {len(groups):,} clones with >= {MINCL} cells "
          f"({sum(x.size for x in groups):,} cells); {tot_pairs:,} within-clone pairs "
          f"= {100*tot_pairs/max(allp,1):.0f}% of all {allp:,}", flush=True)
else:
    ix = np.flatnonzero(g == CL)
    if ix.size > NSUB:
        ix = np.sort(rng.choice(ix, size=NSUB, replace=False))
    groups = [ix]
    print(f"{arm} clone {CL}: {int((g == CL).sum()):,} cells, using {ix.size:,}", flush=True)
cells = np.concatenate(groups)
n = cells.size
codes = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"][cells]
Rc, RPc = R[cells], RP[cells]
off = np.concatenate([[0], np.cumsum([x.size for x in groups])])
print(f"{arm} clone {CL}: {int((g == CL).sum()):,} cells, using {n:,}; "
      f"{n*(n-1)//2:,} pairs; {NSPLIT} tape splits", flush=True)
print(f"  gamma check: max |sum_c r_cz| over tapes in this clone = "
      f"{np.abs(R[g == CL].sum(0)).max():.2e}", flush=True)

# ⚠ deciles up to 0.9 then finer, because the top decile's MEAN relatedness (4.42
# on Mouse2 c76) is well short of the top of the range (~5.5): the closest pairs
# were being averaged away with merely-close ones.
QS = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .98, .995, 1.0])
EDGES = None                      # set from split 0's relatedness, then reused
# pair indices, restricted to WITHIN-clone pairs
ii, jj = [], []
for b_ in range(len(groups)):
    lo, hi = off[b_], off[b_ + 1]
    a_, c_ = np.triu_indices(hi - lo, 1)
    ii.append(a_ + lo); jj.append(c_ + lo)
iu = (np.concatenate(ii), np.concatenate(jj))
print(f"  {iu[0].size:,} within-clone pairs extracted", flush=True)
curves, nulls, refs, relmeans, cnts = [], [], [], [], []
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    agree = np.zeros((n, n), np.int16)
    denom = np.zeros((n, n), np.int16)
    for a_ in A:
        v0 = codes[:, a_, 0] >= 0
        both = v0[:, None] & v0[None, :]
        denom += both
        for d in range(6):
            cdv = codes[:, a_, d]
            ok = cdv >= 0
            agree += ((cdv[:, None] == cdv[None, :]) & ok[:, None] & ok[None, :])
    rel = np.where(denom > 0, agree / np.maximum(denom, 1), np.nan)
    RB, RPB = Rc[:, Bh], RPc[:, Bh]
    sim = (RPB @ RPB.T) / Bh.size            # PEARSON -- the primary statistic
    cov = (RB @ RB.T) / Bh.size              # raw covariance, kept for comparability
    pm = rs.permutation(n)
    simN = (RPB[pm] @ RPB[pm].T) / Bh.size
    r_, s_, sn_, c_ = rel[iu], sim[iu], simN[iu], cov[iu]
    keep = np.isfinite(r_)
    r_, s_, sn_, c_ = r_[keep], s_[keep], sn_[keep], c_[keep]
    if EDGES is None:                        # quantile bins, from split 0, reused
        EDGES = np.unique(np.quantile(r_, QS))
        EDGES[0] -= 1e-9; EDGES[-1] += 1e-9
    bi = np.clip(np.searchsorted(EDGES, r_, side="right") - 1, 0, len(EDGES) - 2)
    nb = len(EDGES) - 1
    cnt = np.bincount(bi, minlength=nb).astype(float)
    mean = lambda v: np.where(cnt > 0, np.bincount(bi, weights=v, minlength=nb) / np.maximum(cnt, 1), np.nan)
    curves.append(mean(s_)); nulls.append(mean(sn_)); relmeans.append(mean(r_))
    cnts.append(cnt); refs.append((s_.mean(), c_.mean()))
    print(f"  split {s+1}/{NSPLIT}: |A|={A.size} |B|={Bh.size}  all-pairs Pearson "
          f"{s_.mean():+.4f}  raw cov {c_.mean():+.2e}  [{time.time()-t0:.0f}s]", flush=True)

C = np.array(curves); N = np.array(nulls); RM = np.array(relmeans); CNT = np.array(cnts).mean(0)
allref = float(np.mean([r[0] for r in refs])); allcov = float(np.mean([r[1] for r in refs]))
# ⚠ the zero-sum constraint sum_c r_cz = 0 applies to the RAW residuals, not to
# r/SD, because SD varies by cell -- so the Pearson all-pairs mean is NOT forced
# negative.  It is the empirical reference level, and the curve crosses it, not 0.
print(f"\n  reference: all-pairs mean Pearson similarity = {allref:+.5f}"
      f"  (the level the curve crosses; not forced to 0)")
print(f"             all-pairs raw covariance          = {allcov:+.3e}")
D = C - N
print(f"\n  {'relatedness bin':>17}{'mean rel':>10}{'pairs':>12}{'observed':>11}"
      f"{'null':>11}{'obs-null':>11}{'sd/splits':>11}{'t':>7}")
for i in range(C.shape[1]):
    sd = np.nanstd(D[:, i], ddof=1); mu = np.nanmean(D[:, i])
    t = mu / (sd / np.sqrt(np.isfinite(D[:, i]).sum())) if sd > 0 else np.nan
    print(f"  {EDGES[i]:>7.3f}-{EDGES[i+1]:<8.3f}{np.nanmean(RM[:, i]):>10.3f}{CNT[i]:>12,.0f}"
          f"{np.nanmean(C[:, i]):>+11.4f}{np.nanmean(N[:, i]):>+11.4f}{mu:>+11.4f}"
          f"{sd:>11.4f}{t:>7.1f}")
TAG = f"{arm}_" + ("pooled" if POOLED else f"c{CL}")
np.savez_compressed(RES / f"variogram_{TAG}.npz", edges=EDGES, obs=C, null=N,
                    relmeans=RM, counts=CNT, allref=allref, allcov=allcov,
                    n_cells=n, nsplit=NSPLIT)
print(f"\nwrote results/variogram_{TAG}.npz")
