#!/usr/bin/env python3
r"""
FIG C, MODEL-FREE.  No fitted model, no residuals, no stratification, no mask.

THE SENTENCE THE FIGURE MAKES
  "If j of your 20 closest relatives are missing this tape, you are missing it X%
   of the time.  If j of 20 arbitrary cells from your clone are, you are missing it
   Y% of the time."
Nothing is estimated anywhere in that statement.

WHY MODEL-FREE IS BETTER HERE (Justin, 2026-09-12).  The previous version binned
held-out entries by the fitted p~ of M3.  Two faults, both fatal to the caption:
  * p~ is an ESTIMATE, so stratifying on it does not hold the truth fixed.  A
    second 20-cell sample adds real information about the true clone-tape rate even
    when no lineage structure exists.  Verified by simulation: in a world with NO
    lineage structure the "arbitrary clone-mates" curve still rose 14.9% -> 28.2%.
    ⇒ the old caption implied the grey line should be flat.  It cannot be.
  * with a 50% training mask the denominator m ~ Binomial(20, 0.5) VARIES, so
    "all your neighbours lack it" meant "all m observed ones", and reaching that
    with m = 4 is ~15,600x easier by chance than with m = 10.

BOTH FAULTS GO AWAY HERE.  There is no estimate to condition on, and with no mask
the denominator is exactly 20, so x is a genuine count.

⚠ THE MASK IS DROPPED, AND THAT IS LEGITIMATE.  The mask exists so the LADDER's
likelihood comparison is honest.  This figure fits nothing and makes no likelihood
claim, so there is no overfitting to protect against.  Self-exclusion stays -- it
is structural, not a convention.

⚠⚠ THE A/B TAPE SPLIT STAYS, AND MUST.  Relatedness is read from half A and dropout
from half B.  Without it, two cells that both lack tapes 1-50 would look related
BECAUSE their dropout matches, and the figure would be predicting the answer from
itself.

⚠ WHAT THIS FIGURE DOES NOT CONTROL, stated rather than buried: alpha_c is no
longer removed, so if nearest relatives systematically share capture quality with
the cell, j_rel could predict for a technical reason.  corr(relatedness, joint
capture) runs -0.108 (Subclone) to +0.278 (Mouse3), so the association is real and
inconsistent in sign.  It is excluded ELSEWHERE -- `76` shows the gradient survives
in 20 of 20 joint-capture strata, and the ladder adjusts for alpha_c explicitly and
still finds +8.31 nats on Subclone.  ⇒ division of labour: this figure is the
communication object, the ladder is the controlled measurement.  Say so in the
caption rather than implying the figure does everything.

Usage: 81_modelfree_neighbours.py <arm> [--mincl 100] [--k 20] [--nsplit 5] [--seed 11]
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
MINCL, KNN = arg("--mincl", 100, int), arg("--k", 20, int)
NSPLIT, SEED = arg("--nsplit", 5, int), arg("--seed", 11, int)

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
X = (~Y).astype(np.uint8)                       # 1 = missing
n, nC = X.shape[0], int(g.max()) + 1
o = np.argsort(g, kind="stable"); X, g, codes_all = X[o], g[o], codes_all[o]
off = np.concatenate([[0], np.cumsum(np.bincount(g))]); csz = np.diff(off)
print(f"{arm}: {n:,} cells, {nC} clones (>= {MINCL} cells, sizes {csz.min()}-{csz.max()}), "
      f"k = {KNN}, {NSPLIT} tape splits", flush=True)

# counts[series][clone][j] and hits[series][clone][j], j = 0..k
cnt = {t: np.zeros((nC, KNN + 1), np.int64) for t in ("rel", "rnd")}
hit = {t: np.zeros((nC, KNN + 1), np.int64) for t in ("rel", "rnd")}
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    XB = np.ascontiguousarray(X[:, Bh])          # dropout: half B only
    for b_ in range(nC):
        lo, hi = off[b_], off[b_ + 1]; m = hi - lo
        cb = codes_all[lo:hi]
        agree = np.zeros((m, m), np.int16); denom = np.zeros((m, m), np.int16)
        for a_ in A:                             # relatedness: half A only
            v0 = cb[:, a_, 0] >= 0
            denom += v0[:, None] & v0[None, :]
            for d in range(6):
                cdv = cb[:, a_, d]; ok = cdv >= 0
                agree += (cdv[:, None] == cdv[None, :]) & ok[:, None] & ok[None, :]
        rel = np.where(denom > 0, agree / np.maximum(denom, 1), -1.0)
        np.fill_diagonal(rel, -np.inf)           # never its own neighbour
        keff = min(KNN, m - 1)
        nbb = np.argpartition(-rel, keff - 1, axis=1)[:, :keff]
        if keff < KNN:
            nbb = nbb[:, np.arange(KNN) % keff]
        del rel, agree, denom
        loc = np.arange(m)
        dr = rs.integers(0, m - 1, size=(m, KNN)); dr += (dr >= loc[:, None])
        Xb = XB[lo:hi]
        for tag, IDX in (("rel", nbb), ("rnd", dr)):
            j = Xb[IDX].sum(1).astype(np.int16)   # [cells, |B|], values 0..KNN
            for v in range(KNN + 1):
                mm = j == v
                cnt[tag][b_, v] += int(mm.sum())
                hit[tag][b_, v] += int(Xb[mm].sum())
    print(f"  split {s+1}/{NSPLIT} [{time.time()-t0:.0f}s]", flush=True)

out = dict(arm=arm, mincl=MINCL, k=KNN, nsplit=NSPLIT, n_cells=int(n),
           n_clones=int(nC), clone_sizes=[int(v) for v in csz],
           counts={t: cnt[t].tolist() for t in cnt}, hits={t: hit[t].tolist() for t in hit})
(RES / f"modelfree_{arm}_mincl{MINCL}_k{KNN}.json").write_text(json.dumps(out))
print(f"\n  observed dropout rate by NUMBER of the {KNN} neighbours missing the tape")
print(f"  {'j':>4}" + "".join(f"{lab:>14}" for lab in
      ("relatives %", "n", "arbitrary %", "n")))
for v in range(KNN + 1):
    cr, hr = cnt["rel"][:, v].sum(), hit["rel"][:, v].sum()
    cd, hd = cnt["rnd"][:, v].sum(), hit["rnd"][:, v].sum()
    print(f"  {v:>4}"
          f"{(f'{100*hr/cr:.1f}' if cr else '-'):>14}{cr:>14,}"
          f"{(f'{100*hd/cd:.1f}' if cd else '-'):>14}{cd:>14,}")
print(f"\nwrote results/modelfree_{arm}_mincl{MINCL}_k{KNN}.json")
