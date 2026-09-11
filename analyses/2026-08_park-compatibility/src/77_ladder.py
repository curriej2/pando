#!/usr/bin/env python3
r"""
THE NESTED LADDER -- how much of dropout does each explanation buy, out of sample?

Spec agreed with Justin 2026-09-11.  Read README "The nested ladder" first: two
rounds of correction are recorded there and both fixed something wrong.

THE ARGUMENT.  Every rung is the next objection an audience raises, and the last
one is the punchline.  Everything is fitted INSIDE A CLONE, so clone, batch,
harvest and mouse are held fixed by construction rather than by trusting a
parameter to have absorbed them.

  M0  logit p = mu                  dropout happens at some rate
  M1  logit p = mu_C                some clones are worse than others
  M2  logit p = mu_C + alpha_c      inside a clone, some cells are badly captured
  M3  logit p = alpha_c + beta_z    inside a clone, some tapes are badly recovered
                                    (clone-wide silencing lives HERE)
  M4  logit p = alpha_c + beta_z + w*u_cz    ...and your relatives still predict you

M3 properly contains M1 and M2 because beta_z absorbs mu_C.

⚠⚠ IDENTIFIABILITY: one constraint PER CLONE, sum_{c in C} alpha_c = 0.  Adding d
to every alpha and subtracting it from every beta changes no probability, and d is
free separately in each clone, so a single global constraint leaves n_clones - 1
ridges unpinned.

⚠⚠ TWO-STAGE ON PURPOSE.  M0-M3 are fresh MLE fits.  Then M3 is FROZEN, u is built
from its residuals, and M4 fits ONLY w on a fixed offset.  (i) removes a
circularity -- u is built from M3's residuals, so M3 must not then move; (ii) it is
CONSERVATIVE, since a constrained optimum cannot exceed an unconstrained one, so
if alpha/beta already absorbed lineage signal then M3 keeps the credit and the
lineage rung understates itself; (iii) "one parameter" becomes literal.
⇒ Delta_4 is a LOWER BOUND on the likelihood-ratio gain.  Report it as one.

⚠⚠ THE NULL IS "k RANDOM CLONE-MATES, SELF EXCLUDED", not a residual permutation.
The M3 score equation forces sum_{c in C} r_cz = 0, so any self-excluding
neighbour set carries -1/(n_C - 1) of the cell's own residual.  74's null permutes
residual rows across fixed neighbour slots, which produces a uniform subset of ALL
n_C cells INCLUDING c, expectation 0, so it does not carry that term and obs-null
fails to cancel it.  Harmless at n_C = 3387 (-0.0003) and fatal at n_C = 31
(-0.0333): the first run of this ladder returned a NEGATIVE lineage rung on Pre-TX
for exactly this reason.  The self-excluding null also happens to be the question
we mean: are the NEAREST relatives more informative than ARBITRARY clone-mates?

WHY HELD-OUT ENTRIES AND NOT CELLS.  M3 carries thousands of parameters and
in-sample likelihood always rises when parameters are added, so an in-sample
ladder measures flexibility, not truth.  Only a real held-out set lets a rung
LOSE, and that possibility is what makes the comparison mean anything.

THE SCORE EQUATIONS ARE THE ENGINE.  At the M3 optimum both margins are matched
EXACTLY on training entries, so the model cannot be surprised by how many tapes a
cell is missing nor by how many cells are missing a tape.  The only thing left is
WHICH cells are missing WHICH tapes -- the pattern, not the abundance.

Usage: 77_ladder.py <arm> [--mincl 100] [--nsplit 5] [--nmask 2] [--k 20]
       [--frac 0.5] [--sdfloor 0.01] [--seed 11] [--tag ""]
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
LO, HI = 1e-12, 1 - 1e-12


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
MINCL = arg("--mincl", 100, int)
NSPLIT, NMASK = arg("--nsplit", 5, int), arg("--nmask", 2, int)
KNN = arg("--k", 20, int)
FRAC = arg("--frac", 0.5, float)
PFLOOR = arg("--sdfloor", 0.01, float)
SEED = arg("--seed", 11, int)
TAG = arg("--tag", "", str)

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
# ⚠ the prefix cache is built on the ALREADY tape-filtered rows (44_prefix_codes6.py
# applies THR itself), so `sel` must NOT be applied to it -- only the clone filter.
codes_all = np.load(RES / f"prefix_codes6_{arm}.npz", allow_pickle=False)["codes"]
assert codes_all.shape[0] == keep.size, (codes_all.shape[0], keep.size)
codes_all = codes_all[keep]
X = (~Y).astype(np.float64)                      # 1 = missing
n, nC = X.shape[0], int(g.max()) + 1
order = np.argsort(g, kind="stable")
X, g, codes_all = X[order], g[order], codes_all[order]
off = np.concatenate([[0], np.cumsum(np.bincount(g))])
csz = np.diff(off)
print(f"{arm}: {n:,} cells in {nC} clones with >= {MINCL} cells "
      f"(sizes {csz.min()}-{csz.max()}, median {int(np.median(csz))})", flush=True)
print(f"  k = {KNN}, {NSPLIT} tape splits x {NMASK} entry masks, train fraction {FRAC}", flush=True)

SDF = np.sqrt(PFLOOR * (1 - PFLOOR))
sig = lambda e: 1.0 / (1.0 + np.exp(-e))


def ll(Xv, pv):
    p = np.clip(pv, LO, HI)
    return float(np.sum(Xv * np.log(p) + (1 - Xv) * np.log1p(-p)))


def fit_m3(Xb, Mb, iters=80):
    """alpha_c + beta_z inside one clone, on TRAINING entries only.
    Newton steps alternating the two blocks.  alpha centred each sweep -- that is
    the per-clone identifiability constraint, and it changes no probability."""
    m, q = Xb.shape
    cw = Mb.sum(1); tw = Mb.sum(0)
    r0 = np.clip((Mb * Xb).sum(0) / np.maximum(tw, 1), 1e-3, 1 - 1e-3)
    be = np.log(r0 / (1 - r0)); al = np.zeros(m)
    for _ in range(iters):
        p = sig(al[:, None] + be[None, :])
        gr = (Mb * (Xb - p)).sum(1); he = (Mb * p * (1 - p)).sum(1)
        al = np.clip(al + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -12, 12)
        al -= al.mean()                           # sum_{c in C} alpha_c = 0
        p = sig(al[:, None] + be[None, :])
        gr = (Mb * (Xb - p)).sum(0); he = (Mb * p * (1 - p)).sum(0)
        be = np.clip(be + np.clip(gr / np.maximum(he, 1e-9), -2, 2), -12, 12)
    return al, be, cw, tw


def fit_w(Xv, ev, uv, iters=80):
    """One scalar on a FIXED offset.  w = 0 recovers M3 exactly."""
    w = 0.0
    for _ in range(iters):
        q = sig(ev + w * uv)
        gr = float(np.dot(uv, Xv - q)); he = float(np.dot(uv * uv, q * (1 - q)))
        if he < 1e-12:
            break
        st = gr / he; w += float(np.clip(st, -2, 2))
        if abs(st) < 1e-11:
            break
    return w


runs = []
t0 = time.time()
for s in range(NSPLIT):
    rs = np.random.default_rng([SEED, s])
    perm = rs.permutation(K)
    A, Bh = np.sort(perm[:K // 2]), np.sort(perm[K // 2:])
    XB = X[:, Bh]; nb = Bh.size
    # --- neighbours depend on A only, so build once per tape split ------------
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
        np.fill_diagonal(rel, -np.inf)            # never its own neighbour
        keff = min(KNN, m - 1)
        nbb = np.argpartition(-rel, keff - 1, axis=1)[:, :keff]
        if keff < KNN:                            # pad by repetition, mean unchanged
            nbb = nbb[:, np.arange(KNN) % keff]
        NBR[lo:hi] = nbb + lo
        del rel, agree, denom
    for r in range(NMASK):
        rm = np.random.default_rng([SEED, s, r, 99])
        M = (rm.random((n, nb)) < FRAC)
        Tr, Va = M, ~M
        if Va.sum() == 0 or Tr.sum() == 0:
            continue
        P = {j: np.zeros((n, nb)) for j in range(5)}
        # M0 -- one global rate
        P[0][:] = np.clip((Tr * XB).sum() / Tr.sum(), 1e-6, 1 - 1e-6)
        for b_ in range(nC):
            lo, hi = off[b_], off[b_ + 1]
            Xb, Mb = XB[lo:hi], Tr[lo:hi].astype(np.float64)
            # M1 -- one rate per clone
            P[1][lo:hi] = np.clip((Mb * Xb).sum() / max(Mb.sum(), 1), 1e-6, 1 - 1e-6)
            # M2 -- one rate per cell (mu_C + alpha_c reparameterised)
            cr = (Mb * Xb).sum(1) / np.maximum(Mb.sum(1), 1)
            P[2][lo:hi] = np.clip(cr, 1e-6, 1 - 1e-6)[:, None]
            # M3 -- alpha_c + beta_z
            al, be, cw, tw = fit_m3(Xb, Mb)
            P[3][lo:hi] = np.clip(sig(al[:, None] + be[None, :]), 1e-6, 1 - 1e-6)
        # --- M4: FREEZE M3, build u from its residuals, fit w alone -----------
        Pt = P[3]
        ETA = np.log(Pt / (1 - Pt))
        RP = (XB - Pt) / np.maximum(np.sqrt(Pt * (1 - Pt)), SDF)
        Mf = Tr.astype(np.float64)
        res = {}
        # ⚠⚠ THE NULL MUST EXCLUDE SELF.  The M3 score equation forces
        # sum_{c in C} r_cz = 0 over training entries, so a neighbour set that
        # EXCLUDES c carries -1/(n_C - 1) of c's own residual -- a leak of the
        # entry into its own predictor.  Permuting residual ROWS across fixed
        # neighbour SLOTS (what 74 does) yields a uniform k-subset of ALL n_C
        # cells, c included, whose expectation is the clone mean = 0.  That null
        # carries NO leak, so obs - null does not cancel it.  Verified by
        # simulation: observed slope -0.0500 at n_C = 21 against theory
        # -1/(n_C-1) = -0.0500, permutation null +0.0005, self-excluding null
        # -0.0500.  Negligible at n_C = 3387 (-0.0003), which is why large clones
        # were unaffected and Pre-TX (median clone 31) went NEGATIVE.
        # ⇒ the null is "k RANDOM clone-mates, self excluded" -- which is also the
        #   question we actually mean: are the NEAREST relatives more informative
        #   than ARBITRARY ones?
        rnull = np.random.default_rng([SEED, s, r, 1234])
        RND = np.zeros_like(NBR)
        for b_ in range(nC):
            lo, hi = off[b_], off[b_ + 1]; m = hi - lo
            loc = np.arange(m)
            draw = rnull.integers(0, m - 1, size=(m, KNN))
            draw += (draw >= loc[:, None])          # skip self, keeping it uniform
            RND[lo:hi] = draw + lo
        for kind, IDX in (("obs", NBR), ("null", RND)):
            mz = Mf[IDX].sum(1)                    # neighbours with a TRAIN entry
            u = (Mf[IDX] * RP[IDX]).sum(1) / np.maximum(mz, 1.0)
            u[mz == 0] = 0.0
            w = fit_w(XB[Tr], ETA[Tr], u[Tr])
            p4 = np.clip(sig(ETA + w * u), LO, HI)
            res[kind] = (w, ll(XB[Va], p4[Va]) / Va.sum())
        L = [ll(XB[Va], P[j][Va]) / Va.sum() for j in range(4)]
        L.append(res["obs"][1])
        runs.append(dict(split=s, mask=r, ll=L, ll_null4=res["null"][1],
                         w=res["obs"][0], w_null=res["null"][0],
                         n_val=int(Va.sum()), n_train=int(Tr.sum())))
        print(f"  split {s} mask {r}: l0 {L[0]:.4f} l1 {L[1]:.4f} l2 {L[2]:.4f} "
              f"l3 {L[3]:.4f} l4 {L[4]:.4f} | l4_null {res['null'][1]:.4f} "
              f"w {res['obs'][0]:+.3f} [{time.time()-t0:.0f}s]", flush=True)

Lm = np.array([r["ll"] for r in runs])
L4n = np.array([r["ll_null4"] for r in runs])
NAME = ["M0 global rate", "M1 + clone", "M2 + cell", "M3 + tape (in clone)", "M4 + lineage"]
print(f"\n  {len(runs)} replicates ({NSPLIT} splits x {NMASK} masks)")
print(f"\n  {'rung':>24}{'mean loglik':>13}{'delta':>10}{'sd':>8}"
      f"{'nats/cell':>11}{'R2':>9}")
for j in range(5):
    d = Lm[:, j] - (Lm[:, j - 1] if j else 0)
    dd = Lm[:, j] - Lm[:, j - 1] if j else np.zeros(len(runs))
    r2 = 1 - Lm[:, j] / Lm[:, 0]
    print(f"  {NAME[j]:>24}{Lm[:,j].mean():>13.5f}"
          f"{(dd.mean() if j else 0):>+10.5f}{(dd.std(ddof=1) if j and len(runs)>1 else 0):>8.5f}"
          f"{(K*dd.mean() if j else 0):>+11.3f}{r2.mean():>9.4f}")
d4 = Lm[:, 4] - Lm[:, 3]; d4n = L4n - Lm[:, 3]
diff = d4 - d4n
Q = 1 - Lm[:, 4] / Lm[:, 3]
Qn = 1 - L4n / Lm[:, 3]
print(f"\n  ⭐ the lineage rung, observed MINUS null (the only quotable form):")
print(f"     delta_4 obs  {d4.mean():+.5f}   null {d4n.mean():+.5f}   "
      f"obs-null {diff.mean():+.5f} +- {diff.std(ddof=1)/np.sqrt(len(runs)) if len(runs)>1 else 0:.5f}")
print(f"     nats per cell (x{K})          obs-null {K*diff.mean():+.3f}")
print(f"     Q = share of REMAINING deviance removed: obs {100*Q.mean():.2f}%  "
      f"null {100*Qn.mean():.2f}%  obs-null {100*(Q-Qn).mean():.2f}%")
print(f"     w  obs {np.mean([r['w'] for r in runs]):+.3f}   "
      f"null {np.mean([r['w_null'] for r in runs]):+.3f}")
out = dict(arm=arm, mincl=MINCL, k=KNN, nsplit=NSPLIT, nmask=NMASK, frac=FRAC,
           n_cells=int(n), n_clones=int(nC), clone_sizes=[int(v) for v in csz],
           K=int(K), runs=runs)
sfx = f"_mincl{MINCL}_k{KNN}" + (f"_{TAG}" if TAG else "")
(RES / f"ladder_{arm}{sfx}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/ladder_{arm}{sfx}.json")
