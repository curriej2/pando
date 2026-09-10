#!/usr/bin/env python3
r"""
================================================================================
 Per-combo nulls on the SOFT statistic: significance and completeness, separated
================================================================================

Successor to 53_percombo.py.  Same clades, same cells, same p~, same 1,000
permutations (same SEED, so the two runs are PAIRED permutation by permutation).
Three things change:

 1. THE DETECTOR IS THE ONE-SIDED SCORE TEST z = (k - E)/sqrt(V).
    53 scored a POINT hypothesis pinned at pi = 1-eps, which fuses "is this
    non-exchangeable?" with "is the loss total?".  Detecting with it and then
    reporting pi_hat ~ 1.00 is partly circular: the statistic SELECTED for it.
    So detection must be agnostic to completeness.

    ⚠⚠ CORRECTION (2026-09-09).  The first version of this script used
    Lambda_soft as that agnostic detector.  THAT WAS WRONG, and the Mouse3 run
    proved it: of 340 events, 226 were clades of 50-199 cells with k - E = +0.27
    cells and z = 0.08 -- no rate elevation whatsoever.  The reason is that
    H0 and H1_soft are NON-NESTED (H0 is a product of DIFFERENT Bernoullis,
    H1_soft a product of IDENTICAL ones), so they differ in two respects at
    once -- the level AND the homogeneity -- and Lambda_soft decomposes as

        Lambda_soft = m*KL(pi_hat || ebar)  +  [ l(ebar*1) - l(ptilde) ]
                    =      T_elev          +        T_misfit
                       (what we want)         (dispersion of ptilde in the
                                               clade; E_H0[.] = -SUM_c
                                               KL(p_c||ebar) <= 0, grows with
                                               m and with the spread)

    Real clades are MORE homogeneous in capture quality than a random subset of
    their clone (lineage correlates with capture, and clade membership requires
    reaching depth d on the anchor), so T_misfit^obs beats T_misfit^perm and a
    positive margin appears with T_elev = 0.  The permutation was correctly
    calibrated; the ALTERNATIVE was a mixture.  Verified: on the 226 artefacts
    T_elev carried 0.1% of the median margin.

    THE FIX is a NESTED one-parameter family -- shift the log-odds:

        H1(delta): X_c ~ Bern( sigma(eta_c + delta) ),  eta_c = logit(ptilde_c)
                                                              = alpha_c + beta_z
                                                                + gamma_{C,z}

    delta = 0 IS H0 exactly, in the interior; delta -> +inf is total loss; the
    per-cell heterogeneity is kept in BOTH models, so T_misfit cannot arise.
    Score and information at delta = 0:

        S = dl/ddelta|_0 = SUM_c (X_c - ptilde_c) = k - E
        I = -d2l/ddelta2|_0 = SUM_c ptilde_c(1-ptilde_c) = V     (= Var_H0[S])
        z = S/sqrt(I) = (k - E)/sqrt(V)

    i.e. simply (observed - expected)/(SD of the count).  E_H0[S] = 0 and
    Var_H0[S] = V hold EXACTLY, for any m and any heterogeneity, so the
    heterogeneity is absorbed into the SCALE and never into the LOCATION.  And
    since p(1-p) is concave, a homogeneous clade has the LARGEST V hence the
    SMALLEST |z| -- the artefact direction is suppressed, not amplified.

    ⚑ ALL THREE detectors are accumulated in this one scan (z, Lambda_soft,
    Lambda_hard) from the same four bincounts and the SAME permutations, so the
    merge can put the fix and the artefact side by side on identical draws.

 2. pi_hat IS RETAINED per combo (as the raw count k; pi_hat = k/m).

 3. NO PERMUTATIONS ARE HELD OUT.  53 carved HOLD=3 draws out of B in advance to
    calibrate the margin's own null -- which forced the choice of HOLD before
    launching.  Unnecessary: calibration draws can be APPENDED afterwards.  Let
    M = max over permutations 0..B-1.  Compare the real data to M, and later
    compare FRESH permutations B, B+1, ... to the SAME stored M.  Under H0,
    swapping "observed" for "draw B+j" leaves M untouched (it depends only on
    0..B-1), so both are "a value outside the max-set minus the max of that
    set" -- exactly the exchangeability the hold-out design bought, with the
    number of calibration draws now a POST-HOC decision costing one scan each
    (~0.1% of the run).  Run mode --calib does exactly that.

BOTH STATISTICS ARE ACCUMULATED (mx_s and mx_h) so the hard-vs-soft choice is
also deferred to the merge -- three extra bytes per slot.  And by the identity

    Lambda_hard = Lambda_soft - m * KL(pi_hat || 1-eps)          [exact]

the observed Lambda_hard is recoverable for ANY eps from (Lambda_soft, k, m),
so eps is not a commitment either.

--------------------------------------------------------------------------------
 THE STATISTICS, per (clade S, tape z), from FOUR bincounts
--------------------------------------------------------------------------------
    k = SUM_{c in S} X_cz                                (missing count)
    E = SUM_{c in S} p~_cz                               (expected count)
    V = SUM_{c in S} p~_cz (1 - p~_cz)                   (variance of the count)
    U = SUM_{c in S} [ X log p~ + (1-X) log(1-p~) ]      (H0 log-likelihood)

    pi_hat = k/m,  ebar = E/m
    z           = (k - E) / sqrt(V)                       <-- THE DETECTOR
    Lambda_hard = k log(1-eps) + (m-k) log(eps) - U       <-- the ATTRIBUTOR
    Lambda_soft = k log(k/m) + (m-k) log((m-k)/m) - U     <-- kept for comparison

ONE-SIDEDNESS IS AUTOMATIC and needs no truncation.  The margin is
z_obs - max_b z_b, and by exchangeability P(z_obs > max of B draws) = 1/(B+1)
under H0, so only combos whose observed score exceeds EVERY permutation earn a
positive margin -- the upper tail, by construction.  (The previous version had
to map the wrong side to -inf, which poisoned the mean/sd accumulators; that
whole complication is gone.)  A called event is additionally required to have
z_obs > 0, since negative z is a DEFICIT of missingness and nothing predicts
heritable excess recovery.

⚠ V > 0 always, since p~ is clipped to [1e-6, 1-1e-6]; guarded anyway.

Usage:
  accumulate : 62_percombo_soft.py <arm> [--nperm 1000] [--permpart i/N]
                                   [--maxd 6] [--seed S] [--eps 0.01]
  calibrate  : 62_percombo_soft.py <arm> --calib 12 [--maxd 6] [--seed S]
  observed   : 62_percombo_soft.py <arm> --nperm 0   (obs only; for validation)

Output: results/percombo_soft/percombo_soft_{arm}_d{D}_p{i}.npz   (accumulators)
        results/percombo_soft/percombo_soft_{arm}_d{D}_cal.npz    (calib draws)
================================================================================
"""
import sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4


def arg(flag, default, cast=str):
    return cast(sys.argv[sys.argv.index(flag) + 1]) if flag in sys.argv else default


arm = sys.argv[1]
EPS = arg("--eps", 0.01, float)
NPERM = arg("--nperm", 1000, int)
SEED = arg("--seed", 20260903, int)
PART = arg("--permpart", "1/1", str)
MAXD = arg("--maxd", 6, int)
NCAL = arg("--calib", 0, int)
DEPTHS = list(range(1, MAXD + 1))
LOG1ME, LOGE = np.log(1 - EPS), np.log(EPS)

# ---------------------------------------------------------------- cells, p~
z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
cache = RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4 else f"prefix_codes_{arm}.npz")
codes = np.load(cache, allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0], (codes.shape, Y.shape)
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
miss = ~Y
print(f"{arm}: {n:,} cells, {K} tapes, {g_clone.max()+1:,} clones, "
      f"depths 1-{MAXD}, eps={EPS}", flush=True)


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
eta = np.log(P / (1 - P))                       # third margin: per (clone, tape)
order0 = np.argsort(g_clone, kind="stable")
bnd0 = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])
for a_, b_ in zip(bnd0[:-1], bnd0[1:]):
    ix = order0[a_:b_]
    e = eta[ix]; tgt = miss[ix].sum(0).astype(float); gam = np.zeros(K)
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)

MISS = miss.astype(np.float64)                              # k
UU = np.where(miss, np.log(P), np.log1p(-P))                # U
PP = P                                                      # E
PQ = P * (1.0 - P)                                          # V

# ------------------------------------------------- block layout, built ONCE
# `codes` and `g_clone` are never permuted, so every block's clade partition --
# and therefore each combo's global slot index -- is identical in the observed
# data and in every permutation.  Byte-identical to 53_percombo.py.
blocks, off = [], 0
for d in DEPTHS:
    for a_ in range(K):
        cd = codes[:, a_, d - 1]
        ok = cd >= 0
        if ok.sum() < MIN_CLADE:
            continue
        idx = np.flatnonzero(ok)
        key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
        _, sub = np.unique(key, return_inverse=True)
        G = int(sub.max() + 1)
        ns = np.bincount(sub, minlength=G)
        keep = np.flatnonzero(ns >= MIN_CLADE).astype(np.int32)
        if keep.size == 0:
            continue
        valid = np.ones((keep.size, K), bool)
        valid[:, a_] = False                     # never score the anchor itself
        blocks.append((a_, idx.astype(np.int32), sub.astype(np.int32), G, off,
                       valid.ravel(), keep, ns[keep].astype(np.int32)))
        off += keep.size * K
NC = off
SIZE = np.concatenate([np.repeat(b[7], K) for b in blocks])
assert SIZE.max() < 65536, f"clade of {SIZE.max()} cells overflows uint16"
VALID = np.concatenate([b[5] for b in blocks])
print(f"  {len(blocks):,} blocks, {NC:,} combo slots "
      f"({100*VALID.sum()/max(NC,1):.1f}% scorable)", flush=True)


def lam_all(mm, uu, pp, pq, out_z, out_s, out_h, out_k, out_e=None, out_v=None):
    """Four bincounts per (block, tape) -> the score z, Lambda_soft,
    Lambda_hard and the missing count k, for every kept combo slot."""
    for a_, idx, sub, G, o, _, keep, ns in blocks:
        nk = keep.size
        k = np.empty((G, K)); U = np.empty((G, K))
        E = np.empty((G, K)); Vv = np.empty((G, K))
        mi, ui, pi_, qi_ = mm[idx], uu[idx], pp[idx], pq[idx]
        for t in range(K):
            k[:, t] = np.bincount(sub, weights=mi[:, t], minlength=G)
            U[:, t] = np.bincount(sub, weights=ui[:, t], minlength=G)
            E[:, t] = np.bincount(sub, weights=pi_[:, t], minlength=G)
            Vv[:, t] = np.bincount(sub, weights=qi_[:, t], minlength=G)
        k = k[keep]; U = U[keep]; E = E[keep]; Vv = Vv[keep]
        m = ns[:, None].astype(float)
        # the detector: (observed - expected) / SD of the count
        zz = (k - E) / np.sqrt(np.maximum(Vv, 1e-12))
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = np.where(k > 0, k * np.log(np.maximum(k, 1e-300) / m), 0.0)
            t2 = np.where(m - k > 0, (m - k) * np.log(np.maximum(m - k, 1e-300) / m), 0.0)
        ls = t1 + t2 - U
        lh = k * LOG1ME + (m - k) * LOGE - U
        sl = slice(o, o + nk * K)
        out_z[sl] = zz.ravel()
        out_s[sl] = ls.ravel()
        out_h[sl] = lh.ravel()
        out_k[sl] = k.ravel()
        if out_e is not None:
            out_e[sl] = E.ravel(); out_v[sl] = Vv.ravel()


# ⚠ the observed pass is not needed in calib mode (28 B/slot; 3.0 GB on Pre-TX)
if NCAL == 0:
    obs_z = np.empty(NC, np.float32); obs_s = np.empty(NC, np.float32)
    obs_h = np.empty(NC, np.float32); obs_k = np.empty(NC, np.float64)
    obs_e = np.empty(NC, np.float32); obs_v = np.empty(NC, np.float32)
    lam_all(MISS, UU, PP, PQ, obs_z, obs_s, obs_h, obs_k, obs_e, obs_v)

order = np.argsort(g_clone, kind="stable")
bnds = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])


def perm_inv(b):
    """Within-clone permutation of cell ROWS.  Whole rows travel, so a cell's
    missingness and its p~ move together and all three margins are preserved."""
    rb = np.random.default_rng([SEED, b])
    pm = order.copy()
    for a_, b_ in zip(bnds[:-1], bnds[1:]):
        pm[a_:b_] = rb.permutation(pm[a_:b_])
    inv = np.empty_like(pm); inv[order] = pm
    return inv


outd = RES / "percombo_soft"; outd.mkdir(exist_ok=True)
tag = f"{arm}_d{MAXD}"
cur_z = np.empty(NC, np.float32); cur_s = np.empty(NC, np.float32)
cur_h = np.empty(NC, np.float32); cur_k = np.empty(NC, np.float64)

# ---------------------------------------------------------------- CALIB MODE
if NCAL > 0:
    # Fresh permutations B .. B+NCAL-1 -- disjoint from the accumulated set, so
    # they are exchangeable with the observed data given max over 0..B-1.
    # ⚠ MEMORY.  Pre-TX has 105.9M slots, so one (NCAL x NC) float32 is 5.1 GB.
    # z is the detector and needs every draw; Lambda_soft/Lambda_hard are only ever
    # used for the side-by-side comparison table, so NCMP draws of those are plenty.
    NCMP = min(NCAL, 3)
    ZZ = np.empty((NCAL, NC), np.float32)
    LS = np.empty((NCMP, NC), np.float32)
    LH = np.empty((NCMP, NC), np.float32)
    KK = np.empty((NCMP, NC), np.uint16)              # k only for the plane's null cloud
    t0 = time.time()
    for j in range(NCAL):
        b = NPERM + j
        iv = perm_inv(b)
        lam_all(MISS[iv], UU[iv], PP[iv], PQ[iv], cur_z, cur_s, cur_h, cur_k)
        ZZ[j] = cur_z
        if j < NCMP:
            LS[j] = cur_s; LH[j] = cur_h; KK[j] = cur_k.astype(np.uint16)
        print(f"    calib draw {j} (perm {b}) [{time.time()-t0:.0f}s]", flush=True)
    fn = outd / f"percombo_soft_{tag}_cal.npz"
    np.savez_compressed(fn, z=ZZ, lam_s=LS, lam_h=LH, k=KK,
                        meta=np.array([NPERM, NCAL, SEED, MAXD, NC, EPS, NCMP], np.float64))
    print(f"wrote {fn.relative_to(ROOT)}  ({NCAL} calibration draws, {time.time()-t0:.0f}s)")
    sys.exit(0)

# ------------------------------------------------------------ ACCUMULATE MODE
ip, npart = (int(x) for x in PART.split("/"))
assert 1 <= ip <= npart, PART
lo, hi = ((ip - 1) * NPERM) // npart, (ip * NPERM) // npart
mx_z = np.full(NC, -np.inf, np.float32)
mx_s = np.full(NC, -np.inf, np.float32)
mx_h = np.full(NC, -np.inf, np.float32)
n_ge_z = np.zeros(NC, np.int16); n_ge_h = np.zeros(NC, np.int16)
s1 = np.zeros(NC, np.float32); s2 = np.zeros(NC, np.float32)
print(f"  part {ip}/{npart}: permutations {lo}..{hi-1} of B={NPERM}", flush=True)
t0 = time.time()
for j, b in enumerate(range(lo, hi)):
    iv = perm_inv(b)
    lam_all(MISS[iv], UU[iv], PP[iv], PQ[iv], cur_z, cur_s, cur_h, cur_k)
    np.maximum(mx_z, cur_z, out=mx_z)
    np.maximum(mx_s, cur_s, out=mx_s)
    np.maximum(mx_h, cur_h, out=mx_h)
    n_ge_z += (cur_z >= obs_z)
    n_ge_h += (cur_h >= obs_h)
    s1 += cur_z; s2 += cur_z * cur_z          # z is always finite: clean sums
    if (j + 1) % 25 == 0 or j == 0:
        print(f"    perm {b} [{j+1} done, {time.time()-t0:.0f}s]", flush=True)

fn = outd / f"percombo_soft_{tag}_p{ip:03d}.npz"
np.savez_compressed(
    fn,
    obs_z=obs_z, obs_s=obs_s, obs_h=obs_h,
    obs_k=obs_k.astype(np.uint16), obs_e=obs_e, obs_v=obs_v,
    size=SIZE.astype(np.uint16), valid=VALID,
    mx_z=np.where(np.isfinite(mx_z), mx_z, -1e30).astype(np.float32),
    mx_s=np.where(np.isfinite(mx_s), mx_s, -1e30).astype(np.float32),
    mx_h=np.where(np.isfinite(mx_h), mx_h, -1e30).astype(np.float32),
    n_ge_z=n_ge_z, n_ge_h=n_ge_h, s1=s1, s2=s2,
    meta=np.array([NPERM, npart, ip, SEED, hi - lo, MAXD, NC, EPS], np.float64))
print(f"wrote {fn.relative_to(ROOT)}  ({hi-lo} permutations, {time.time()-t0:.0f}s)")
