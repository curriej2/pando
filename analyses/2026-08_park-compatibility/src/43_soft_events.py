#!/usr/bin/env python3
r"""
================================================================================
 Hard vs soft: complete losses, or graded rate shifts too?
================================================================================

40_event_catalogue.py asks a POINT question -- "was this tape lost entirely on
this clade's stem?" -- by pinning H1 at P(missing) = 1 - eps with NO free
parameter.  That is deliberately strict, and it makes real partial silencing
invisible.  Worked, for a 20-cell clade the model expects 30% missing:

    16/20 missing (a 0.30 -> 0.80 shift)  Lambda_hard =  +2.11 nats  -> rejected
    20/20 missing (a 0.30 -> 1.00 shift)  Lambda_hard = +23.88 nats  -> kept

An 11x difference in evidence for a modest difference in severity.

THE SOFT VARIANT gives the clade its own rate pi, fitted:

    Lambda_soft = [k log(k/m) + (m-k) log((m-k)/m)]  -  [SUM_miss log p~ + SUM_pres log(1-p~)]

Same 20-cell examples: 16/20 -> +10.68 nats (now detectable), 20/20 -> +24.08.

WHY pi_hat = k/m.  Under H1_soft every cell of the clade shares one probability pi.
Maximising k log pi + (m-k) log(1-pi) gives k/pi = (m-k)/(1-pi) => pi_hat = k/m --
the observed fraction missing.  A likelihood-ratio test gives the alternative its
BEST shot, so the fitted value is what goes in.  The hard test, by contrast, is a
POINT hypothesis (one specific claim: total loss), which is why it has no fitted
parameter and why Lambda_soft >= Lambda_hard always: H1_hard is the special case
pi = 1-eps of H1_soft, and pi_hat maximises over all pi.

⚠ CORRECTION to what I said last session: Wilks does NOT apply here.  H0 (per-cell
p~, no free parameters) is nested inside H1_soft only when every p~ in the clade is
equal, which it is not.  The models are non-nested, Lambda_soft can even be
NEGATIVE (if the varying-p~ model tracks the pattern better than any constant pi),
and there is no chi^2_1 reference.  The permutation null is not merely prudent, it
is the only calibration available.

ONE-SIDED.  pi_hat can fall BELOW the expected rate -- a clade with fewer missing
than predicted.  Nothing predicts heritable excess RECOVERY, so the catalogue keeps
only pi_hat > expected, and the count in the opposite direction is a null
calibration check (it should be ~0), exactly as the significantly-negative tapes
served in script 37.

EPSILON.  ⚠ Second correction: "fitting eps" was the wrong word.  Any estimate from
called events is selected on having few present cells, so it is biased downward.
What is honest is (a) a SENSITIVITY SWEEP over eps, reported, and (b) the pi_hat
distribution from the soft run read descriptively, with the selection bias stated.
Both are produced here.

BOTH STATISTICS COME FROM ONE SCAN.  Per (clade, tape) we need four group sums:
    k = SUM miss,  A = SUM miss*log(p~),  B = SUM (1-miss)*log(1-p~),  E = SUM p~
    Lambda_hard = k log(1-eps) + (m-k) log eps - (A+B)
    Lambda_soft = k log(k/m)  + (m-k) log((m-k)/m) - (A+B)

Usage: 43_soft_events.py <arm> [--nperm 5] [--lam 10]
Output: results/soft_events_{arm}.json
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
DEPTHS, MIN_CLADE = [1, 2, 3, 4], 4
EPS_SWEEP = [0.005, 0.01, 0.02, 0.05]

arm = sys.argv[1]
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 5
LAM0 = float(sys.argv[sys.argv.index("--lam") + 1]) if "--lam" in sys.argv else 10.0
rng = np.random.default_rng(20260903)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s], clone[s]
codes = np.load(RES / f"prefix_codes_{arm}.npz", allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
miss = ~Y
print(f"{arm}: {n:,} cells, {K} tapes, {Gc:,} clones", flush=True)


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
# the third margin -- per (clone, tape), so Lambda measures WITHIN-clone structure
eta = np.log(P / (1 - P))
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

MISS = miss.astype(float)
ARR1 = MISS * np.log(P)                      # -> A
ARR2 = (1 - MISS) * np.log1p(-P)             # -> B
anchors = np.arange(K)


def scan(mi, a1, a2, pp):
    """Return (lam_soft, lam_hard[eps], pi_hat, exp_rate, size) over all candidates."""
    LS, LH, PI, ER, SZ = [], {e: [] for e in EPS_SWEEP}, [], [], []
    for d in DEPTHS:
        for a_ in anchors:
            cd = codes[:, a_, d - 1]
            ok = cd >= 0
            if ok.sum() < MIN_CLADE:
                continue
            idx = np.flatnonzero(ok)
            key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
            _, sub = np.unique(key, return_inverse=True)
            G = sub.max() + 1
            m = np.bincount(sub, minlength=G).astype(float)
            k = np.empty((G, K)); A = np.empty((G, K)); B = np.empty((G, K)); E = np.empty((G, K))
            for t in range(K):
                k[:, t] = np.bincount(sub, weights=mi[idx, t], minlength=G)
                A[:, t] = np.bincount(sub, weights=a1[idx, t], minlength=G)
                B[:, t] = np.bincount(sub, weights=a2[idx, t], minlength=G)
                E[:, t] = np.bincount(sub, weights=pp[idx, t], minlength=G)
            mm = m[:, None]
            with np.errstate(divide="ignore", invalid="ignore"):
                t1 = np.where(k > 0, k * np.log(np.maximum(k, 1e-300) / mm), 0.0)
                t2 = np.where(mm - k > 0, (mm - k) * np.log(np.maximum(mm - k, 1e-300) / mm), 0.0)
            ls = t1 + t2 - (A + B)
            pi = k / mm; er = E / mm
            keep = (m >= MIN_CLADE)[:, None] & np.ones((1, K), bool)
            keep[:, a_] = False
            keep &= pi > er                              # one-sided: excess only
            gg, zz = np.nonzero(keep & (ls >= LAM0))
            if gg.size:
                LS.append(ls[gg, zz]); PI.append(pi[gg, zz])
                ER.append(er[gg, zz]); SZ.append(mm[gg, 0])
                for e_ in EPS_SWEEP:
                    LH[e_].append(k[gg, zz] * np.log(1 - e_) +
                                  (mm[gg, 0] - k[gg, zz]) * np.log(e_) - (A + B)[gg, zz])
    cat = lambda L: np.concatenate(L) if L else np.zeros(0)
    return cat(LS), {e: cat(v) for e, v in LH.items()}, cat(PI), cat(ER), cat(SZ)


ls, lh, pi, er, sz = scan(MISS, ARR1, ARR2, P)
print(f"  soft candidates at Lambda >= {LAM0}: {ls.size:,}", flush=True)

order = np.argsort(g_clone, kind="stable")
bnds = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])
null_ls = []
for b in range(NPERM):
    pm = order.copy()
    for a_, b_ in zip(bnds[:-1], bnds[1:]):
        pm[a_:b_] = rng.permutation(pm[a_:b_])
    inv = np.empty_like(pm); inv[order] = pm
    nls, _, _, _, _ = scan(MISS[inv], ARR1[inv], ARR2[inv], P[inv])
    null_ls.append(nls)
    print(f"    perm {b+1}/{NPERM}: {nls.size:,} soft candidates", flush=True)

grid = np.unique(np.round(np.geomspace(LAM0, max(ls.max(), LAM0 * 2), 60), 2))
on = np.array([(ls >= t).sum() for t in grid], float)
nn = np.array([np.mean([(v >= t).sum() for v in null_ls]) for t in grid])
fdrc = np.where(on > 0, nn / np.maximum(on, 1), 1.0)
ok = np.flatnonzero(fdrc <= 0.05)
LAM = float(grid[ok[0]]) if ok.size else float(grid[-1])
FDR = float(fdrc[ok[0]]) if ok.size else float(fdrc[-1])
print(f"  FDR curve: " + "  ".join(f"{t:.0f}n:{100*f:.0f}%" for t, f in zip(grid[::6], fdrc[::6])))
print(f"  => soft threshold {LAM:.1f} nats at FDR {100*FDR:.1f}%", flush=True)

sel = ls >= LAM
out = {"arm": arm, "lambda_soft_threshold": LAM, "fdr": FDR, "n_perm": NPERM,
       "n_soft_candidates": int(sel.sum()), "eps_sweep": {}}
for e_ in EPS_SWEEP:
    hard_pass = (lh[e_] >= LAM0) & sel
    out["eps_sweep"][str(e_)] = {
        "n_hard_within_soft": int(hard_pass.sum()),
        "frac_of_soft": float(hard_pass.mean()) if sel.sum() else 0.0,
        "median_pi_hard": float(np.median(pi[hard_pass])) if hard_pass.sum() else None}
if sel.sum():
    p_ = pi[sel]
    out["pi_quantiles"] = {str(q): float(np.quantile(p_, q))
                           for q in (0.1, 0.25, 0.5, 0.75, 0.9)}
    out["pi_hist"] = np.histogram(p_, bins=np.linspace(0, 1, 21))[0].tolist()
    out["frac_pi_ge_099"] = float((p_ >= 0.99).mean())
    out["frac_pi_lt_090"] = float((p_ < 0.90).mean())
    out["median_expected_rate"] = float(np.median(er[sel]))
    out["median_clade_size"] = float(np.median(sz[sel]))
    print(f"\n  soft events {int(sel.sum()):,}: pi_hat median {np.median(p_):.3f} "
          f"(expected {np.median(er[sel]):.3f}), clade size median {np.median(sz[sel]):.0f}")
    print(f"  pi_hat >= 0.99: {100*out['frac_pi_ge_099']:.1f}%  |  "
          f"pi_hat < 0.90 (PARTIAL): {100*out['frac_pi_lt_090']:.1f}%")
    print("  eps sensitivity -- share of soft events the HARD test also keeps:")
    for e_ in EPS_SWEEP:
        d = out["eps_sweep"][str(e_)]
        print(f"    eps={e_:<6} {d['n_hard_within_soft']:>8,}  "
              f"({100*d['frac_of_soft']:.1f}% of soft)", flush=True)
(RES / f"soft_events_{arm}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/soft_events_{arm}.json")
