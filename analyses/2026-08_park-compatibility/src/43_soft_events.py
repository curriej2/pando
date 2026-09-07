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

--------------------------------------------------------------------------------
 B = 1,000 PERMUTATIONS  (added 2026-09-03)
--------------------------------------------------------------------------------
Same scheme as 40_event_catalogue.py: --permpart i/N runs a slice of the B
permutations and writes only the count vector on a FIXED threshold grid to
results/permparts/permcounts_43_{arm}{_d6}_p{i}.npz, then exits.  Permutation b
is seeded from (SEED, b) so the slices are addable in any order.
46_perm_merge.py pools them; the observed run is redone once with --nullfile.

Usage: 43_soft_events.py <arm> [--nperm 5] [--lam 10] [--maxd 4]
       [--seed S] [--permpart i/N] [--nullfile F]
Output: results/soft_events_{arm}.json
================================================================================
"""
import json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4
EPS_SWEEP = [0.005, 0.01, 0.02, 0.05]

arm = sys.argv[1]
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 5
LAM0 = float(sys.argv[sys.argv.index("--lam") + 1]) if "--lam" in sys.argv else 10.0
# --maxd 6 uses the full-depth prefix cache, so a loss can be localised to the
# finest clade the recorder resolves.  That is the test for whether "partial"
# events are graded silencing or just a complete loss on a clade we scored too
# coarsely: at maximum resolution the second kind must resolve into complete ones.
MAXD = int(sys.argv[sys.argv.index("--maxd") + 1]) if "--maxd" in sys.argv else 4
DEPTHS = list(range(1, MAXD + 1))
SEED = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260903
PART = sys.argv[sys.argv.index("--permpart") + 1] if "--permpart" in sys.argv else None
NULLF = sys.argv[sys.argv.index("--nullfile") + 1] if "--nullfile" in sys.argv else None
TARGET = float(sys.argv[sys.argv.index("--target") + 1]) if "--target" in sys.argv else 0.05
# ⚠⚠ --budget X selects each stratum's threshold by an ABSOLUTE expected-false
# count (null mean <= X per stratum) instead of a RATE.  Measured 2026-09-07:
# at a 2-nat floor the 5% rate admits 1,313-47,898 expected false CANDIDATES per
# arm, and the FDR is computed on candidates while the reported quantity is
# events after the overlap collapse -- a gap that did not matter at 0.01% FDR
# and is fatal at 5%.  An absolute budget bounds the false EVENT count too,
# since dedup can only ever reduce it.
BUDGET = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else None
tag = ("" if MAXD == 4 else f"_d{MAXD}") + (f"_lam{LAM0:g}" if LAM0 != 10.0 else "") \
      + (f"_b{BUDGET:g}" if BUDGET is not None else "")

# fixed threshold grid -- see 40_event_catalogue.py for why it cannot depend on
# the data.  Identical definition, so the two catalogues are directly comparable.
GRID = np.unique(np.round(np.geomspace(LAM0, 1.0e7, 600), 4))


# clade-size strata -- see 40_event_catalogue.py for why these are required and
# not a refinement: a fixed Lambda floor is a clade-size filter in disguise,
# since a fully-missing m-cell clade at expected rate p caps at ~m*log(1/p) nats.
SIZE_EDGES = np.array([4, 6, 10, 20, 50, 200, 10 ** 9])
NB = len(SIZE_EDGES) - 1
SIZE_LABEL = ["4-5", "6-9", "10-19", "20-49", "50-199", "200+"]


def strat_of(sizes):
    return np.clip(np.searchsorted(SIZE_EDGES, sizes, side="right") - 1, 0, NB - 1)


def hist_add(acc, vals, sizes):
    """Accumulate #{Lambda in bin} per size stratum, streaming."""
    bi = np.clip(np.searchsorted(GRID, vals, side="right") - 1, 0, GRID.size - 1)
    acc += np.bincount(strat_of(sizes) * GRID.size + bi,
                       minlength=NB * GRID.size).reshape(NB, GRID.size)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s], clone[s]
cache = RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4 else f"prefix_codes_{arm}.npz")
codes = np.load(cache, allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
miss = ~Y
print(f"{arm}: {n:,} cells, {K} tapes, {Gc:,} clones, depths 1-{MAXD} "
      f"({cache.name})", flush=True)


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


def scan(mi, a1, a2, pp, counts_only=False):
    """Return (lam_soft, lam_hard[eps], pi_hat, exp_rate, size, depth, counts).

    ⚠ counts_only streams: it accumulates the stratified histogram and keeps NO
    per-candidate arrays.  At a 4-nat floor Pre-TX yields millions of candidates
    per scan and the permutation parts would otherwise hold nine float arrays of
    that length, 1,000 times over."""
    acc = np.zeros((NB, GRID.size), dtype=np.int64)
    LS, LH, PI, ER, SZ, DP = [], {e: [] for e in EPS_SWEEP}, [], [], [], []
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
                hist_add(acc, ls[gg, zz], mm[gg, 0])
            if gg.size and not counts_only:
                LS.append(ls[gg, zz]); PI.append(pi[gg, zz])
                ER.append(er[gg, zz]); SZ.append(mm[gg, 0])
                DP.append(np.full(gg.size, d, dtype=float))
                for e_ in EPS_SWEEP:
                    LH[e_].append(k[gg, zz] * np.log(1 - e_) +
                                  (mm[gg, 0] - k[gg, zz]) * np.log(e_) - (A + B)[gg, zz])
    cat = lambda L: np.concatenate(L) if L else np.zeros(0)
    ge = np.cumsum(acc[:, ::-1], axis=1)[:, ::-1]        # #{Lambda >= t} per stratum
    return (cat(LS), {e: cat(v) for e, v in LH.items()}, cat(PI), cat(ER),
            cat(SZ), cat(DP), ge)


order = np.argsort(g_clone, kind="stable")
bnds = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])


def perm_inv(b):
    """Within-clone row permutation for permutation index b, seeded from (SEED, b)."""
    rb = np.random.default_rng([SEED, b])
    pm = order.copy()
    for a_, b_ in zip(bnds[:-1], bnds[1:]):
        pm[a_:b_] = rb.permutation(pm[a_:b_])
    inv = np.empty_like(pm)
    inv[order] = pm
    return inv


if PART is not None:
    ip, npart = (int(x) for x in PART.split("/"))
    assert 1 <= ip <= npart, PART
    lo, hi = ((ip - 1) * NPERM) // npart, (ip * NPERM) // npart
    idx = list(range(lo, hi))
    print(f"  part {ip}/{npart}: permutations {lo}..{hi-1} of B={NPERM}", flush=True)
    rows, t0 = [], time.time()
    for j, b in enumerate(idx):
        inv = perm_inv(b)
        *_, cn = scan(MISS[inv], ARR1[inv], ARR2[inv], P[inv], counts_only=True)
        rows.append(cn)
        print(f"    perm {b}: {cn[:,0].sum():,} soft candidates "
              f"[{j+1}/{len(idx)}, {time.time()-t0:.0f}s]", flush=True)
    outd = RES / "permparts"
    outd.mkdir(exist_ok=True)
    fn = outd / f"permcounts_43_{arm}{tag}_p{ip:03d}.npz"
    np.savez_compressed(fn, grid=GRID, counts=np.array(rows, dtype=np.int64),
                        size_edges=SIZE_EDGES,
                        perm_index=np.array(idx, dtype=np.int64),
                        meta=np.array([NPERM, npart, ip, SEED, LAM0, MAXD], dtype=float))
    print(f"wrote {fn.relative_to(ROOT)}  ({len(idx)} permutations, "
          f"{time.time()-t0:.0f}s)")
    sys.exit(0)

ls, lh, pi, er, sz, dp, obs_c = scan(MISS, ARR1, ARR2, P)
obs_c = obs_c.astype(float)
grid = GRID
print(f"  soft candidates at Lambda >= {LAM0}: {ls.size:,}  (by stratum " +
      " ".join(f"{SIZE_LABEL[i]}:{obs_c[i,0]:,.0f}" for i in range(NB)) + ")", flush=True)

if NULLF is not None:
    zn = np.load(NULLF, allow_pickle=False)
    assert np.array_equal(zn["grid"], GRID), "null file was built on a different grid"
    assert np.array_equal(zn["size_edges"], SIZE_EDGES), "null file used different strata"
    null_counts = zn["counts"].astype(float)             # (B, NB, G)
    print(f"  null: {null_counts.shape[0]:,} pooled permutations from "
          f"{Path(NULLF).name}", flush=True)
else:
    rows = []
    for b in range(NPERM):
        inv = perm_inv(b)
        *_, cn = scan(MISS[inv], ARR1[inv], ARR2[inv], P[inv], counts_only=True)
        rows.append(cn)
        print(f"    perm {b+1}/{NPERM}: {cn[:,0].sum():,.0f} soft candidates", flush=True)
    null_counts = np.array(rows, dtype=float)
B = int(null_counts.shape[0])
assert null_counts.shape[1:] == (NB, grid.size), null_counts.shape
nul_c = null_counts.mean(0)


def qcurve(obs, nul):
    """BH-style monotonisation: running minimum UP the grid (see 40 for why)."""
    raw = np.where(obs > 0, nul / np.maximum(obs, 1.0), 1.0)
    return np.minimum(np.minimum.accumulate(raw), 1.0)


QS = np.array([qcurve(obs_c[i], nul_c[i]) for i in range(NB)])
js = np.full(NB, -1, dtype=int)
for i in range(NB):
    okc = (np.flatnonzero(nul_c[i] <= BUDGET) if BUDGET is not None
           else np.flatnonzero(QS[i] <= TARGET))
    if okc.size:
        js[i] = int(okc[0])
LAMv = np.array([grid[js[i]] if js[i] >= 0 else np.inf for i in range(NB)])
print(f"\n  per-stratum soft threshold at FDR <= {100*TARGET:g}%:")
print(f"   {'clade':>8} {'candidates':>12} {'threshold':>10} {'FDR%':>8} {'null mean':>10}")
for i in range(NB):
    t_ = f"{LAMv[i]:.1f}" if js[i] >= 0 else "NONE"
    f_ = f"{100*QS[i,js[i]]:.4f}" if js[i] >= 0 else "--"
    n_ = nul_c[i, js[i]] if js[i] >= 0 else nul_c[i, 0]
    print(f"   {SIZE_LABEL[i]:>8} {obs_c[i,0]:>12,.0f} {t_:>10} {f_:>8} {n_:>10.2f}")

on, nn = obs_c.sum(0), nul_c.sum(0)
fdrc = qcurve(on, nn)
okg = np.flatnonzero(fdrc <= TARGET)
j0 = int(okg[0]) if okg.size else len(grid) - 1
LAM = float(grid[j0]); FDR = float(fdrc[j0])
use = [i for i in range(NB) if js[i] >= 0]
Cobs = float(sum(obs_c[i, js[i]] for i in use))
Cnull = (sum(null_counts[:, i, js[i]] for i in use) if use else np.zeros(B))
p_thresh = (1 + int((Cnull >= Cobs).sum())) / (B + 1)
p_floor = (1 + int((null_counts[:, :, 0].sum(1) >= on[0]).sum())) / (B + 1)
show = np.unique(np.clip(np.searchsorted(
    grid, [LAM0, 6, 10, 15, 20, 30, 50, 100, 300]), 0, len(grid) - 1))
print(f"\n  pooled FDR curve: " + "  ".join(
    f"{grid[i]:.0f}n:{100*fdrc[i]:.3g}%" for i in show))
print(f"  => pooled soft threshold {LAM:.1f} nats at FDR {100*FDR:.3g}%")
print(f"  global permutation p (B={B:,}): {p_thresh:.4g} at the stratum thresholds "
      f"({Cobs:,.0f} obs vs null max {Cnull.max():,.0f}), {p_floor:.4g} at the "
      f"{LAM0:g}-nat floor", flush=True)

sel = ls >= LAMv[strat_of(sz)]           # each clade against its own size stratum
out = {"arm": arm, "max_depth": MAXD, "lambda_soft_threshold": LAM, "fdr": FDR,
       "n_perm": B, "seed": SEED,
       "null_source": (Path(NULLF).name if NULLF else "inline"),
       "n_soft_candidates": int(sel.sum()), "eps_sweep": {}, "by_depth": {},
       "p_global_at_threshold": p_thresh, "p_global_at_scan_floor": p_floor,
       "p_floor_resolution": 1.0 / (B + 1),
       "null_mean_at_threshold": float(nn[j0]),
       "null_max_at_threshold": float(null_counts.sum(1)[:, j0].max()),
       "null_mean_at_scan_floor": float(nn[0]),
       "null_max_at_scan_floor": float(null_counts.sum(1)[:, 0].max()),
       "target_fdr": TARGET, "size_labels": SIZE_LABEL,
       "size_edges": SIZE_EDGES.tolist(),
       "per_stratum": {SIZE_LABEL[i]: {
           "candidates": float(obs_c[i, 0]),
           "threshold": (float(LAMv[i]) if js[i] >= 0 else None),
           "fdr": (float(QS[i, js[i]]) if js[i] >= 0 else None)} for i in range(NB)},
       "n_selected_at_stratum_thresholds": int(sel.sum()),
       "fdr_curve": {"lambda": grid.tolist(), "observed": on.tolist(),
                     "null": nn.tolist(), "fdr": fdrc.tolist()}}
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
    print("\n  pi_hat by the DEPTH of the clade -- if 'partial' events are just clades")
    print("  scored too coarsely, pi_hat must rise with depth:")
    print("   depth   events   median pi   >=0.99    <0.90   median size")
    for d_ in DEPTHS:
        m_ = sel & (dp == d_)
        if m_.sum() < 20:
            continue
        p_d = pi[m_]
        out["by_depth"][str(d_)] = {
            "n": int(m_.sum()), "median_pi": float(np.median(p_d)),
            "frac_ge_099": float((p_d >= 0.99).mean()),
            "frac_lt_090": float((p_d < 0.90).mean()),
            "median_size": float(np.median(sz[m_])),
            "median_expected": float(np.median(er[m_]))}
        r_ = out["by_depth"][str(d_)]
        print(f"   {d_:>5} {r_['n']:>8,} {r_['median_pi']:>11.3f} "
              f"{100*r_['frac_ge_099']:>8.1f}% {100*r_['frac_lt_090']:>7.1f}% "
              f"{r_['median_size']:>13.0f}", flush=True)
    print("\n  eps sensitivity -- share of soft events the HARD test also keeps:")
    for e_ in EPS_SWEEP:
        d = out["eps_sweep"][str(e_)]
        print(f"    eps={e_:<6} {d['n_hard_within_soft']:>8,}  "
              f"({100*d['frac_of_soft']:.1f}% of soft)", flush=True)
(RES / f"soft_events_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/soft_events_{arm}{tag}.json")
