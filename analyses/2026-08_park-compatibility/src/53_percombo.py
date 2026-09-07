#!/usr/bin/env python3
r"""
================================================================================
 Per-combo nulls: where does the real labelling fall, IN NATS, for each combo?
================================================================================

Justin's question (2026-09-07), and it is the right one: rather than pooling
candidates and thresholding a global Lambda, score EACH (clade, tape) combo
against its OWN 1,000 permutations and report the margin in nats.

WHY THIS IS BETTER THAN THE STRATIFIED GLOBAL THRESHOLD IT REPLACES

 1. It conditions on everything, automatically.  The whole clade-size-strata
    apparatus in 40_event_catalogue.py exists to stop a 5-cell clade being
    judged against a null dominated by 500-cell clades.  A per-combo null
    conditions on clade size, on that clade's own p_tilde values, and on its
    clone's composition -- at full resolution, not in six hand-drawn bins.
 2. The Lambda CAP stops being a problem.  "A 9-cell Pre-TX clade cannot exceed
    12.4 nats" only bit because 12.4 was compared to a global threshold.  The cap
    applies to that clade's null too: a combo whose own permutations top out at
    6 nats and which scores 11 is decisive, whatever any global cut says.
 3. No candidate-vs-event FDR mismatch.  Each combo is scored on its own terms
    instead of picking a global cut and then deduplicating.

⚠ WHY IT WAS NOT DONE FIRST.  The permutation parts stored count vectors
AGGREGATED OVER COMBOS, throwing combo identity away.  That was right for the
aggregate FDR and wrong for this.  A combo does have a stable identity --
(depth, anchor, subclade code, tape), all derived from `codes`, which the
permutation never touches -- so per-combo accumulators are well defined.

⚠ REPORT THE MARGIN IN NATS, NOT A Z-SCORE.  For a small clone the null is
genuinely coarse: a 4-cell clade in a 6-cell clone has only C(6,4) = 15 distinct
compositions, so 1,000 permutations sample 15 values and sd can be exactly 0.
A z-score explodes there.  margin = Lambda_obs - max_b Lambda_b degrades
gracefully and is reported as the primary statistic; mean/sd are context.

⚠ MULTIPLICITY STILL NEEDS AN AGGREGATE, and this is how it is calibrated.
With up to 105M combos some will beat their own null by chance, so the margin
needs its own null distribution ACROSS combos -- which accumulators alone cannot
give, since the margin is a function of the whole permutation ensemble.

THE FIX: HOLD OUT the first HOLD permutations as PSEUDO-OBSERVED.  Accumulators
run over the remaining B - HOLD only.  Then for every combo

    m_obs  = Lambda_obs      - max_{b >= HOLD} Lambda_b
    m_null = Lambda_{b<HOLD} - max_{b >= HOLD} Lambda_b        (HOLD draws)

are margins of a candidate over the max of the SAME permutation set, so under
H0 they are exactly exchangeable -- no approximation, no parametric null, and
nothing censored at 1/(B+1).  HOLD x NC null draws then calibrate the aggregate
FDR on the margin: FDR(t) = (#m_null >= t)/HOLD / (#m_obs >= t).  That is script
54; this script produces its inputs.

WHAT IS ACCUMULATED, per combo, streaming (never 1,000 values per combo):
    n_ge   #{b : Lambda_b >= Lambda_obs}      -> exceedance count
    mx     max_b Lambda_b                     -> the margin, Lambda_obs - mx
    s1     sum_b Lambda_b                     -> mean
    s2     sum_b Lambda_b^2                   -> sd
    pseudo Lambda for the HOLD held-out permutations -> the margin's own null

Usage: 53_percombo.py <arm> [--nperm 1000] [--permpart i/N] [--seed S]
                      [--eps 0.01] [--maxd 4]
Output: results/percombo/percombo_{arm}_p{i}.npz   (obs + the four accumulators)
================================================================================
"""
import sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4

arm = sys.argv[1]
EPS = float(sys.argv[sys.argv.index("--eps") + 1]) if "--eps" in sys.argv else 0.01
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 1000
SEED = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260903
PART = sys.argv[sys.argv.index("--permpart") + 1] if "--permpart" in sys.argv else "1/1"
MAXD = int(sys.argv[sys.argv.index("--maxd") + 1]) if "--maxd" in sys.argv else 4
DEPTHS = list(range(1, MAXD + 1))
HOLD = int(sys.argv[sys.argv.index("--hold") + 1]) if "--hold" in sys.argv else 3

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
cache = RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4 else f"prefix_codes_{arm}.npz")
codes = np.load(cache, allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
miss = ~Y
print(f"{arm}: {n:,} cells, {K} tapes, {Gc:,} clones, depths 1-{MAXD}", flush=True)


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
W = np.where(miss, np.log(1 - EPS) - np.log(P), np.log(EPS) - np.log1p(-P))

# ---- block layout, built ONCE.  `codes` and `g_clone` are never permuted, so
# every block's clade partition -- and therefore each combo's global index -- is
# identical in the observed data and in every permutation.
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
        # ⚠ COMPACT the slot space: keep only clades that can be scored at all.
        # Allocating G*K per block wasted 79% of the array on Mouse3 (26.2M slots
        # for 5.5M scorable combos) because most clades fall below MIN_CLADE --
        # and on Pre-TX that waste would have been ~490M slots, ~24 GB.
        keep = np.flatnonzero(ns >= MIN_CLADE).astype(np.int32)
        if keep.size == 0:
            continue
        valid = np.ones((keep.size, K), bool)
        valid[:, a_] = False                     # never score the anchor itself
        blocks.append((a_, idx.astype(np.int32), sub.astype(np.int32), G, off,
                       valid.ravel(), keep, ns[keep].astype(np.int32)))
        off += keep.size * K
NC = off
print(f"  {len(blocks):,} blocks, {NC:,} combo slots "
      f"({100*sum(b[5].sum() for b in blocks)/max(NC,1):.1f}% scorable)", flush=True)


def lam_all(Wm, out):
    """Fill `out` (length NC) with Lambda for every kept combo slot."""
    for a_, idx, sub, G, o, _, keep, _ns in blocks:
        L = np.empty((G, K))
        Wi = Wm[idx]
        for k in range(K):
            L[:, k] = np.bincount(sub, weights=Wi[:, k], minlength=G)
        out[o:o + keep.size * K] = L[keep].ravel()


obs = np.empty(NC, dtype=np.float32)
lam_all(W, obs)

order = np.argsort(g_clone, kind="stable")
bnds = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])


def perm_inv(b):
    rb = np.random.default_rng([SEED, b])
    pm = order.copy()
    for a_, b_ in zip(bnds[:-1], bnds[1:]):
        pm[a_:b_] = rb.permutation(pm[a_:b_])
    inv = np.empty_like(pm); inv[order] = pm
    return inv


ip, npart = (int(x) for x in PART.split("/"))
assert 1 <= ip <= npart, PART
lo, hi = ((ip - 1) * NPERM) // npart, (ip * NPERM) // npart
# float32 throughout: Lambda is O(1-100), so 1,000 accumulated adds carry a
# relative error ~1e-4 on the mean -- irrelevant next to the coarseness of the
# null itself, and it halves a 5 GB footprint on Pre-TX.
n_ge = np.zeros(NC, dtype=np.int16)
mx = np.full(NC, -np.inf, dtype=np.float32)
s1 = np.zeros(NC, dtype=np.float32)
s2 = np.zeros(NC, dtype=np.float32)
cur = np.empty(NC, dtype=np.float32)
# permutations b < HOLD are PSEUDO-OBSERVED: kept aside, never accumulated, and
# used in script 54 as exchangeable null draws of the margin statistic.
pseudo, n_acc = {}, 0
print(f"  part {ip}/{npart}: permutations {lo}..{hi-1} of B={NPERM} "
      f"(hold-out: {[b for b in range(lo, hi) if b < HOLD]})", flush=True)
t0 = time.time()
for j, b in enumerate(range(lo, hi)):
    lam_all(W[perm_inv(b)], cur)
    if b < HOLD:
        pseudo[b] = cur.copy()
        print(f"    perm {b}: HELD OUT as pseudo-observed", flush=True)
        continue
    n_ge += (cur >= obs)
    np.maximum(mx, cur, out=mx)
    s1 += cur
    s2 += cur * cur
    n_acc += 1
    if n_acc % 25 == 0 or n_acc == 1:
        print(f"    perm {b} [{n_acc} accumulated, {time.time()-t0:.0f}s]", flush=True)

outd = RES / "percombo"; outd.mkdir(exist_ok=True)
fn = outd / f"percombo_{arm}{'' if MAXD == 4 else f'_d{MAXD}'}_p{ip:03d}.npz"
np.savez_compressed(
    fn, obs=obs, n_ge=n_ge,
    mx=np.where(np.isfinite(mx), mx, -1e30).astype(np.float32), s1=s1, s2=s2,
    valid=np.concatenate([b[5] for b in blocks]),
    size=np.concatenate([np.repeat(b[7], K) for b in blocks]).astype(np.int32),
    meta=np.array([NPERM, npart, ip, SEED, n_acc, MAXD, NC, HOLD], dtype=np.float64),
    hold_ids=np.array(sorted(pseudo), dtype=np.int64),
    **{f"pseudo{b}": pseudo[b] for b in sorted(pseudo)})
print(f"wrote {fn.relative_to(ROOT)}  ({n_acc} accumulated, "
      f"{len(pseudo)} held out, {time.time()-t0:.0f}s)")
