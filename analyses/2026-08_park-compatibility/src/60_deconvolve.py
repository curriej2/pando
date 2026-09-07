#!/usr/bin/env python3
r"""
================================================================================
 60 -- DECONVOLVING tape -> symbol from many clades at once
================================================================================

⚠⚠ WHY 59 CANNOT ANSWER THIS ON ITS OWN.  59's statistic depends on the CLADE,
not on the tape: it contrasts the clade's insert composition with the rest of its
clone, and the lost tape enters only by being excluded from the pool.  So a clade
that lost k integrations yields ONE depletion vector carrying all k depletions
superposed.  Measured 2026-09-07: Subclone clone 6's depth-1 clade shows a LADDER
of six depleted symbols (Lambda 85, 36, 24, 21, 13, 9), and its three catalogue
"events" are one clade scored three times -- identical theta, identical Lambda.
Reading a map off single clades is therefore impossible in principle, not merely
underpowered.

--------------------------------------------------------------------------------
 THE MODEL
--------------------------------------------------------------------------------
Unit g is a scored group (a clade, or a whole clone for the clone-wide layer).
X[g, z] = 1 if tape z is called lost in g.  For each symbol s,

    y[g, s] = log-odds shift of s inside g vs outside   ~   SUM_z X[g,z] * B[z,s]

fitted by WEIGHTED least squares, weight = the inverse variance of y[g,s].  Under
the co-integration hypothesis B[z, .] is nonzero at exactly one symbol s(z), so
the recovered argmax_s (-B[z,s]) IS the map -- and the map, not any single
coefficient, is the evidence.

Thousands of units with different, overlapping loss sets is what identifies B.
⚠ Losses CLUSTER (rho_tape = 0.25), so tapes that are always lost together are
not separable.  This is a property of the design matrix, not of the method, so
the script reports the conditioning honestly: effective rank of X, and a per-tape
variance inflation factor.  A tape with a high VIF gets an unreliable call and is
excluded from the concordance test rather than quietly averaged in.

--------------------------------------------------------------------------------
 y AND ITS WEIGHT: site-stratified Mantel-Haenszel
--------------------------------------------------------------------------------
Composition varies by site (site 6's TV vs pooled is ~2x sites 1-5), and clade
writes sit at deeper sites than clone writes, so every contrast is stratified by
site.  59 profiles a full likelihood per symbol; that is too slow for thousands of
units, so here the same estimand is taken in closed form:

    OR_MH = SUM_j R_j / SUM_j S_j,   R_j = a_j (B_j - b_j) / n_j
                                     S_j = b_j (A_j - a_j) / n_j

with the Robins-Breslow-Greenland variance.  y = log OR_MH, weight = 1 / var.
⚠ Both estimate the same site-stratified common odds ratio; MH is used for the
scan and the likelihood fit is kept for the hand-inspected examples in 59, so the
two should agree there -- that agreement is worth checking, not assuming.

--------------------------------------------------------------------------------
 LAYERS
--------------------------------------------------------------------------------
 clade      units are the distinct (clone, anchor, depth) clades of the event
            catalogue; contrast is against the rest of the clone.  The lineage
            claim, and the internal control is tight.
 clonewide  units are whole clones with a clone-wide called loss; contrast is
            against all OTHER clones in the arm.  ⚠ Weaker control (clones differ
            by founder effect) but far more units, the loss predates the founder
            so every within-clone write is post-loss, and it is the only layer
            with power in Pre-TX, whose clades run 71-95 cells and returned
            nothing above the null in 4 of 5 cases.

⚠⚠ THE COMPARISON GROUP MUST BE POOLED FROM PER-CLONE EXTRACTIONS, NOT EXTRACTED
AS ONE GROUP.  Handing extract_writes every cell outside clone c makes it assess
polymorphism ACROSS clones, where almost every parent has several children, so
the identity-by-descent dedup collapses and nearly every entry counts as a write
-- the exact pseudoreplication this whole direction was rebuilt to avoid.  So
writes are extracted once per clone and the contrast is assembled by arithmetic:

    inside(c)  = T_c        - SUM_{z in drop} T_cz
    outside(c) = T - T_c    - SUM_{z in drop} (T_z - T_cz)

with T the arm total, T_c the clone's table, T_z a tape's table over all clones,
and T_cz the (clone, tape) cross term, kept only for the tapes actually dropped.

--------------------------------------------------------------------------------
 ⚠⚠ --split: THE ONLY CLEAN REPRODUCIBILITY TEST, AND WHY CROSS-ARM IS NOT ONE
--------------------------------------------------------------------------------
The README treats five-way cross-arm concordance as the killer test, on the
grounds that all 166 TargetBCs are identical across arms.  They are -- but the
ARMS ARE NOT INDEPENDENT SAMPLES.  Measured 2026-09-07 on ClonalBC:

    Initial x Mouse1  285 shared clones = 96.3% of Mouse 1's
    Initial x Mouse2  210 = 96.3%      Initial x Mouse3  143 = 95.3%
    Mouse1  x Mouse2   93 = 42.7%      Mouse2  x Mouse3   47 = 31.3%

The mice were transplanted from the pre-TX pool, silencing is HERITABLE and
predates transplantation, so the same clone measured in two arms carries the same
losses.  Cross-arm agreement can therefore arise with no mechanism at all -- it is
the same events re-measured, not a replication.  (Only Subclone is near-disjoint
from the mice: 1, 4 and 1 shared clones.)

⇒ --split h/N fits the map on a random 1/N of the arm's CLONES.  Two disjoint
splits share no cell and no clone and no lineage, so agreement between them tests
the one thing that matters: whether z -> s(z) is a property of the LINE rather
than of the clones it was estimated from.  Same seed, different h, disjoint sets.

Usage: 60_deconvolve.py <arm> [--layer clade|clonewide|both] [--ridge 1.0]
                        [--poly codes|codes_or_none] [--maxunits N]
                        [--split h/N] [--splitseed S]
Outputs results/depmap_{arm}_{layer}.npz   (B, se, VIF, X diagnostics)
        results/depmap_{arm}_{layer}.json  (the recovered map + conditioning)
================================================================================
"""
import gzip, json, sys, time, collections
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_w = __import__("57_writes")
load_arm, extract_writes, clade_of = _w.load_arm, _w.extract_writes, _w.clade_of

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
NSITE = 6

arm = sys.argv[1]
LAYER = sys.argv[sys.argv.index("--layer") + 1] if "--layer" in sys.argv else "clade"
RIDGE = float(sys.argv[sys.argv.index("--ridge") + 1]) if "--ridge" in sys.argv else 1.0
POLY = sys.argv[sys.argv.index("--poly") + 1] if "--poly" in sys.argv else "codes"
MAXU = int(sys.argv[sys.argv.index("--maxunits") + 1]) if "--maxunits" in sys.argv else 10 ** 9
SPLIT = sys.argv[sys.argv.index("--split") + 1] if "--split" in sys.argv else None
SPLITSEED = int(sys.argv[sys.argv.index("--splitseed") + 1]) if "--splitseed" in sys.argv else 7
tag = LAYER + (f"_s{SPLIT.replace('/', 'of')}" if SPLIT else "")

Y, clone, codes, S, syms0, tapes = load_arm(arm)
cmask, remap, syms = _w.conforming(syms0)
K, A = len(tapes), len(syms)
cvals = _w.clone_values(clone)
print(f"{arm}: {len(clone):,} cells, {K} tapes, {A} conforming symbols, layer={LAYER}",
      flush=True)


def mh(a, b):
    """Site-stratified Mantel-Haenszel log-OR and its RBG variance.

    a, b : (J, A) inside / outside counts of each symbol.  Returns (y, var).
    """
    A_j, B_j = a.sum(1, keepdims=True), b.sum(1, keepdims=True)
    a1, a0 = a, A_j - a
    b1, b0 = b, B_j - b
    n = A_j + B_j
    R, Sm = a1 * b0 / n, b1 * a0 / n
    Rs, Ss = R.sum(0), Sm.sum(0)
    ok = (Rs > 0) & (Ss > 0)
    y = np.where(ok, np.log(np.maximum(Rs, 1e-12) / np.maximum(Ss, 1e-12)), 0.0)
    P, Q = (a1 + b0) / n, (b1 + a0) / n
    v = np.where(ok, (P * R).sum(0) / (2 * np.maximum(Rs, 1e-12) ** 2)
                 + ((P * Sm).sum(0) + (Q * R).sum(0))
                 / (2 * np.maximum(Rs * Ss, 1e-12))
                 + (Q * Sm).sum(0) / (2 * np.maximum(Ss, 1e-12) ** 2), np.inf)
    return y, np.where(ok, v, np.inf)


def tab(t_, j_, s_):
    return np.bincount(j_ * A + s_, minlength=NSITE * A).reshape(NSITE, A).astype(float)


def score_unit(ins, out, drop):
    ti, ji, si = _w.filter_writes(*extract_writes(codes[ins], S[ins], Y[ins],
                                                  drop_tapes=drop, poly=POLY), remap)
    to, jo, so = _w.filter_writes(*extract_writes(codes[out], S[out], Y[out],
                                                  drop_tapes=drop, poly=POLY), remap)
    if len(ti) < 100 or len(to) < 100:
        return None
    y, v = mh(tab(ti, ji, si), tab(to, jo, so))
    return y, v, len(ti), len(to)


# ---------------------------------------------------------------------------
#  PER-CLONE WRITE TABLES -- the clone-wide layer's contrast is built from these
# ---------------------------------------------------------------------------
CW = {}
if LAYER in ("clonewide", "both"):
    T_all = np.zeros((NSITE, A))
    T_tape = np.zeros((K, NSITE, A))
    T_cl, T_cl_tape = {}, {}
    _cw_keep = None
    if SPLIT:
        _h, _N = (int(x) for x in SPLIT.split("/"))
        _allc = np.unique(clone)
        _lab = np.random.default_rng(SPLITSEED).permutation(len(_allc)) % _N
        _cw_keep = set(_allc[_lab == _h].tolist())
    for cl in np.unique(clone):
        if _cw_keep is not None and int(cl) not in _cw_keep:
            continue
        g = clone == cl
        if g.sum() < 2:
            continue
        ti, ji, si = _w.filter_writes(*extract_writes(codes[g], S[g], Y[g], poly=POLY),
                                      remap)
        if len(ti) == 0:
            continue
        tc = tab(ti, ji, si)
        T_cl[int(cl)] = tc
        T_all += tc
        np.add.at(T_tape, (ti, ji, si), 1.0)
        T_cl_tape[int(cl)] = (ti, ji, si)
    print(f"  clone-wide layer: {len(T_cl):,} clones, "
          f"{int(T_all.sum()):,} writes pooled", flush=True)


def score_clonewide(cl, drop):
    tc = T_cl.get(int(cl))
    if tc is None:
        return None
    ti, ji, si = T_cl_tape[int(cl)]
    dm = np.isin(ti, drop)
    cross = tab(ti[dm], ji[dm], si[dm])          # T_cz, summed over dropped tapes
    ins = tc - cross
    out = T_all - tc - (T_tape[list(drop)].sum(0) - cross)
    if ins.sum() < 100 or out.sum() < 100:
        return None
    y, v = mh(ins, out)
    return y, v, int(ins.sum()), int(out.sum())


# ---------------------------------------------------------------------------
#  UNITS
# ---------------------------------------------------------------------------
units = []
if LAYER in ("clade", "both"):
    lost = collections.defaultdict(set)
    for r in _w._events(arm):
        lost[(int(r["clone"]), int(r["anchor_tape"]), int(r["depth"]))].add(int(r["tape"]))
    for (ci, an, d), tp in lost.items():
        units.append(("clade", int(cvals[ci]), an, d, sorted(tp)))
if LAYER in ("clonewide", "both"):
    cw = collections.defaultdict(set)
    with gzip.open(RES / f"clonewide_{arm}.tsv.gz", "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            r = dict(zip(hdr, ln.rstrip("\n").split("\t")))
            cw[int(r["clone"])].add(int(r["tape"]))
    for ci, tp in cw.items():
        units.append(("clonewide", int(cvals[ci]), -1, -1, sorted(tp)))
if SPLIT:
    h, N = (int(x) for x in SPLIT.split("/"))
    allc = np.unique(clone)
    lab = np.random.default_rng(SPLITSEED).permutation(len(allc)) % N
    keep = set(allc[lab == h].tolist())
    before = len(units)
    units = [u for u in units if u[1] in keep]
    print(f"  split {h}/{N} (seed {SPLITSEED}): {len(keep):,}/{len(allc):,} clones, "
          f"{len(units):,}/{before:,} units", flush=True)
units = units[:MAXU]
print(f"  {len(units):,} units "
      f"({sum(1 for u in units if u[0]=='clade')} clade, "
      f"{sum(1 for u in units if u[0]=='clonewide')} clone-wide)", flush=True)

rowsY, rowsV, rowsX, meta = [], [], [], []
t0 = time.time()
for i, (kind, cl, an, d, tp) in enumerate(units):
    drop = tuple(tp) + ((an,) if an >= 0 else ())
    if kind == "clade":
        ins, out = clade_of(codes, clone, cl, an, d)
        # ⚠ the catalogue's clade for this (clone, anchor, depth) is the largest
        # code group; clade_of defaults to it, which is what the scan scored.
        r = score_unit(ins, out, drop)
        if r is None:
            continue
        y, v, wi, wo = r
    else:
        r = score_clonewide(cl, tp)
        if r is None:
            continue
        y, v, wi, wo = r
    x = np.zeros(K); x[list(tp)] = 1.0
    rowsY.append(y); rowsV.append(v); rowsX.append(x)
    meta.append({"kind": kind, "clone": cl, "anchor": an, "depth": d,
                 "n_lost": len(tp), "W_in": wi, "W_out": wo})
    if (i + 1) % 200 == 0:
        print(f"    {i+1}/{len(units)} units  [{time.time()-t0:.0f}s]", flush=True)

Ym, Vm, Xfull = np.array(rowsY), np.array(rowsV), np.array(rowsX)
G = Ym.shape[0]
# ⚠ Only tapes that are actually LOST somewhere are estimable, and a tape seen in
# one or two units is estimable in name only.  Fitting the full 166 columns on an
# arm with fewer units than tapes (Mouse 3: 155 units) is underdetermined, and
# the ridge would hide that rather than report it.  So restrict the design to
# covered tapes and mark the rest uncallable, explicitly.
MIN_COV = 3
covfull = Xfull.sum(0)
cols = np.where(covfull >= MIN_COV)[0]
assert len(cols) > 0, f"{arm}: no tape lost in >= {MIN_COV} units"
if len(cols) >= G:
    cols = cols[np.argsort(-covfull[cols])[:max(G // 2, 1)]]
Xm = Xfull[:, cols]
Kf = len(cols)
print(f"  {G:,} usable units; {Kf}/{K} tapes with coverage >= {MIN_COV} are fitted; "
      f"X density {Xm.mean():.4f}", flush=True)

# ---- conditioning, reported not assumed -----------------------------------
sv = np.linalg.svd(Xm, compute_uv=False)
eff_rank = int((sv > sv.max() * 1e-10).sum())
XtX = Xm.T @ Xm + RIDGE * np.eye(Kf)
XtXi = np.linalg.inv(XtX)
viff = np.diag(XtXi) * np.diag(XtX)         # variance inflation per fitted tape
vif = np.full(K, np.inf); vif[cols] = viff
print(f"  design: rank {eff_rank}/{Kf}, cond {sv.max()/max(sv.min(),1e-12):.1f}, "
      f"VIF median {np.median(viff):.2f} max {viff.max():.1f}", flush=True)

# ---- weighted least squares, one symbol at a time -------------------------
Wm = np.where(np.isfinite(Vm) & (Vm > 0), 1.0 / np.maximum(Vm, 1e-6), 0.0)
B = np.zeros((K, A)); SE = np.zeros((K, A))
for s in range(A):
    w = Wm[:, s]
    Xw = Xm * w[:, None]
    M = Xm.T @ Xw + RIDGE * np.eye(Kf)
    Mi = np.linalg.inv(M)
    B[cols, s] = Mi @ (Xw.T @ Ym[:, s])
    SE[cols, s] = np.sqrt(np.maximum(np.diag(Mi), 0.0))
Z = np.where(SE > 0, B / np.maximum(SE, 1e-12), 0.0)

# ---- the recovered map ----------------------------------------------------
best = np.argmin(Z, axis=1)                  # most depleted symbol per tape
zbest = Z[np.arange(K), best]
rest = Z.copy(); rest[np.arange(K), best] = 0.0
second = np.min(rest, axis=1)
gap = second - zbest                         # >0 means the top is a clear outlier
cov = covfull
good = np.zeros(K, dtype=bool); good[cols] = True
good &= (cov >= 5) & (vif < 5)
uniq = len(set(best[good].tolist()))
res = {
    "arm": arm, "layer": LAYER, "tag": tag, "split": SPLIT,
    "poly": POLY, "ridge": RIDGE,
    "n_units": int(G), "n_tapes": int(K), "n_tapes_fitted": int(Kf),
    "min_coverage": MIN_COV, "n_symbols": int(A),
    "design": {"rank": eff_rank, "cond": float(sv.max() / max(sv.min(), 1e-12)),
               "vif_median": float(np.median(vif)), "vif_max": float(vif.max())},
    "callable_tapes": int(good.sum()),
    "distinct_symbols_called": uniq,
    "expected_distinct_166_of_256": 122.0,
    "z_best_median": float(np.median(zbest[good])) if good.any() else None,
    "gap_median": float(np.median(gap[good])) if good.any() else None,
    "map": [{"tape": str(tapes[t]), "idx": int(t), "symbol": str(syms[best[t]]),
             "z": float(zbest[t]), "gap": float(gap[t]), "beta": float(B[t, best[t]]),
             "n_units": int(cov[t]), "vif": float(vif[t]), "callable": bool(good[t])}
            for t in range(K)],
}
np.savez_compressed(RES / f"depmap_{arm}_{tag}.npz", B=B, SE=SE, Z=Z, X=Xfull,
                    cols=cols, vif=vif, cov=cov, symbols=syms, tapes=tapes,
                    good=good)
(RES / f"depmap_{arm}_{tag}.json").write_text(json.dumps(res, indent=1))
print(f"{arm}/{LAYER}: {good.sum()}/{K} tapes callable, {uniq} distinct symbols "
      f"(122 predicted for 166 draws from 256); median z {res['z_best_median']}, "
      f"median gap {res['gap_median']}", flush=True)
print("  strongest 10: " + "  ".join(
    f"{tapes[t]}->{syms[best[t]]}(z={zbest[t]:.1f})"
    for t in np.where(good)[0][np.argsort(zbest[good])[:10]]), flush=True)
