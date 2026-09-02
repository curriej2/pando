#!/usr/bin/env python3
r"""
================================================================================
 Is tape loss heritable BELOW the clone?  -- the relatedness gradient
================================================================================

Clone is a coarse grouping: Subclone's median clone is 997 cells and Mouse2's
signal is carried by one clone of 3,387, so a loss occurring partway down such a
tree is diluted by averaging over the whole clone.  Worse, colony identity in the
Subclone arm is also a culture batch, so clone-level agreement has a non-lineage
explanation.  Both problems have one resolution: look FINER than the clone.

  a Dollo loss on the tree  -> subclades agree MORE than the clone as a whole,
                               and more the deeper (closer) the subclade
  a flat well/batch effect  -> subclades agree exactly as much as random subsets
                               of the clone; no gradient

GROUPING, and why it is not circular.  A subclade is the set of cells sharing the
depth-d prefix of one ANCHOR tape z, within a clone.  Missingness is then measured
on every tape EXCEPT z.  Cross-tape only -- the same independence the compatibility
work rests on (D.4 fact 1).  Averaged over all 166 anchors, swept over d = 1..4.

THE NULL IS ANALYTIC.  Under random assignment of a clone's cells to subgroups of
the observed sizes, for any residual vector r,

    E[ sum_{c != c' in g} r_c r_c' ] = m(m-1)/(n(n-1)) * sum_{c != c' in C} r_c r_c'

because every unordered pair of the clone is equally likely to land together.  So
the permutation MEAN is a size-weighted rescaling of the clone-level pair sum, with
no simulation.  We take

    rho_sub  = SUM_z,C A / SUM_z,C D                    (observed, subclade groups)
    rho_null = SUM_z,C E[A] / SUM_z,C E[D]              (same cells, clone-level)
    excess   = rho_sub - rho_null                       (agreement BEYOND the clone)

--verify runs explicit within-clone permutations on a subsample of anchors and
checks them against the analytic null, since the whole result rests on it.

Usage: 35_lineage_depth.py <arm> [--screen] [--depths 1,2,3,4] [--verify N]
Output: results/lineage_depth_{arm}{_screen}.json, results/prefix_codes_{arm}.npz
================================================================================
"""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_COUNT, N_SITES, MAX_D = 1000, 6, 4
MIN_CLONE = 10
JUNK = 2

arm = sys.argv[1]
SCREEN = "--screen" in sys.argv
DEPTHS = ([int(x) for x in sys.argv[sys.argv.index("--depths") + 1].split(",")]
          if "--depths" in sys.argv else [1, 2, 3, 4])
VERIFY = int(sys.argv[sys.argv.index("--verify") + 1]) if "--verify" in sys.argv else 0
rng = np.random.default_rng(20260902)


# ---------------------------------------------------------------- prefix codes
def build_codes():
    """codes[c, z, d] = integer id of the depth-(d+1) prefix of tape z, -1 if the
    cell does not reach that depth.  Same cell set and row order as script 23's
    cache (ClonalBC present AND >= THR recovered tapes, in CSV order)."""
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    d0 = xi_all[arm]
    keep = {k for k, v in d0["xi"].items() if v * d0["n_edits"] >= MIN_COUNT}
    sym_id = {s: i for i, s in enumerate(sorted(keep))}
    clone_of = {}
    with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["ClonalBC"] and r["ClonalBC"] != "None":
                clone_of[r["CellID"]] = r["ClonalBC"]
    with open(DATA / f"{arm}_EditTable_filtered.csv", newline="") as fh:
        rr = csv.reader(fh); hdr = next(rr)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
        tapes = sorted(blocks)
        idxs = [[i for _, i in sorted(blocks[tp])] for tp in tapes]
        K = len(tapes)
        rows = []
        for row in rr:
            if clone_of.get(row[0]) is None:
                continue
            v = row[1:]
            a = np.full((K, MAX_D), -1, dtype=np.int16)
            ntapes = 0
            for t, cols in enumerate(idxs):
                seen = False
                for j, i in enumerate(cols):
                    val = v[i]
                    if val != "None":
                        seen = True
                    if j < MAX_D:
                        sid = sym_id.get(val, -1)
                        if sid < 0:
                            break
                        a[t, j] = sid
                    elif val not in keep:
                        break
                if seen:
                    ntapes += 1
            if ntapes < THR[arm]:
                continue
            rows.append(a)
    S = np.array(rows, dtype=np.int16)               # (n, K, MAX_D) symbol ids
    n, K, _ = S.shape
    A = len(sym_id) + 1
    codes = np.full((n, K, MAX_D), -1, dtype=np.int32)
    prev = np.zeros((n, K), dtype=np.int64)
    ok = np.ones((n, K), dtype=bool)
    for d in range(MAX_D):
        ok &= S[:, :, d] >= 0
        prev = prev * A + (S[:, :, d].astype(np.int64) + 1)
        flat = np.where(ok, prev, -1).ravel()
        uniq, inv = np.unique(flat, return_inverse=True)
        c = inv.reshape(n, K).astype(np.int32)
        codes[:, :, d] = np.where(ok, c, -1)
    return codes


cache = RES / f"prefix_codes_{arm}.npz"
if cache.exists():
    codes = np.load(cache, allow_pickle=False)["codes"]
else:
    codes = build_codes()
    np.savez_compressed(cache, codes=codes)
    print(f"built prefix codes -> {cache.name}", flush=True)

# ---------------------------------------------------------------- features
z = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, dep, trm, clone, sample = (z["recovered"], z["depth"].astype(np.int16),
                              z["term"], z["clone"], z["sample"])
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, dep, trm, clone = Y[sel], dep[sel], trm[sel], clone[sel]
assert codes.shape[0] == Y.shape[0], "prefix codes misaligned with the cell set"
if SCREEN:
    good = ~np.load(RES / f"collision_flags_{arm}.npz", allow_pickle=False)["flag"]
    Y, dep, trm, clone, codes = Y[good], dep[good], trm[good], clone[good], codes[good]
    print(f"collision screen: dropped {int((~good).sum()):,} cells")
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
print(f"{arm}{' +screen' if SCREEN else ''}: {n:,} cells, {K} tapes, {Gc:,} clones",
      flush=True)


def fit_twoway(M, W, iters=40):
    cm = (M * W).sum(0) / np.maximum(W.sum(0), 1)
    cm = np.clip(cm, 1e-3, 1 - 1e-3)
    beta = np.log(cm / (1 - cm)); alpha = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))
            if ax == 0:
                g = (W * (M - p)).sum(1); h = (W * p * (1 - p)).sum(1)
                alpha = np.clip(alpha + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
                alpha -= alpha.mean()
            else:
                g = (W * (M - p)).sum(0); h = (W * p * (1 - p)).sum(0)
                beta = np.clip(beta + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :]))), 1e-6, 1 - 1e-6)


families = {"missing": ((~Y).astype(float), np.ones_like(Y, dtype=float))}
for L in (5, 6):                       # T0 calibration, marginals nearest missing
    families[f"depth>={L}"] = ((dep >= L).astype(float),
                               (Y & ((dep >= L) | (trm != JUNK))).astype(float))
resid = {}
for name, (M, W) in families.items():
    P = fit_twoway(M, W)
    r = np.where(W > 0, M - P, 0.0).astype(np.float64)
    v = np.where(W > 0, P * (1 - P), 0.0).astype(np.float64)
    resid[name] = (r, r * r, np.sqrt(v), v)
    print(f"  {name:10s} marginal {float((M*W).sum()/W.sum()):.3f}", flush=True)


def pair_sums(x2, x1, lab, G, w=None):
    """per feature: SUM_g w_g [ (sum_g x1)^2 - sum_g x2 ]."""
    out = np.empty(x1.shape[1])
    for k in range(x1.shape[1]):
        S1 = np.bincount(lab, weights=x1[:, k], minlength=G)
        S2 = np.bincount(lab, weights=x2[:, k], minlength=G)
        t = S1 * S1 - S2
        out[k] = t.sum() if w is None else (t * w).sum()
    return out


out = {"arm": arm, "screen": SCREEN, "n_cells": int(n), "depths": {}}
for d in DEPTHS:
    acc = {nm: dict(A=0.0, D=0.0, EA=0.0, ED=0.0) for nm in families}
    n_anchor = 0
    for zi in range(K):
        cd = codes[:, zi, d - 1]
        valid = cd >= 0
        if valid.sum() < 2 * MIN_CLONE:
            continue
        cl = g_clone[valid]
        # keep only clones that are big enough AND actually split by this anchor
        key = cl.astype(np.int64) * (cd.max() + 2) + cd[valid]
        _, sub = np.unique(key, return_inverse=True)
        csz = np.bincount(cl); ssz = np.bincount(sub)
        nsubs = np.bincount(cl, weights=(np.bincount(sub, minlength=sub.max()+1) > 0)[sub] * 0.0)
        # a clone contributes only if it has >=2 subgroups and >=MIN_CLONE cells
        first = np.zeros(sub.max() + 1, dtype=np.int64)
        first[sub] = cl
        ngrp = np.bincount(first, minlength=csz.size)
        good_clone = (csz >= MIN_CLONE) & (ngrp >= 2)
        m = good_clone[cl]
        if m.sum() < 2 * MIN_CLONE:
            continue
        idx = np.flatnonzero(valid)[m]
        cl2 = g_clone[idx]
        _, cl2i = np.unique(cl2, return_inverse=True)
        _, sub2 = np.unique(key[m], return_inverse=True)
        Gs, Gk = sub2.max() + 1, cl2i.max() + 1
        # analytic within-clone permutation weight  c_C = SUM_g m(m-1) / n(n-1)
        ms = np.bincount(sub2).astype(float)
        owner = np.zeros(Gs, dtype=np.int64); owner[sub2] = cl2i
        num = np.bincount(owner, weights=ms * (ms - 1), minlength=Gk)
        ns = np.bincount(cl2i).astype(float)
        wC = num / np.maximum(ns * (ns - 1), 1e-9)
        n_anchor += 1
        for nm, (r, r2, sv, v) in resid.items():
            A = pair_sums(r2[idx], r[idx], sub2, Gs)
            D = pair_sums(v[idx], sv[idx], sub2, Gs)
            EA = pair_sums(r2[idx], r[idx], cl2i, Gk, w=wC)
            ED = pair_sums(v[idx], sv[idx], cl2i, Gk, w=wC)
            A[zi] = D[zi] = EA[zi] = ED[zi] = 0.0        # never the anchor itself
            acc[nm]["A"] += A.sum(); acc[nm]["D"] += D.sum()
            acc[nm]["EA"] += EA.sum(); acc[nm]["ED"] += ED.sum()
    rec = {"n_anchors": n_anchor}
    for nm, a in acc.items():
        rs = a["A"] / a["D"] if a["D"] else float("nan")
        rn = a["EA"] / a["ED"] if a["ED"] else float("nan")
        rec[nm] = {"rho_sub": rs, "rho_null": rn, "excess": rs - rn}
        print(f"  depth {d}  {nm:10s} rho_sub {rs:+.4f}  clone-null {rn:+.4f}  "
              f"EXCESS {rs - rn:+.4f}   ({n_anchor} anchors)", flush=True)
    out["depths"][str(d)] = rec

# ---------------------------------------------------------------- verify the null
# The entire result is observed-minus-analytic-null, so check the analytic null
# against explicit within-clone permutations on a handful of anchors.
if VERIFY:
    d = DEPTHS[0]
    r, r2, sv, v = resid["missing"]
    tested, rows = 0, []
    for zi in range(K):
        if tested >= 8:
            break
        cd = codes[:, zi, d - 1]
        valid = cd >= 0
        if valid.sum() < 2 * MIN_CLONE:
            continue
        cl = g_clone[valid]
        key = cl.astype(np.int64) * (cd.max() + 2) + cd[valid]
        _, sub = np.unique(key, return_inverse=True)
        first = np.zeros(sub.max() + 1, dtype=np.int64); first[sub] = cl
        ngrp = np.bincount(first, minlength=np.bincount(cl).size)
        good = (np.bincount(cl) >= MIN_CLONE) & (ngrp >= 2)
        m = good[cl]
        if m.sum() < 2 * MIN_CLONE:
            continue
        idx = np.flatnonzero(valid)[m]
        _, cl2i = np.unique(g_clone[idx], return_inverse=True)
        _, sub2 = np.unique(key[m], return_inverse=True)
        Gs, Gk = sub2.max() + 1, cl2i.max() + 1
        ms = np.bincount(sub2).astype(float)
        owner = np.zeros(Gs, dtype=np.int64); owner[sub2] = cl2i
        ns = np.bincount(cl2i).astype(float)
        wC = np.bincount(owner, weights=ms * (ms - 1), minlength=Gk) / np.maximum(ns * (ns - 1), 1e-9)
        EA = pair_sums(r2[idx], r[idx], cl2i, Gk, w=wC); EA[zi] = 0.0
        order = np.argsort(cl2i, kind="stable")
        bounds = np.concatenate([[0], np.cumsum(np.bincount(cl2i))])
        sims = []
        for _ in range(VERIFY):
            perm = order.copy()
            for a_, b_ in zip(bounds[:-1], bounds[1:]):
                perm[a_:b_] = rng.permutation(perm[a_:b_])
            inv = np.empty_like(perm); inv[order] = perm     # shuffle within clone
            A = pair_sums(r2[idx][inv], r[idx][inv], sub2, Gs); A[zi] = 0.0
            sims.append(A.sum())
        rows.append((zi, float(np.mean(sims)), float(np.std(sims, ddof=1)), float(EA.sum())))
        tested += 1
    print(f"\n  --- analytic null check, depth {d}, {VERIFY} permutations per anchor")
    for zi, mu, sd, ea in rows:
        print(f"      anchor {zi:>4}  permuted mean {mu:+.5g} +/- {sd:.3g}   "
              f"analytic {ea:+.5g}   ratio {mu/ea if ea else float('nan'):.4f}")
    rat = [mu / ea for _, mu, _, ea in rows if ea]
    print(f"      mean ratio {np.mean(rat):.4f} (1.0 = analytic null exact)")
    out["null_check"] = {"depth": d, "n_perm": VERIFY,
                         "mean_ratio": float(np.mean(rat)),
                         "rows": [{"anchor": int(a), "perm_mean": b, "perm_sd": c,
                                   "analytic": e} for a, b, c, e in rows]}

tag = "_screen" if SCREEN else ""
(RES / f"lineage_depth_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/lineage_depth_{arm}{tag}.json")
