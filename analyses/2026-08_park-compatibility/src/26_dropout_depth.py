#!/usr/bin/env python3
r"""
================================================================================
 Is dropout INFORMATIVE?  -- coupling between recovery and edit depth
================================================================================

Dropout that is independent of the editing state is a nuisance; dropout that is
correlated with how much a tape has recorded biases every downstream time-like
quantity (lambda, branch lengths, node ages), because the tapes we fail to see
are not a random sample of the tapes that exist.

TWO MARGINS, TWO CONFOUNDS.

(1) PER TAPE.  x = fraction of cells recovering tape z; y = mean depth given
    recovered.  A positive slope is the row-A9 prediction: pegRNA and TargetBC are
    co-integrated, so a silenced locus is both unreadable AND unedited -- one
    latent cause, two symptoms.
    CONFOUND: a low-recovery tape is seen mostly in HIGH-capture cells, which may
    themselves have deeper tapes.  Control: recompute y using only cells with
    R_c >= HI_CELL, where every tape is seen in a comparable cell population.

(2) PER CELL.  x = R_c; y = mean depth over that cell's recovered tapes.  A rising
    slope means the cells we see worst also look YOUNGEST, so naive marginalisation
    drags lambda and node times downward.
    CONFOUND: a low-R_c cell recovers preferentially the EASY tapes, and easy tapes
    may differ in depth.  Control: recompute y over a fixed reference set of tapes
    (the TOP_TAPES most-recovered ones), identical for every cell.

DEPTH is the prefix length from 10_build_characters.py: leading sites carrying an
alphabet symbol.  Junk-terminated tapes are reported separately -- a junk value
censors depth for a reason unrelated to editing, so it must not be pooled with
UNEDITED (a genuine "had not got that far").

Output: results/dropout_depth.json
================================================================================
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ARMS = ["Initial", "Subclone", "Mouse1", "Mouse2", "Mouse3"]
ABSENT, UNEDITED, JUNK, COMPLETE = 0, 1, 2, 3
HI_CELL, TOP_TAPES = 140, 30


def spearman(x, y):
    n = len(x)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    r = float(np.corrcoef(rx, ry)[0, 1])
    z = 0.5 * np.log((1 + r) / (1 - r)) if abs(r) < 1 else float("inf")
    se = 1 / np.sqrt(n - 3)
    lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    return {"rho": r, "n": int(n), "ci": [float(lo), float(hi)]}


out = {}
for a in ARMS:
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y, dep, trm = z["recovered"], z["depth"].astype(float), z["term"]
    n, K = Y.shape
    ok = Y & (trm != JUNK)              # recovered AND depth not censored by junk
    Rc = Y.sum(1)
    rate = Y.mean(0)                    # per-tape recovery
    d = {"n_cells": n, "junk_frac_of_recovered": float((trm == JUNK).sum() / Y.sum())}

    # ---------------- (1) per tape
    with np.errstate(invalid="ignore"):
        depth_t = np.where(ok.sum(0) > 0, dep.sum(0, where=ok) / np.maximum(ok.sum(0), 1), np.nan)
    hi = Rc >= HI_CELL
    okh = ok[hi]
    depth_t_hi = np.where(okh.sum(0) > 20,
                          dep[hi].sum(0, where=okh) / np.maximum(okh.sum(0), 1), np.nan)
    m = np.isfinite(depth_t)
    d["per_tape"] = {"raw": spearman(rate[m], depth_t[m]),
                     "rate_range": [float(rate.min()), float(rate.max())],
                     "depth_range": [float(np.nanmin(depth_t)), float(np.nanmax(depth_t))]}
    mh = np.isfinite(depth_t_hi)
    d["per_tape"]["controlled"] = spearman(rate[mh], depth_t_hi[mh])
    d["per_tape"]["n_hi_cells"] = int(hi.sum())
    d["per_tape"]["points"] = {"rate": rate.tolist(), "depth": np.where(m, depth_t, np.nan).tolist(),
                               "depth_hi": np.where(mh, depth_t_hi, np.nan).tolist()}

    # ---------------- (2) per cell
    nok = ok.sum(1)
    mean_depth = np.where(nok > 0, dep.sum(1, where=ok) / np.maximum(nok, 1), np.nan)
    ref = np.argsort(rate)[::-1][:TOP_TAPES]          # fixed easy-tape reference set
    okr = ok[:, ref]
    nref = okr.sum(1)
    depth_ref = np.where(nref >= 5,
                         dep[:, ref].sum(1, where=okr) / np.maximum(nref, 1), np.nan)
    s = np.isfinite(mean_depth)
    sr = np.isfinite(depth_ref)
    d["per_cell"] = {"raw": spearman(Rc[s].astype(float), mean_depth[s]),
                     "controlled": spearman(Rc[sr].astype(float), depth_ref[sr]),
                     "n_ref_tapes": TOP_TAPES}
    edges = [0, 20, 40, 60, 80, 100, 120, 140, 167]
    prof = []
    for lo, h_ in zip(edges[:-1], edges[1:]):
        sel = (Rc >= lo) & (Rc < h_)
        if sel.sum() < 30: continue
        prof.append({"lo": lo, "hi": h_, "n": int(sel.sum()),
                     "mean_depth": float(np.nanmean(mean_depth[sel])),
                     "mean_depth_ref": float(np.nanmean(depth_ref[sel])),
                     "frac_full": float(np.nanmean((dep[sel] == 6)[ok[sel]])) })
    d["per_cell"]["profile"] = prof
    out[a] = d

    pt, pc = d["per_tape"], d["per_cell"]
    print(f"\n{a}  ({n:,} cells, junk-censored {100*d['junk_frac_of_recovered']:.1f}% of recovered)")
    print(f"  per tape : rho raw {pt['raw']['rho']:+.3f} "
          f"[{pt['raw']['ci'][0]:+.3f},{pt['raw']['ci'][1]:+.3f}]   "
          f"controlled (R_c>={HI_CELL}, n={pt['n_hi_cells']:,} cells) "
          f"{pt['controlled']['rho']:+.3f} "
          f"[{pt['controlled']['ci'][0]:+.3f},{pt['controlled']['ci'][1]:+.3f}]")
    print(f"  per cell : rho raw {pc['raw']['rho']:+.3f}   "
          f"controlled (top-{TOP_TAPES} tapes) {pc['controlled']['rho']:+.3f} "
          f"[{pc['controlled']['ci'][0]:+.3f},{pc['controlled']['ci'][1]:+.3f}]")
    print("  R_c bin -> mean depth (all recovered / reference tapes / frac full):")
    for b in prof:
        print(f"    {b['lo']:3d}-{b['hi']:3d}  n={b['n']:6,d}  "
              f"{b['mean_depth']:.3f} / {b['mean_depth_ref']:.3f} / {100*b['frac_full']:.1f}%")

(RES / "dropout_depth.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/dropout_depth.json")
