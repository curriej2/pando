#!/usr/bin/env python3
r"""
FEASIBILITY / POWER for the A9 lineage question -- NOT the test itself.

Question to be answered later: do related cells lose the SAME tapes?

What limits us: a cell enters only if it (i) clears the arm's tape filter and
(ii) carries a ClonalBC -- and Fig 3 showed barcode loss tracks R_c, so the
clones we can work in are built from high-capture cells.  This script counts what
survives, and computes the null spread of the statistic the test will use, so the
figure is planned against real power rather than hope.

STATISTIC THE TEST WILL USE.  Within clone C, for tape z, let
    a_{C,z} = fraction of C's cells in which tape z is entirely unrecovered.
Under dropout that is independent given the marginals, a_{C,z} is an average of
independent Bernoullis with cell-specific rates (1 - alpha_c beta_z), so

    Var_null(a_{C,z}) = (1/n_C^2) sum_{c in C} pi_c (1 - pi_c),   pi_c = P(miss)

which for a clone of n_C cells is roughly p_miss(1-p_miss)/n_C.  Heritable loss
puts excess mass at a = 0 and a = 1.  We report that null sd per clone-size band
so we know the smallest detectable excess.

Also counted: INFORMATIVE absence characters -- (clone, tape) pairs where the
absence set is neither empty nor everything (3 <= |A| <= n_C - 3), since only
those can be tested for compatibility against the edit characters.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ARMS = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
BANDS = [(10, 20), (20, 50), (50, 100), (100, 500), (500, 10**9)]

out = {}
for a in ARMS:
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y, cl, sm = z["recovered"], z["clone"], z["sample"]
    keep = (Y.sum(1) >= THR[a]) & (cl >= 0)
    Y, cl, sm = Y[keep], cl[keep], sm[keep]
    n, k = Y.shape
    ids, inv = np.unique(cl, return_inverse=True)
    sizes = np.bincount(inv)
    d = {"cells_usable": int(n), "clones": int(ids.size),
         "median_clone": float(np.median(sizes)), "max_clone": int(sizes.max())}
    print(f"\n=== {a}: {n:,} usable cells in {ids.size:,} clones "
          f"(median {np.median(sizes):.0f}, max {sizes.max():,})")
    for lo, hi in BANDS:
        sel = np.where((sizes >= lo) & (sizes < hi))[0]
        if sel.size == 0:
            continue
        cells = int(sizes[sel].sum())
        # informative absence characters + null sd of a_{C,z}, over these clones
        ninf, nulls = 0, []
        for ci in sel:
            m = inv == ci
            nc = int(m.sum())
            miss = ~Y[m]
            cnt = miss.sum(0)
            ninf += int(((cnt >= 3) & (cnt <= nc - 3)).sum())
            pi = miss.mean(1)                       # per-cell miss propensity
            nulls.append(np.sqrt((pi * (1 - pi)).sum()) / nc)
        d[f"{lo}-{hi if hi < 10**8 else 'inf'}"] = {
            "clones": int(sel.size), "cells": cells,
            "informative_absence_chars": ninf,
            "null_sd_a": float(np.mean(nulls))}
        print(f"  size {lo:>4}-{'inf' if hi>10**8 else hi:<4} : {sel.size:>5} clones, "
              f"{cells:>7,} cells, {ninf:>8,} informative absence characters, "
              f"null sd of a_C,z = {np.mean(nulls):.4f}")
    # how much of each clone sits in one sample -- the batch confound for the test
    top = []
    for ci in np.where(sizes >= 10)[0]:
        s = sm[inv == ci]
        top.append(np.bincount(s).max() / s.size)
    if top:
        d["clone_sample_purity_median"] = float(np.median(top))
        print(f"  clones >=10 cells: median share in a single sample = "
              f"{np.median(top):.2f}  (1.00 => batch fully confounded with clone)")
    out[a] = d

(RES / "lineage_feasibility.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/lineage_feasibility.json")
