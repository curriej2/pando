#!/usr/bin/env python3
r"""
21_f4_followup.py -- the targeted follow-up to the failed rung-1 validation (approved by Justin
2026-09-30; README "Session 20").

WHY.  18_validate_dropout.py FAILED at seed 30092026: two trips that were one fluctuation counted
twice (F4o duplicates F4t where censoring is rare, z correlation 0.98-1.00), and an unexplained
excess in F4o on SLOW (10 of 10 W-bins positive, obs - exp +0.0004 ... +0.010).  This asks one
question: is that excess a defect in the generator / closed form, or chance?

PAIRS.  One tip from each side of every internal node x, per tape z; W_x = W(t_x) is the fraction
of all editing done by the split, mu_z = r_z Lam_T the tape's total expected edits.

CHECKS -- each POOLED into one pre-declared cell per fixture (all nodes and tapes of a unit):
  F4d  1[D1 = 0] 1[D2 = 0]  vs  e^{-mu_z (2 - W_x)}      censoring alone, no mask (NEW)
  F4t  X1 X2                vs  pi1 pi2                   the mask alone
  F4o  O1 O2                vs  the joint closed form of 18 (censoring x mask)
  The 10 W-bins are reported as secondary, outside the verdict.

SAMPLE.  All 840 stored n = 210 library trees (42 rho x theta cells x 20) x 5 panel draws, on the
SLOW and CIS fixtures of 18.  Unit = (panel, tree): given the panel, trees are independent and the
closed forms are exact conditional expectations, so residuals from different panels pool.
se across units.  New seed.
⚠ DIAGNOSTIC, outside the verdict (disclosed before the run): seed 30092026's exact SLOW panel
replayed on the same 840 trees, reported separately -- tests whether the excess belongs to that
panel's tapes.

READING RULE (fixed before the run): no-effect value z = 0.  |z| <= 3 on all six pooled cells
(3 checks x 2 fixtures) => no defect; the excess was chance.  Otherwise a defect: F4d failing
points at the editing layer or the formula, F4d passing with F4o failing at the mask x censoring
combination -- stop and trace the code.

usage:  21_f4_followup.py [--seed S] [--quick]
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("valdrop", _HERE / "18_validate_dropout.py")
vd = importlib.util.module_from_spec(_spec)
sys.modules["valdrop"] = vd
_spec.loader.exec_module(vd)
dl, ed = vd.dl, vd.ed

RES = _HERE.parent / "results"
N_PANELS, NBIN = 5, 10
REPLAY_SEED = 30092026            # the failed run; SLOW was fixture index 2 there


def tree_set(per_cell):
    cells = sorted({(float(p.name.split("_rho")[1].split("_")[0]),
                     float(p.name.split("_th")[1].split("_")[0]))
                    for p in ed.LIB.glob(f"trees_n{vd.N_TIPS}_rho*_th*_r*.npz")})
    trees = []
    for rho, th in cells:
        tt = ed.load_library_trees(vd.N_TIPS, rho, th)
        assert len(tt) >= per_cell, (rho, th, len(tt))
        trees += tt[:per_cell]
    return cells, trees


def one_unit(parent, branch, rec, beta, s_alpha, seeds):
    """Pooled and binned residuals of F4d / F4t / F4o for one (panel, tree)."""
    n, K = vd.N_TIPS, vd.K
    sdec, smask, spair = seeds
    D, _, out, _ = dl.decorate_tips(parent, branch, n, rec, np.random.default_rng(sdec))
    md = dl.mask_draws(n, K, np.random.default_rng(smask))
    X, pi = dl.technical_mask(md, beta, s_alpha)
    O = ~dl.observe(X, D)
    nodes, c1, c2 = vd.split_pairs(out["parent"], n, np.random.default_rng(spair))
    Wa = rec.W_of(out["t"][nodes])[:, None]
    mu = (rec.Lam_T * rec.r)[None, :]
    ee = np.exp(-mu)
    both0 = np.exp(-mu * (2 - Wa))
    p1, p2 = pi[c1], pi[c2]
    obs = {"F4d": ((D[c1] == 0) & (D[c2] == 0)).astype(float),
           "F4t": (X[c1] & X[c2]).astype(float),
           "F4o": (O[c1] & O[c2]).astype(float)}
    exp_ = {"F4d": np.broadcast_to(both0, obs["F4d"].shape),
            "F4t": p1 * p2,
            "F4o": 1 - (1 - p1) * (1 - ee) - (1 - p2) * (1 - ee)
                   + (1 - p1) * (1 - p2) * (1 - 2 * ee + both0)}
    wb = np.minimum((Wa[:, 0] * NBIN).astype(int), NBIN - 1)
    res = {}
    for c in obs:
        r = obs[c] - exp_[c]
        res[c] = (r.mean(), obs[c].mean(), exp_[c].mean(),
                  [(r[wb == b].mean(), r[wb == b].size) if (wb == b).any() else (np.nan, 0)
                   for b in range(NBIN)])
    return res


def summarise(units):
    """units: list of res dicts -> pooled z per check and binned z (se across units)."""
    out = {}
    for c in ("F4d", "F4t", "F4o"):
        r = np.array([u[c][0] for u in units])
        o = np.array([u[c][1] for u in units]); e = np.array([u[c][2] for u in units])
        z = float(r.mean() / (r.std(ddof=1) / np.sqrt(r.size)))
        bins = []
        for b in range(NBIN):
            rb = np.array([u[c][3][b][0] for u in units if u[c][3][b][1] > 0])
            bins.append({"bin": b, "units": int(rb.size),
                         "z": float(rb.mean() / (rb.std(ddof=1) / np.sqrt(rb.size)))
                         if rb.size > 20 else None,
                         "mean_resid": float(rb.mean()) if rb.size else None})
        out[c] = {"units": int(r.size), "obs": float(o.mean()), "exp": float(e.mean()),
                  "mean_resid": float(r.mean()), "se": float(r.std(ddof=1) / np.sqrt(r.size)),
                  "z": z, "bins": bins}
    return out


def run_panel(fx, pdraw, trees, seed, tag):
    u, r = dl.panel_tapes(pdraw, fx["panel"])
    beta, _ = dl.panel_beta(pdraw, u, fx["panel"], vd.MASK)
    rec = ed.Recorder(k=vd.K, Lam_T=fx["Lam_T"], r=r, shares=((0.0,),), **fx["rate"])
    units = []
    for ti, (parent, branch) in enumerate(trees):
        seeds = [np.random.SeedSequence([seed, s, *tag, ti]) for s in (2, 3, 4)]
        units.append(one_unit(parent, branch, rec, beta, vd.MASK.s_alpha, seeds))
    mu = rec.Lam_T * rec.r
    return units, {"closed_tapes": int(u.sum()), "mu_min": float(mu.min()),
                   "mu_max": float(mu.max()), "P_unedited_max": float(np.exp(-mu.min()))}


def main():
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 1102026
    quick = "--quick" in sys.argv
    t0 = time.perf_counter()
    cells, trees = tree_set(2 if quick else 20)
    npan = 2 if quick else N_PANELS
    result = {"settings": {"seed": seed, "n": vd.N_TIPS, "tree_cells": cells,
                           "trees": len(trees), "panels_per_fixture": npan,
                           "replay_seed": REPLAY_SEED}, "fixtures": {}, "replay": {}}
    for name in ("SLOW", "CIS"):
        fx = vd.FIXTURES[name]
        ci = list(vd.FIXTURES).index(name)
        allu, pinfo, per_panel = [], [], []
        for p in range(npan):
            pdraw = dl.panel_draws(vd.K, np.random.default_rng(
                np.random.SeedSequence([seed, 11, ci, p])))
            units, info = run_panel(fx, pdraw, trees, seed, (ci, p))
            s = summarise(units)
            info.update({c: {"z": s[c]["z"], "mean_resid": s[c]["mean_resid"]} for c in s})
            pinfo.append(info)
            allu += units
            print(f"{name} panel {p}: {info}  ({time.perf_counter() - t0:.0f} s)", flush=True)
        result["fixtures"][name] = {"pooled": summarise(allu), "per_panel": pinfo}
    # diagnostic: the failed run's SLOW panel, replayed on the same trees
    fx = vd.FIXTURES["SLOW"]
    pdraw = dl.panel_draws(vd.K, np.random.default_rng(np.random.SeedSequence([REPLAY_SEED, 1, 2])))
    units, info = run_panel(fx, pdraw, trees, seed, (99, 0))
    result["replay"] = {"panel_info": info, "pooled": summarise(units)}
    zs = [abs(result["fixtures"][f]["pooled"][c]["z"]) for f in ("SLOW", "CIS")
          for c in ("F4d", "F4t", "F4o")]
    result["summary"] = {"worst_pooled_abs_z": float(max(zs)),
                         "verdict": "NO DEFECT" if max(zs) <= 3.0 else "DEFECT"}
    result["seconds"] = round(time.perf_counter() - t0, 1)
    (RES / ("f4_followup_quick.json" if quick else "f4_followup.json")).write_text(
        json.dumps(result, indent=1, default=float))

    print(f"\n{'fixture':8s} {'check':4s} {'units':>6s} {'obs':>9s} {'exp':>9s} {'obs-exp':>10s} {'z':>6s}")
    for f in ("SLOW", "CIS"):
        for c in ("F4d", "F4t", "F4o"):
            s = result["fixtures"][f]["pooled"][c]
            print(f"{f:8s} {c:4s} {s['units']:6d} {s['obs']:9.5f} {s['exp']:9.5f} "
                  f"{s['mean_resid']:+10.6f} {s['z']:+6.2f}")
    for c in ("F4d", "F4t", "F4o"):
        s = result["replay"]["pooled"][c]
        print(f"{'REPLAY':8s} {c:4s} {s['units']:6d} {s['obs']:9.5f} {s['exp']:9.5f} "
              f"{s['mean_resid']:+10.6f} {s['z']:+6.2f}   (diagnostic)")
    print(f"\nworst pooled |z| {result['summary']['worst_pooled_abs_z']:.2f}  =>  "
          f"{result['summary']['verdict']}   ({result['seconds']} s)")


if __name__ == "__main__":
    main()
