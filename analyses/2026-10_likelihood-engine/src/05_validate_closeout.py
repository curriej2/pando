#!/usr/bin/env python3
r"""
05_validate_closeout.py -- closing out the engine validation: the tree-prior checks V8a and V8+.

Approved by Justin 2026-10-12 ("go ahead with the engine validation"; proposal in README, Session 25).
03_validate_engine.py is left untouched as the record of the 2026-10-06 run. V5' (gradients against
an improved finite-difference reference) is NOT here: a probe showed the proposed fixed step
(1e-4) does not resolve the failing entry, so V5' returns to Justin as a revised proposal.

V8a  (hard) the ENGINE's prior density on one branching height, integrated numerically from 0 to x,
     reproduces the simulator's CDF bds_G(x) -- the function V8 tests against the library. Closes
     the link "V8 tests bds_G, the engine uses log_tree_prior".                       tol 1e-8 abs
V8+  (statistical) fresh forward-simulated trees (2,000 per cell; seeds disjoint from the library):
     the prior's probability-integral transform u = G(h) of every branching height is uniform --
     decile occupancy (minus 0.1) averaged over trees, and the root's PIT G(h_max)^(n-1) - 0.5,
     each as mean / se across trees. The failing 2026-10-06 cell plus four neighbours.
     Reading: |z| <= 3 PASS; 3-5 marginal (report, re-draw); > 5 defect. At 2,000 trees the
     20-tree deficit (z ~ -3), if real, would return at z ~ -30.
     A regression first: the same statistic on the LIBRARY trees of the failing cell must reproduce
     the recorded z (-3.00, -3.03), so the fresh-tree numbers are the identical statistic.

usage: 05_validate_closeout.py [--trees 2000] [--procs 5]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import pathlib
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.integrate import quad

if "CUDA_VISIBLE_DEVICES" not in os.environ:
    os.environ.setdefault("JAX_PLATFORMS", "cpu")
HERE = pathlib.Path(__file__).resolve().parent
SIM = HERE.parents[1] / "2026-09_simulator" / "src"
RES = HERE.parent / "results"
FRESH = RES / "v8_fresh_trees"                      # .npz, gitignored


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


bd = load("bdtree", SIM / "01_bdtree.py")
ed = load("editing", SIM / "14_editing.py")

FAILING = (64, 0.25, 0.3)
CELLS = [FAILING, (64, 0.1, 0.3), (64, 0.25, 0.7), (25, 0.25, 0.3), (210, 0.25, 0.3)]
V8_GRID = [(n, rho, th) for n in (25, 64, 210) for rho in (0.002, 0.1, 0.25) for th in (0.3, 0.7)]
SEED0 = 61012026          # library uses 20260921 + 1,000,003 i; the offset is not a multiple of 1,000,003


# ---------------------------------------------------------------------------------------------------
# the V8 statistic, exactly as in 03_validate_engine.check_V8
# ---------------------------------------------------------------------------------------------------

def v8_stat(trees, n, rho, th):
    b, d, _ = bd.rate_params(n, rho, th)
    dec, rootp = [], []
    for parent, branch in trees:
        t = ed.node_times(parent, branch, n)
        h = 1.0 - t[n:]
        u = bd.bds_G(h, b, d, rho, 1.0)
        dec.append(np.histogram(u, bins=10, range=(0, 1))[0] / u.size - 0.1)
        rootp.append(float(bd.bds_G(h.max(), b, d, rho, 1.0)) ** (n - 1) - 0.5)
    dec, rootp = np.array(dec), np.array(rootp)
    z = lambda v: float(v.mean() / (v.std(ddof=1) / math.sqrt(v.size)))
    out = {f"decile{q}": z(dec[:, q]) for q in range(10)}
    out["root_pit"] = z(rootp)
    out["_mean_dev"] = [float(x) for x in dec.mean(axis=0)]
    out["_trees"] = int(dec.shape[0])
    return out


# ---------------------------------------------------------------------------------------------------
# V8+: fresh trees
# ---------------------------------------------------------------------------------------------------

def fresh_cell(args):
    (n, rho, th), R = args
    b, d, T = bd.rate_params(n, rho, th)
    f = FRESH / f"trees_n{n}_rho{rho:g}_th{th:g}_R{R}.npz"
    if f.exists():
        z = np.load(f)
        return (n, rho, th), list(zip(z["parent"].astype(np.int64), z["branch"])), "cached"
    t0 = time.perf_counter()
    trees, seeds = [], []
    for i in range(R):
        s = SEED0 + 1_000_003 * i
        tr, _ = bd.simulate_tree(b, d, T, rho, n, seed=s)
        assert tr is not None, f"no acceptance at seed {s}"
        bl = np.zeros(tr.parent.size)
        nz = tr.parent >= 0
        bl[nz] = tr.time[nz] - tr.time[tr.parent[nz]]
        trees.append((tr.parent.astype(np.int64), bl))
        seeds.append(s)
    np.savez_compressed(f, parent=np.array([p for p, _ in trees]), branch=np.array([q for _, q in trees]),
                        seed=np.array(seeds))
    return (n, rho, th), trees, f"{time.perf_counter() - t0:.0f} s"


# ---------------------------------------------------------------------------------------------------
# V8a: the engine's own density integrates to bds_G
# ---------------------------------------------------------------------------------------------------

def check_V8a():
    eng = load("engine", HERE / "02_engine.py")
    import jax
    import jax.numpy as jnp

    class _D:                                       # log_tree_prior reads only tip_obs.shape[0]
        tip_obs = np.zeros((1, 1))

    @jax.jit
    def logdens(h, b, d, rho):
        return eng.log_tree_prior(jnp.stack([jnp.ones(()), 1.0 - h]), b, d, rho, _D)

    worst, rows = 0.0, []
    xs = np.linspace(0.05, 1.0, 20)
    for n, rho, th in V8_GRID:
        b, d, _ = bd.rate_params(n, rho, th)
        dens = lambda h: float(np.exp(logdens(h, b, d, rho)))
        acc, prev = 0.0, 0.0
        for x in xs:                                # cumulative, interval by interval
            acc += quad(dens, prev, x, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
            prev = x
            e = abs(acc - float(bd.bds_G(x, b, d, rho, 1.0)))
            worst = max(worst, e)
            rows.append({"n": n, "rho": rho, "theta": th, "x": float(x), "engine_cdf": acc, "err": e})
    return {"configs": len(rows), "worst_abs": worst, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", type=int, default=2000)
    ap.add_argument("--procs", type=int, default=5)
    a = ap.parse_args()
    FRESH.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    # regression: the statistic on the library trees of the failing cell
    reg = v8_stat(ed.load_library_trees(*FAILING), *FAILING)
    rec = json.loads((RES / "validate_engine.json").read_text())["detail"]["V8"]
    key = "n64_rho0.25_th0.3:"
    reg_dev = max(abs(reg[f"decile{q}"] - rec[key + f"decile{q}"]) for q in range(10))
    print(f"regression on library trees: max |z - recorded z| = {reg_dev:.2e} "
          f"(decile0 {reg['decile0']:.2f}, decile1 {reg['decile1']:.2f})", flush=True)
    assert reg_dev < 1e-9, "statistic does not reproduce the recorded V8 values"

    v8a = check_V8a()
    print(f"V8a done: worst |engine CDF - bds_G| = {v8a['worst_abs']:.2e} over {v8a['configs']}", flush=True)

    with Pool(a.procs) as pool:
        got = pool.map(fresh_cell, [(c, a.trees) for c in CELLS])
    v8p = {}
    for cell, trees, how in got:
        s = v8_stat(trees, *cell)
        v8p[f"n{cell[0]}_rho{cell[1]:g}_th{cell[2]:g}"] = s
        zs = [s[f"decile{q}"] for q in range(10)] + [s["root_pit"]]
        print(f"V8+ {cell}: {s['_trees']} trees ({how}), worst |z| {max(abs(z) for z in zs):.2f}; "
              f"decile0 {s['decile0']:+.2f}, decile1 {s['decile1']:+.2f}", flush=True)

    zf = [abs(v[k]) for v in v8p.values() for k in v if not k.startswith("_")]
    zfail = [abs(v8p["n64_rho0.25_th0.3"][k]) for k in ("decile0", "decile1")]
    verdict = ("DEFECT" if max(zfail) > 5 or max(zf) > 5 else
               "MARGINAL" if max(zf) > 3 else "PASS")
    summary = {"V8a": {"configs": v8a["configs"], "worst_abs": v8a["worst_abs"], "tol": 1e-8,
                       "pass": v8a["worst_abs"] <= 1e-8},
               "V8plus": {"cells": len(v8p), "stat_cells": len(zf), "worst_abs_z": max(zf),
                          "failing_cell_decile01_abs_z": zfail, "verdict": verdict},
               "regression_max_dev": reg_dev, "trees_per_cell": a.trees, "seed0": SEED0,
               "seconds": round(time.perf_counter() - t0, 1)}
    (RES / "validate_closeout.json").write_text(json.dumps(
        {"summary": summary, "detail": {"V8a": v8a, "V8plus": v8p, "regression": reg}}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
