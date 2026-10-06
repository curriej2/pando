#!/usr/bin/env python3
r"""
04_time_engine.py -- timing of the JAX likelihood engine vs clone size (DESCRIPTIVE ONLY).

Approved 2026-10-06 (PROPOSAL PART 2, ../CLAUDE.md). For one library tree per size, simulated rung-2
data (silencing on), it reports seconds per evaluation of all k tapes' log-likelihoods and per
gradient (w.r.t. node times and every parameter), on whatever device JAX sees (CPU or GPU), double
precision. Compile time is reported separately. Nothing here enters any verdict.

usage: 04_time_engine.py [--k 30] [--sizes 25,210,1976,11081]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("validate_engine", HERE / "03_validate_engine.py")
val = importlib.util.module_from_spec(spec)
sys.modules["validate_engine"] = val
spec.loader.exec_module(val)
import jax                                  # noqa: E402
import jax.numpy as jnp                     # noqa: E402

CELL = {25: (0.1, 0.3), 210: (0.1, 0.3), 1976: (0.1, 0.3), 11081: (0.05, 0.3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--sizes", default="25,210,1976,11081")
    ap.add_argument("--reps", type=int, default=5)
    a = ap.parse_args()
    val.K_TAPES = a.k
    dev = jax.devices()[0]
    rows = []
    for n in [int(s) for s in a.sizes.split(",")]:
        rho, th = CELL[n]
        parent, branch = val.ed.load_library_trees(n, rho, th)[0]
        tsim = time.perf_counter()
        sim = val.simulate(parent, branch, n, val.SETTINGS["sil"], [4, n])
        tsim = time.perf_counter() - tsim
        prm, structure, data, _ = val.engine_inputs(sim, parent, n)
        ll = val.loglik_fn(structure)
        t = jnp.asarray(sim["t"])
        g = jax.jit(jax.grad(lambda tt, pp: jnp.sum(ll(tt, pp, data)), argnums=(0, 1)))
        c0 = time.perf_counter(); ll(t, prm, data).block_until_ready(); c_ll = time.perf_counter() - c0
        c0 = time.perf_counter(); jax.block_until_ready(g(t, prm)); c_g = time.perf_counter() - c0
        e_ll, e_g = [], []
        for _ in range(a.reps):
            s0 = time.perf_counter(); ll(t, prm, data).block_until_ready(); e_ll.append(time.perf_counter() - s0)
            s0 = time.perf_counter(); jax.block_until_ready(g(t, prm)); e_g.append(time.perf_counter() - s0)
        row = {"n": n, "k": a.k, "device": str(dev), "simulate_s": tsim, "compile_loglik_s": c_ll,
               "compile_grad_s": c_g, "loglik_s": float(np.median(e_ll)), "grad_s": float(np.median(e_g)),
               "loglik_total_nats": float(jnp.sum(ll(t, prm, data)))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    tag = "gpu" if dev.platform == "gpu" else "cpu"
    out = HERE.parent / "results" / f"timing_{tag}_k{a.k}.json"
    out.write_text(json.dumps(rows, indent=1))
    print(f"\n{'n':>6s} {'loglik s':>10s} {'grad s':>10s} {'compile s':>10s}   device {dev}")
    for r in rows:
        print(f"{r['n']:6d} {r['loglik_s']:10.4f} {r['grad_s']:10.4f} {r['compile_loglik_s'] + r['compile_grad_s']:10.1f}")


if __name__ == "__main__":
    main()
