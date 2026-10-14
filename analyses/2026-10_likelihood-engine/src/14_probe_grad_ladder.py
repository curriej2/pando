#!/usr/bin/env python3
r"""
14_probe_grad_ladder.py -- step ladder for the one V5 outlier (validation figure 4, panel b).

Approved by Justin 2026-10-13 ("you can make this figure"; option ii: re-run the ladder rather than plot
numbers typed from the README). Reproduces EXACTLY the V5 entry that failed on 2026-10-06 -- tree 16, the
"const" setting (Gam_T = 0, no silencing), d(J)/d(Gam_T) -- by replaying 03_validate_engine.check_V5's set-up
at seed 6102026, asserts the recorded autodiff and h = 1e-5 values, then evaluates central finite
differences at h = 1e-3 ... 1e-8. Nothing in the engine or the check is changed.
Output: results/v5_ladder.json.  Env: envs/pando-jax.  Seconds.
"""
from __future__ import annotations

import importlib.util
import json
import os
import pathlib
import sys

import numpy as np

if "CUDA_VISIBLE_DEVICES" not in os.environ:
    os.environ.setdefault("JAX_PLATFORMS", "cpu")
HERE = pathlib.Path(__file__).resolve().parent
RES = HERE.parent / "results"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v = load("validate_engine", HERE / "03_validate_engine.py")
jax, jnp, ed, bd, eng, ref = v.jax, v.jnp, v.ed, v.bd, v.eng, v.ref

SEED, TREE = 6102026, 16
REC_AD, REC_FD = 2679.07139773103, 2679.7352913490613          # validate_engine.json, V5 row (tree 16, Gam_T)

# --- replay check_V5's set-up for tree i = 16 (same cells list, setting rotation, seeds) ---
cells = [(25, 0.1, 0.3), (64, 0.25, 0.7), (25, 0.002, 0.7), (64, 0.1, 0.3)]
names = list(v.SETTINGS)
i = TREE
n, rho, th = cells[i % len(cells)]
parent, branch = ed.load_library_trees(n, rho, th)[i // len(cells)]
st = names[i % len(names)]
sim = v.simulate(parent, branch, n, v.SETTINGS[st], [SEED, 5, i])
prm, structure, data, _ = v.engine_inputs(sim, parent, n)
b, d, _ = bd.rate_params(n, rho, th)
J = eng.make_objective(v.loglik_fn(structure))
x0 = jnp.asarray(ref.x_from_times(sim["t"], parent, n))
theta0 = dict(prm, b=jnp.asarray(b), d=jnp.asarray(d))


def f(x, th_):
    pr = {kk: vv for kk, vv in th_.items() if kk not in ("b", "d")}
    return J(x, pr, data, (th_["b"], th_["d"], rho))


def at_gam(g):
    th2 = dict(theta0)
    th2["Gam_T"] = theta0["Gam_T"] + g
    return float(f(x0, th2))


ad = float(jax.grad(f, argnums=1)(x0, theta0)["Gam_T"])
fd5 = (at_gam(1e-5) - at_gam(-1e-5)) / 2e-5
print(f"tree {i}: n={n} rho={rho} theta={th} setting={st} Gam_T={float(theta0['Gam_T'])}")
print(f"autodiff {ad!r} (recorded {REC_AD!r});  central fd h=1e-5 {fd5!r} (recorded {REC_FD!r})")
assert abs(ad - REC_AD) <= 1e-9 * abs(REC_AD), "autodiff does not reproduce the recorded V5 value"
assert abs(fd5 - REC_FD) <= 1e-7 * abs(REC_FD), "h = 1e-5 difference does not reproduce the recorded V5 value"

ladder = []
for h in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
    fd = (at_gam(h) - at_gam(-h)) / (2 * h)
    row = {"h": h, "fd": fd, "rel_err_vs_autodiff": abs(fd - ad) / max(abs(fd), 1.0) if np.isfinite(fd) else None}
    ladder.append(row)
    print(f"  h={h:.0e}  fd={fd!r}  rel err vs autodiff={row['rel_err_vs_autodiff']}")

out = {"tree": i, "n": n, "rho": rho, "theta": th, "setting": st, "Gam_T": float(theta0["Gam_T"]),
       "J_at_truth": at_gam(0.0), "autodiff": ad, "recorded_autodiff": REC_AD, "recorded_fd_h1e-5": REC_FD,
       "tol_V5": 1e-5, "err_definition": "|fd - autodiff| / max(|fd|, 1), as in V5", "ladder": ladder}
(RES / "v5_ladder.json").write_text(json.dumps(out, indent=1))
print(f"wrote {RES / 'v5_ladder.json'}")
