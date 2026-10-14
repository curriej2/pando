#!/usr/bin/env python3
r"""
06_validate_sciphy.py -- V6: our engine against SciPhy's own likelihood, on the base model.

Approved by Justin 2026-10-12 (engine-validation closeout; proposal in README, Session 25). The only
check written by someone else: it catches a misconception the engine and the simulator might share,
for the part of the model SciPhy has.

BASE MODEL = what SciPhy computes: constant editing rate per tape (clock = r_z Lam_T per experiment),
one constant composition xi, no silencing (Gam_T = 0), no technical dropout (alpha_c = -800, so
pi = 0 exactly), founder blank at the origin t = 0, tips at t = 1.
SciPhy side: SciPhyTreeLikelihood from azwaans/SciPhy (commit 2a49194), source unmodified, compiled
against BEAST 2.7.7; one likelihood per tape (src/06_sciphy_harness/V6Harness.java).

Two codings of a blank (unedited) tape on our side, both must agree with SciPhy:
  V6a  observed blank (SciPhy's own coding: the tip IS the empty string);
  V6b  blank reported MISSING (Park's coding); with pi = 0 the missing-tip vector is 1[d = 0], which is
       the same probability, but those cells now enter through the floating (all-missing-subtree) path.
Tolerance: relative 1e-8 per tape.

usage: 06_validate_sciphy.py [--configs 50] [--seed 12102026]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pathlib
import subprocess
import sys
import time

import numpy as np

if "CUDA_VISIBLE_DEVICES" not in os.environ:
    os.environ.setdefault("JAX_PLATFORMS", "cpu")
HERE = pathlib.Path(__file__).resolve().parent
RES = HERE.parent / "results"
TOOLS = pathlib.Path("/data1/choij10/justin/tools")
JAVA = TOOLS / "beast_dl/beast/jre/bin/java"
CP = ":".join([str(TOOLS / "beast_dl/beast/lib/packages/BEAST.base.jar"), str(TOOLS / "beast_dl/beast/lib/launcher.jar"),
               str(TOOLS / "sciphy_classes"), str(TOOLS / "v6_harness_classes")])   # harness compiled from src/06_sciphy_harness
WORK = pathlib.Path("/data1/choij10/justin/claude_scratch/v6_inputs")      # regenerable, not committed


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


v = load("validate_engine", HERE / "03_validate_engine.py")    # same simulator objects and engine inputs
jnp, dl, ed, sl = v.jnp, v.dl, v.ed, v.sl

CELLS = [(n, rho, th) for n in (25, 64, 210) for rho in (0.002, 0.1, 0.25) for th in (0.3, 0.7)]
LAMS = (1.0, 2.5, 4.0, 5.5, 8.0)       # edits/tape over the experiment: sparse ... past saturation (N = 6)
BASE = dict(rate=v.CONST, comp=None, Hbar=0.0, sil=((0, 1), (1,)))


def newick(parent, t, n):
    nn = len(parent)
    kids = [[] for _ in range(nn)]
    for c in range(nn - 1):
        kids[parent[c]].append(c)

    def rec(x):
        s = f"t{x}" if x < n else "(" + ",".join(rec(c) for c in kids[x]) + ")"
        return s if x == nn - 1 else f"{s}:{repr(float(t[x] - t[parent[x]]))}"
    return rec(nn - 1) + ";"


def one_config(i, seed):
    n, rho, th = CELLS[i % len(CELLS)]
    lam = LAMS[i % len(LAMS)]
    parent, branch = ed.load_library_trees(n, rho, th)[i // len(CELLS)]
    rng = np.random.default_rng([seed, i])
    pdraw = dl.panel_draws(v.K_TAPES, rng)
    _, r = dl.panel_tapes(pdraw, v.PANEL)
    rec = v.make_recorder(BASE, r, Lam_T=lam)
    D, sym, out, _ = dl.decorate_tips(parent, branch, n, rec, rng)
    t = out["t"]
    k, N = D.shape[1], v.N6
    sp = sl.SilencingParams(0.0, (0, 1), (1,))
    res = {}
    for coding, Y in (("V6a_observed_blank", np.ones_like(D, bool)), ("V6b_blank_as_missing", D > 0)):
        sim = dict(t=t, D=D, sym=sym, Y=Y, alpha=np.full(n, -800.0), beta=np.zeros(k), r=r, rec=rec, sp=sp)
        prm, structure, data, _ = v.engine_inputs(sim, parent, n)
        res[coding] = np.asarray(v.loglik_fn(structure)(jnp.asarray(t), prm, data))
    ck, xi = v.eng.comp_epochs_from_recorder(rec)
    assert xi.shape[0] == 1
    xi = xi[0] / xi[0].sum()
    M = xi.size
    lines = [f"{n} {k} {N} {M}", newick(parent, t, n), " ".join(repr(float(x)) for x in xi)]
    for z in range(k):
        lines.append(repr(float(rec.r[z] * rec.Lam_T)))
        for c in range(n):
            s = [int(sym[c, z, j]) + 1 if j < D[c, z] else 0 for j in range(N)]
            lines.append(",".join(map(str, s)))
    f = WORK / f"config_{i:03d}.txt"
    f.write_text("\n".join(lines) + "\n")
    meta = dict(config=i, n=n, rho=rho, theta=th, Lam_T=lam, M=M,
                full_tape_share=float((D == N).mean()), blank_share=float((D == 0).mean()))
    return f, res, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", type=int, default=50)
    ap.add_argument("--seed", type=int, default=12102026)
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    built = [one_config(i, a.seed) for i in range(a.configs)]
    print(f"engine side done ({time.perf_counter() - t0:.0f} s); running SciPhy", flush=True)
    # teeth: config 0 again with every SciPhy clock rate x (1 + 1e-4) -- must disagree beyond tol
    f0, _, m0 = built[0]
    L = f0.read_text().splitlines()
    for z in range(v.K_TAPES):
        j = 3 + z * (m0["n"] + 1)
        L[j] = repr(float(L[j]) * (1 + 1e-4))
    fneg = WORK / "teeth_rate_1e-4.txt"
    fneg.write_text("\n".join(L) + "\n")
    subprocess.run([str(JAVA), "-Dbeast.base.versionxml=" + str(TOOLS / "beast_dl/beast/lib/packages/BEAST.base.version.xml"),
                    "-cp", CP, "V6Harness"] + [str(f) for f, _, _ in built] + [str(fneg)], check=True)
    neg = np.loadtxt(str(fneg) + ".sciphy_ll")
    teeth = float(np.max(np.abs(built[0][1]["V6a_observed_blank"] - neg) / np.abs(neg)))
    print(f"teeth: rates x (1 + 1e-4) -> worst rel {teeth:.2e} (must exceed 1e-8)", flush=True)

    rows, worst = [], {"V6a_observed_blank": 0.0, "V6b_blank_as_missing": 0.0}
    for f, res, meta in built:
        ls = np.loadtxt(str(f) + ".sciphy_ll")
        for coding, le in res.items():
            rel = np.abs(le - ls) / np.abs(ls)
            worst[coding] = max(worst[coding], float(rel.max()))
            rows.append(dict(meta, coding=coding, worst_rel=float(rel.max()),
                             sum_ll_engine=float(le.sum()), sum_ll_sciphy=float(ls.sum())))
    tol = 1e-8
    summary = {c: {"configs": a.configs, "tape_comparisons": a.configs * v.K_TAPES, "worst_rel": w,
                   "tol": tol, "pass": w <= tol} for c, w in worst.items()}
    summary.update(seed=a.seed, sciphy_commit="2a49194", beast="2.7.7", teeth_rate_1e4_worst_rel=teeth,
                   seconds=round(time.perf_counter() - t0, 1))
    (RES / "validate_sciphy.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=1))
    for c, s in summary.items():
        if isinstance(s, dict):
            print(f"{c:24s} {s['tape_comparisons']:6d} tape comparisons  worst rel {s['worst_rel']:.2e}  "
                  f"tol {tol:.0e}  {'PASS' if s['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
