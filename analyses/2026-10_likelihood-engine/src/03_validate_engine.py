#!/usr/bin/env python3
r"""
03_validate_engine.py -- validation of the likelihood engine (01 reference twin, 02 JAX engine).

Approved by Justin 2026-10-06 ("PROPOSAL (2026-10-06), PART 2", ../CLAUDE.md). A failed check is a
BUG, not a finding. Runs in /data1/choij10/justin/envs/pando-jax.

HARD CHECKS (exact up to floating point; worst relative error vs a stated tolerance)
  V1  twin vs (a) dense pruning on the FULL state space (no slot/lcp logic) and (b) literal
      enumeration of every hidden history; tiny trees, N = 3, M = 3 (41 states incl. S); reads
      SAMPLED from the full-space model, plus constructed patterns (a whole subtree missing, nothing
      read); time-varying rate, two composition epochs, silencing, dropout, censoring.   tol 1e-10
  V2  sum of L over EVERY possible set of reads (each cell: any string of length 1..N, or missing)
      equals 1 -- the likelihood is a proper distribution over what can be observed.      tol 1e-10
  V3  JAX engine vs twin, per-tape ell_z, on simulator data on library trees.              tol 1e-10
  V4  closed-form segment kernel vs matrix exponential of the slot-chain generator.        tol 1e-12
  V5  JAX gradients of the objective vs central finite differences (node coordinates and
      parameters), relative to max(|g|, 1).                                               tol 1e-5
STATISTICAL CHECKS (|z| vs 0, se across trees; PASS iff worst |z| <= 3; one trip -> re-run)
  V7  score at the truth: mean over simulated datasets of d(sum ell)/d(parameter) at the true values.
  V8  tree prior: library branching heights, transformed by the prior's CDF, are uniform -- decile
      occupancy and the root's probability-integral transform, per library cell.
(V6, SciPhy's Java likelihood, is a separate setup step.)

usage: 03_validate_engine.py [--seed S] [--quick]
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import os
import pathlib
import sys
import time

import numpy as np
from scipy.linalg import expm

if "CUDA_VISIBLE_DEVICES" not in os.environ:          # no GPU allocated: use the CPU, skip the CUDA probe
    os.environ.setdefault("JAX_PLATFORMS", "cpu")
HERE = pathlib.Path(__file__).resolve().parent
SIM = HERE.parents[1] / "2026-09_simulator" / "src"
RES = HERE.parent / "results"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ref = load("engine_ref", HERE / "01_engine_ref.py")
eng = load("engine", HERE / "02_engine.py")
dl = load("dropout", SIM / "17_dropout.py")
ed = dl.ed
sl = load("silencing", SIM / "22_silencing.py")
bd = load("bdtree", SIM / "01_bdtree.py")
import jax                                      # noqa: E402  (after eng sets x64)
import jax.numpy as jnp                         # noqa: E402

K_TAPES, N6 = 30, 6


# ===================================================================================================
# helpers: tiny random trees, sampling reads from the FULL-space model
# ===================================================================================================

def random_tree(n, rng):
    nodes, parent, nxt = list(range(n)), -np.ones(2 * n - 1, int), n
    while len(nodes) > 1:
        i, j = rng.choice(len(nodes), 2, replace=False)
        a, b = nodes[i], nodes[j]
        parent[a] = parent[b] = nxt
        nodes = [x for q, x in enumerate(nodes) if q not in (i, j)] + [nxt]
        nxt += 1
    return parent


def random_params(n, k, N, M, rng):
    return ref.Params(
        N=N, Lam_T=rng.uniform(0.5, 4.0), r=rng.uniform(0.5, 1.5, k),
        rate=ref.PWShape([0, rng.uniform(.2, .8), 1], rng.uniform(.1, 2, 2)),
        comp_knots=np.array([0, rng.uniform(.2, .8), 1]), xi=rng.dirichlet(np.ones(M), 2),
        Gam_T=rng.uniform(0, 0.6), sil=ref.PWShape([0, rng.uniform(.2, .8), 1], rng.uniform(.1, 2, 2)),
        alpha=rng.normal(0, .7, n), beta=rng.normal(-.5, .7, k))


def sample_reads_full(parent, t, n, p, M, rng):
    """Draw reads for every tape from the FULL-space kernels (independent of the slot code)."""
    states = ref.full_states(p.N, M)
    nn, k = len(parent), len(p.r)
    reads = [[None] * k for _ in range(n)]
    for z in range(k):
        X = [None] * nn
        Kst = ref.full_kernel(states, ref.branch_epochs(0.0, t[nn - 1], z, p), p.N)
        X[nn - 1] = rng.choice(len(states), p=Kst[0] / Kst[0].sum())
        for v in range(nn - 2, -1, -1):
            Kv = ref.full_kernel(states, ref.branch_epochs(t[parent[v]], t[v], z, p), p.N)
            row = Kv[X[parent[v]]]
            X[v] = rng.choice(len(states), p=row / row.sum())
        for c in range(n):
            s = states[X[c]]
            pi = ref.sigmoid(p.alpha[c] + p.beta[z])
            if s == 'S' or len(s) == 0 or rng.random() < pi:
                reads[c][z] = None
            else:
                reads[c][z] = s
    return reads


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ===================================================================================================
# V1, V2, V4 (twin only)
# ===================================================================================================

def check_V1(rng, quick):
    N, M, worst_full, worst_enum, n_full, n_enum = 3, 3, 0.0, 0.0, 0, 0
    trials = 20 if quick else 100
    for tr in range(trials):
        n = int(rng.integers(3, 7))
        parent = random_tree(n, rng)
        t = ref.times_from_x(rng.normal(0, 1.2, n - 1), parent, n)
        p = random_params(n, 2, N, M, rng)
        sampled = sample_reads_full(parent, t, n, p, M, rng)
        patterns = [sampled]
        if tr % 2 == 0:                                   # a whole subtree missing on tape 0
            sub = [c for c in range(n) if parent[c] == parent[0]]
            patterns.append([[None if (c in sub and z == 0) else sampled[c][z] for z in range(2)]
                             for c in range(n)])
        if tr % 5 == 0:
            patterns.append([[None, None] for _ in range(n)])
        for reads in patterns:
            for z in range(2):
                rz = [reads[c][z] for c in range(n)]
                a = ref.loglik_tape(parent, t, n, rz, z, p)
                b = ref.loglik_tape_full(parent, t, n, rz, z, p, M)
                worst_full = max(worst_full, rel(a, b)); n_full += 1
                if n - 1 <= 3:
                    e = ref.loglik_tape_enumerate(parent, t, n, rz, z, p, M)
                    worst_enum = max(worst_enum, rel(a, e)); n_enum += 1
    return {"vs_full_space": {"configs": n_full, "worst_rel": worst_full},
            "vs_enumeration": {"configs": n_enum, "worst_rel": worst_enum}}


def check_V2(rng, quick):
    out = []
    for (N, M, n) in [(2, 2, 3), (2, 2, 4), (3, 2, 3)]:
        options = [None] + [s for j in range(1, N + 1) for s in itertools.product(range(M), repeat=j)]
        reps = 3 if quick else 10
        for _ in range(reps):
            parent = random_tree(n, rng)
            t = ref.times_from_x(rng.normal(0, 1.2, n - 1), parent, n)
            p = random_params(n, 1, N, M, rng)
            tot = 0.0
            for combo in itertools.product(options, repeat=n):
                tot += math.exp(ref.loglik_tape(parent, t, n, list(combo), 0, p))
            out.append({"N": N, "M": M, "n": n, "sum": tot, "abs_dev": abs(tot - 1.0)})
    return {"configs": len(out), "worst_abs": max(o["abs_dev"] for o in out), "rows": out}


def generator(mu, G, xi_next, anchored, N):
    S = N + 1
    Q = np.zeros((N + 2, N + 2))
    for d in range(N):
        Q[d, d + 1] = mu * (xi_next[d] if anchored else 1.0)
        Q[d, d] -= mu
    for d in range(N + 1):
        if not anchored:
            Q[d, S] += G
        Q[d, d] -= G
    return Q


def check_V4(rng, quick):
    N, worst, cnt = N6, 0.0, 0
    for _ in range(200 if quick else 1000):
        M = 5
        xi = rng.dirichlet(np.ones(M))
        mu, G = rng.uniform(0, 4), rng.uniform(0, 1)
        anchored = rng.random() < 0.5
        lcp = tuple(rng.integers(0, M, N)) if anchored else None
        Kc = ref.segment_kernel(mu, G, xi, lcp, N)
        xi_next = np.array([xi[s] for s in lcp]) if anchored else np.ones(N)
        Ke = expm(generator(mu, G, xi_next, anchored, N))
        cols = slice(0, N + 1) if anchored else slice(0, N + 2)
        worst = max(worst, np.abs(Kc[:N + 1, cols] - Ke[:N + 1, cols]).max()); cnt += 1
    return {"configs": cnt, "worst_abs": worst}


# ===================================================================================================
# simulator-driven checks (V3, V5, V7)
# ===================================================================================================

TV_RATE = dict(rate_knots=(0.0, 0.25, 0.35, 0.6, 1.0), rate_levels=(1.2, 0.0, 1.6, 0.8))
CONST = dict(rate_knots=(0.0, 1.0), rate_levels=(1.0,))
SETTINGS = {
    "const":       dict(rate=CONST, comp=None, Hbar=0.0, sil=((0, 1), (1,))),
    "sil":         dict(rate=CONST, comp=None, Hbar=0.10, sil=((0, 1), (1,))),
    "sil_early":   dict(rate=CONST, comp=None, Hbar=0.10, sil=((0, .4, 1), (2, 1 / 3))),
    "tv_rate_comp": dict(rate=TV_RATE, comp=((0, .4, .7, 1), ((.02,), (.3,), (.02,))),
                         Hbar=0.10, sil=((0, .6, 1), (1 / 3, 2))),
}
PANEL = dl.PanelParams(phi=0.15, r_cl=0.6, s_r=0.3, delta_beta=2.0)
MASK = dl.MaskParams(pbar=0.3, s_beta=1.0, s_alpha=0.7)


def make_recorder(st, r, Lam_T=5.5):
    kw = dict(k=K_TAPES, Lam_T=Lam_T, r=r, **st["rate"])
    if st["comp"] is None:
        kw["shares"] = ((0.0,),)
    else:
        kw["comp_knots"], kw["shares"] = st["comp"]
    return ed.Recorder(**kw)


def simulate(parent, branch, n, st, seed):
    rng = np.random.default_rng(seed)
    pdraw = dl.panel_draws(K_TAPES, rng)
    u, r = dl.panel_tapes(pdraw, PANEL)
    beta, _ = dl.panel_beta(pdraw, u, PANEL, MASK)
    rec = make_recorder(st, r)
    D, sym, out, _ = dl.decorate_tips(parent, branch, n, rec, rng)
    md = dl.mask_draws(n, K_TAPES, rng)
    X, _ = dl.technical_mask(md, beta, MASK.s_alpha)
    sp = sl.SilencingParams(st["Hbar"], st["sil"][0], st["sil"][1])
    so = sl.silence(parent, out["t"], sl.loss_draws(len(parent), K_TAPES, rng), sp)
    Y = sl.observe(X, D, so["H"][:n])
    return dict(t=out["t"], D=D, sym=sym, Y=Y, alpha=MASK.s_alpha * md["a"], beta=beta, r=r,
                rec=rec, sp=sp)


def engine_inputs(sim, parent, n):
    rec, sp = sim["rec"], sim["sp"]
    ck, xi = eng.comp_epochs_from_recorder(rec)
    prm = dict(Lam_T=jnp.asarray(rec.Lam_T), r=jnp.asarray(rec.r), rate_levels=jnp.asarray(rec.omega),
               logxi=jnp.asarray(eng.safe_logxi(xi)), Gam_T=jnp.asarray(sp.Gam_T),
               sil_levels=jnp.asarray(sl.Shape(sp.knots, sp.levels).omega),
               alpha=jnp.asarray(sim["alpha"]), beta=jnp.asarray(sim["beta"]),
               m=jnp.ones(len(parent)))
    structure = (N6, tuple(rec.kappa), tuple(ck), tuple(sl.Shape(sp.knots, sp.levels).kappa))
    data = eng.prepare(parent, n, sim["D"], sim["sym"], sim["Y"], N6)
    twin = ref.Params(N=N6, Lam_T=rec.Lam_T, r=rec.r, rate=ref.PWShape(rec.kappa, rec.omega),
                      comp_knots=ck, xi=xi, Gam_T=sp.Gam_T, sil=ref.PWShape(sp.knots, sp.levels),
                      alpha=sim["alpha"], beta=sim["beta"])
    return prm, structure, data, twin


_LL = {}


def loglik_fn(structure):
    if structure not in _LL:
        _LL[structure] = eng.make_loglik(*structure)
    return _LL[structure]


def check_V3(seed, quick):
    worst, rows, cnt = 0.0, [], 0
    sizes = (25, 64) if quick else (25, 64, 210)
    per = 1 if quick else 5
    names = list(SETTINGS)
    for n in sizes:
        for rho in (0.002, 0.1, 0.25):
            for th in (0.3, 0.7):
                for i, (parent, branch) in enumerate(ed.load_library_trees(n, rho, th)[:per]):
                    st = names[(i + cnt) % len(names)]
                    sim = simulate(parent, branch, n, SETTINGS[st], [seed, 3, n, int(rho * 1e4), i])
                    prm, structure, data, twin = engine_inputs(sim, parent, n)
                    lj = np.asarray(loglik_fn(structure)(jnp.asarray(sim["t"]), prm, data))
                    reads = ref.reads_from_sim(sim["D"], sim["sym"], sim["Y"])
                    lr = ref.loglik(parent, sim["t"], n, reads, twin)
                    w = float(np.max(np.abs(lj - lr) / np.abs(lr)))
                    worst = max(worst, w); cnt += 1
                    rows.append({"n": n, "rho": rho, "theta": th, "tree": i, "setting": st, "worst_rel": w})
    return {"trees": cnt, "tape_comparisons": cnt * K_TAPES, "worst_rel": worst, "rows": rows}


def check_V5(seed, quick):
    worst, rows = 0.0, []
    ntrees = 4 if quick else 20
    cells = [(25, 0.1, 0.3), (64, 0.25, 0.7), (25, 0.002, 0.7), (64, 0.1, 0.3)]
    names = list(SETTINGS)
    for i in range(ntrees):
        n, rho, th = cells[i % len(cells)]
        parent, branch = ed.load_library_trees(n, rho, th)[i // len(cells)]
        st = names[i % len(names)]
        sim = simulate(parent, branch, n, SETTINGS[st], [seed, 5, i])
        prm, structure, data, _ = engine_inputs(sim, parent, n)
        b, d, _ = bd.rate_params(n, rho, th)
        J = eng.make_objective(loglik_fn(structure))
        x0 = jnp.asarray(ref.x_from_times(sim["t"], parent, n))
        theta0 = dict(prm, b=jnp.asarray(b), d=jnp.asarray(d))

        def f(x, th_):
            pr = {kk: vv for kk, vv in th_.items() if kk not in ("b", "d")}
            return J(x, pr, data, (th_["b"], th_["d"], rho))

        gx, gth = jax.grad(f, argnums=(0, 1))(x0, theta0)
        h = 1e-5
        checks = [("x", j) for j in range(min(len(x0), 12))] + \
                 [("Lam_T", None), ("Gam_T", None), ("b", None), ("d", None), ("beta", 0), ("alpha", 0)]
        for name, idx in checks:
            def shifted(eps):
                if name == "x":
                    return f(x0.at[idx].add(eps), theta0)
                th2 = dict(theta0)
                th2[name] = theta0[name] + eps if idx is None else theta0[name].at[idx].add(eps)
                return f(x0, th2)
            fd = float((shifted(h) - shifted(-h)) / (2 * h))
            ad = float(gx[idx]) if name == "x" else float(gth[name] if idx is None else gth[name][idx])
            e = abs(ad - fd) / max(abs(fd), 1.0)
            worst = max(worst, e)
            rows.append({"tree": i, "param": name, "idx": idx, "autodiff": ad, "fd": fd, "err": e})
    return {"trees": ntrees, "gradients": len(rows), "worst": worst, "rows": rows}


def check_V7(seed, quick):
    """Mean score at the truth over datasets simulated on library trees (n = 210)."""
    ndata = 40 if quick else 200
    cells = [(rho, th) for rho in (0.0005, 0.002, 0.02, 0.1, 0.25) for th in (0.3, 0.7)]
    out = {}
    for st in ("sil", "tv_rate_comp"):
        scores = []
        for i in range(ndata):
            rho, th = cells[i % len(cells)]
            parent, branch = ed.load_library_trees(210, rho, th)[i // len(cells)]
            sim = simulate(parent, branch, 210, SETTINGS[st], [seed, 7, i, len(st)])
            prm, structure, data, _ = engine_inputs(sim, parent, 210)
            ll = loglik_fn(structure)
            t = jnp.asarray(sim["t"])
            gt, gp = jax.grad(lambda tt, pp: jnp.sum(ll(tt, pp, data)), argnums=(0, 1))(t, prm)
            scores.append({
                "dlogLam_T": float(gp["Lam_T"] * prm["Lam_T"]),
                "dGam_T": float(gp["Gam_T"]) if st != "const" else 0.0,
                "dbeta_shift": float(jnp.sum(gp["beta"])),
                "dalpha_shift": float(jnp.sum(gp["alpha"])),
                "dt_root": float(gt[-1]),
                "dt_internal_sum": float(jnp.sum(gt[210:]))})
        for key in scores[0]:
            v = np.array([s[key] for s in scores])
            z = float(v.mean() / (v.std(ddof=1) / math.sqrt(v.size))) if v.std() > 0 else 0.0
            out[f"{st}:{key}"] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)), "n": int(v.size), "z": z}
    return out


def check_V8(quick):
    out = {}
    for n in ((25, 64) if quick else (25, 64, 210)):
        for rho in (0.002, 0.1, 0.25):
            for th in (0.3, 0.7):
                trees = ed.load_library_trees(n, rho, th)
                b, d, _ = bd.rate_params(n, rho, th)
                dec, rootp = [], []
                for parent, branch in trees:
                    t = ed.node_times(parent, branch, n)
                    h = 1.0 - t[n:]
                    u = bd.bds_G(h, b, d, rho, 1.0)                  # prior CDF of one branching height
                    dec.append(np.histogram(u, bins=10, range=(0, 1))[0] / u.size - 0.1)
                    rootp.append(float(bd.bds_G(h.max(), b, d, rho, 1.0)) ** (n - 1) - 0.5)
                dec, rootp = np.array(dec), np.array(rootp)
                for q in range(10):
                    v = dec[:, q]
                    out[f"n{n}_rho{rho}_th{th}:decile{q}"] = float(v.mean() / (v.std(ddof=1) / math.sqrt(v.size)))
                out[f"n{n}_rho{rho}_th{th}:root_pit"] = float(rootp.mean() / (rootp.std(ddof=1) / math.sqrt(rootp.size)))
    return out


# ===================================================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=6102026)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    res, t0 = {}, time.perf_counter()
    for name, fn in (("V1", lambda: check_V1(rng, a.quick)), ("V2", lambda: check_V2(rng, a.quick)),
                     ("V4", lambda: check_V4(rng, a.quick)), ("V3", lambda: check_V3(a.seed, a.quick)),
                     ("V5", lambda: check_V5(a.seed, a.quick)), ("V7", lambda: check_V7(a.seed, a.quick)),
                     ("V8", lambda: check_V8(a.quick))):
        tc = time.perf_counter()
        res[name] = fn()
        print(f"{name} done in {time.perf_counter() - tc:.0f} s", flush=True)
    hard = {
        "V1 twin vs full state space": (res["V1"]["vs_full_space"]["configs"], res["V1"]["vs_full_space"]["worst_rel"], 1e-10),
        "V1 twin vs enumeration": (res["V1"]["vs_enumeration"]["configs"], res["V1"]["vs_enumeration"]["worst_rel"], 1e-10),
        "V2 sum over all reads = 1": (res["V2"]["configs"], res["V2"]["worst_abs"], 1e-10),
        "V3 JAX vs twin, per-tape ell": (res["V3"]["tape_comparisons"], res["V3"]["worst_rel"], 1e-10),
        "V4 closed form vs expm": (res["V4"]["configs"], res["V4"]["worst_abs"], 1e-12),
        "V5 autodiff vs finite diff": (res["V5"]["gradients"], res["V5"]["worst"], 1e-5),
    }
    stat = {"V7": {k: v["z"] for k, v in res["V7"].items()}, "V8": res["V8"]}
    hard_ok = all(w <= tol for (_, w, tol) in hard.values())
    zs = [abs(z) for d in stat.values() for z in d.values()]
    trips = [(c, k, z) for c, d in stat.items() for k, z in d.items() if abs(z) > 3.0]
    verdict = "PASS" if hard_ok and not trips else (
        "RERUN (single statistical trip)" if hard_ok and len(trips) == 1 else "FAIL")
    summary = {"hard": {k: {"configs": c, "worst": w, "tol": tol, "pass": w <= tol} for k, (c, w, tol) in hard.items()},
               "stat_cells": len(zs), "stat_worst_abs_z": max(zs), "stat_trips": trips,
               "expected_trips_by_chance": 0.0027 * len(zs), "verdict": verdict,
               "seed": a.seed, "quick": a.quick, "seconds": round(time.perf_counter() - t0, 1)}
    RES.mkdir(parents=True, exist_ok=True)
    fname = "validate_engine_quick.json" if a.quick else "validate_engine.json"
    (RES / fname).write_text(json.dumps({"summary": summary, "detail": res}, indent=1, default=float))
    print(f"\n{'check':34s} {'configs':>8s} {'worst':>10s} {'tol':>8s}  verdict")
    for k, (c, w, tol) in hard.items():
        print(f"{k:34s} {c:8d} {w:10.2e} {tol:8.0e}  {'PASS' if w <= tol else 'FAIL'}")
    for c, d in stat.items():
        w = max(abs(z) for z in d.values())
        print(f"{c + ' (statistical, |z|)':34s} {len(d):8d} {w:10.2f} {3.0:8.1f}  {'PASS' if w <= 3 else 'TRIP'}")
    print(f"\nstatistical cells {len(zs)}, trips {len(trips)} (chance expects {0.0027 * len(zs):.2f}); "
          f"verdict: {verdict}; {summary['seconds']} s")


if __name__ == "__main__":
    main()
