#!/usr/bin/env python3
r"""
15_validate_editing.py -- validate the editing layer (14_editing.py) against the closed forms of
notes/sciphy_notes.md §S4.9.5.  Approved 2026-09-27 ("PROPOSAL (2026-09-24c)", ../CLAUDE.md).

A failed check is a BUG, not a finding.  Every entry is a deviation from a closed form in
standard errors, the se taken ACROSS TREES (cells of one tree share ancestry, so never across
cells or pairs).  The yardstick is calibrated two ways: a half-split of the trees for the
closed-form checks, and route (ii) against itself at a second seed for the route comparison.

CONFIGURATIONS (all r_z = 1 unless stated; N = 6, k = 30, Lam_T = 5.5, Lam_pre = 0)
  CONST  lambda_0 = 1; one signal channel of one symbol at p = 0.05; xi^L = Park mice pooled.
  TV     lambda_0 on knots (0, .25, .35, .6, 1) with levels (1.2, 0, 1.6, 0.8) -- a ZERO-rate
         window at .25-.35 -- and p = .02 / .30 / .02 on knots (0, .4, .7, 1), deliberately not
         aligned with the rate knots.  This is the setting in which check B has teeth.
  RHET   CONST with r_z ~ lognormal(0, 0.5) rescaled to mean 1.  ⚠ ADDED beyond the proposal: the
         per-tape speed is a code path the proposal's r = 1 runs never exercise.  Checks A_r and
         A' only.

CHECKS
  A    tip depth pmf vs min(Poisson(r_z Lam_T), N)                      (tree-blind)
  A'   internal-node depth pmf vs min(Poisson(Lam_z(t_node)), N), in 3 node-time bins
  B    share of the signal channel at slot j among tips of final depth d, vs
       int K_{j,d}(x) p(Lam^{-1}(x)) dx  -- the only check that sees edit ORDER
  B_L  within-lineage composition, pooled over slots, in 10 equal-mass groups of xi^L
  C    route (i) vs route (ii), paired per tree: tip depth pmf, mean edit time by slot,
       signal share by slot
  D    hard assertions: hand-off slot == depth recomputed from the path; depth <= N; edit times
       inside their branch; times non-decreasing along every tape; symbol ids valid; the
       generator's own depth == the recomputed depth; prefix identical up to the LCA's depth
  E    collision at the first slot past the LCA vs <xi(t1), xi(t2)> on the realised times,
       binned by LCA depth (and, under TV, by how many of the two writes fell in the on-window),
       over random TIP PAIRS -- as proposed
  E_u  ⚠ ADDED after the 20-tree smoke run: the same statistic over UNIQUE comparisons.  Most tip
       pairs compare the same two edits (a slot written just below the LCA is shared by every tip
       beneath it), so E rests on a handful of independent rare events per tree and its
       t-statistic is not calibrated.  E_u keys each comparison by (branch that wrote edit 1,
       branch that wrote edit 2, tape) and counts it once.  The weighting depends on structure and
       times only, never on symbols, so both are unbiased; E_u is the one whose se and
       minimum-count rule are honest.
  A_r  (RHET) tip mean depth per tape vs E[min(Poisson(r_z Lam_T), N)], 5 r-quantile bins

READING RULE (proposal): PASS = worst |z| <= 3.0 over every tested cell AND zero D violations.
Fixed BEFORE the 200-tree run: "every tested cell" includes the ADDED checks (A_r, E_u) too.
A cell is tested only if its expected positive and negative counts, summed over trees, are each
>= 100 (otherwise its se is not estimable); excluded cells are reported, not hidden.

usage:  15_validate_editing.py [--seed S] [--quick]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
import time

import numpy as np
from scipy.special import gammainc, gammaln
from scipy.stats import beta as beta_dist
from scipy.stats import poisson

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("editing", _HERE / "14_editing.py")
ed = importlib.util.module_from_spec(_spec)
sys.modules["editing"] = ed
_spec.loader.exec_module(ed)

RES = _HERE.parent / "results"
N_TIPS = 210
TREE_CELLS = [(rho, th) for rho in (0.0005, 0.002, 0.02, 0.1, 0.25) for th in (0.3, 0.7)]
BIG = (1976, 0.25, 0.5)            # cost + D at the largest design size
PAIRS_PER_TREE = 3000
MIN_EXPECTED = 100.0
Z_PASS = 3.0
N_LGROUPS = 10

CONFIGS = {
    "CONST": dict(rate_knots=(0.0, 1.0), rate_levels=(1.0,), comp_knots=(0.0, 1.0),
                  shares=((0.05,),)),
    "TV": dict(rate_knots=(0.0, 0.25, 0.35, 0.6, 1.0), rate_levels=(1.2, 0.0, 1.6, 0.8),
               comp_knots=(0.0, 0.4, 0.7, 1.0), shares=((0.02,), (0.30,), (0.02,))),
}


def trunc_pmf(mu, N):
    """P(min(Poisson(mu), N) = d), d = 0..N; shape mu.shape + (N+1,)."""
    mu = np.asarray(mu, float)
    p = poisson.pmf(np.arange(N), mu[..., None])
    return np.concatenate([p, 1.0 - p.sum(-1, keepdims=True)], axis=-1)


def kernel_signal_share(rec, d, j, grid=200_001):
    """int K_{j,d}(x) p~(x) dx at r = 1, mu = Lam_T (§S4.9.5 B)."""
    mu, N = rec.Lam_T, rec.N
    x = np.linspace(0.0, mu, grid)[1:-1]
    ptil = rec.share_at(rec.W_inv(x / rec.Lam_T)).sum(-1)
    if d < N:
        dens = beta_dist.pdf(x / mu, j, d + 1 - j)
    else:
        dens = np.exp((j - 1) * np.log(x) - x - gammaln(j))
        if N - j > 0:
            dens = dens * gammainc(N - j, mu - x)
    return float(np.trapezoid(dens * ptil, x) / np.trapezoid(dens, x))


class Acc:
    """Per-cell, per-tree record: (residual, observed, expected, expected+, expected-)."""

    def __init__(self):
        self.cells: dict = {}

    def add(self, key, tree, resid, obs, exp, epos, eneg):
        self.cells.setdefault(key, {})[tree] = (float(resid), float(obs), float(exp),
                                                float(epos), float(eneg))

    def summarise(self, split_perm=None):
        rows = []
        for key, per in self.cells.items():
            trees = sorted(per)
            v = np.array([per[t] for t in trees])
            epos, eneg = v[:, 3].sum(), v[:, 4].sum()
            tested = epos >= MIN_EXPECTED and eneg >= MIN_EXPECTED and len(trees) >= 20
            r = v[:, 0]
            z = float(r.mean() / (r.std(ddof=1) / np.sqrt(r.size))) if tested and r.std() > 0 \
                else (0.0 if tested else np.nan)
            znull = np.nan
            if tested and split_perm is not None:
                pos = {t: i for i, t in enumerate(split_perm)}
                h = np.array([pos.get(t, -1) % 2 for t in trees])
                a, b = r[h == 0], r[h == 1]
                if a.size > 5 and b.size > 5:
                    znull = float((a.mean() - b.mean()) /
                                  np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size))
            rows.append({"key": list(key), "n_trees": len(trees), "tested": bool(tested),
                         "obs": float(v[:, 1].mean()), "exp": float(v[:, 2].mean()),
                         "z": z, "z_null_split": znull,
                         "expected_pos": float(epos), "expected_neg": float(eneg)})
        return rows


def ancestors(parent):
    nn = parent.size
    anc = np.zeros((nn, nn), bool)
    for x in range(nn - 1, -1, -1):
        if parent[x] >= 0:
            anc[x] = anc[parent[x]]
        anc[x, x] = True
    return anc


def structural(out, rec, pd, tsym, ttime):
    """Check D on one decorated tree: (violations, comparisons) per assertion."""
    parent, t, t0 = out["parent"], out["t"], out["t0"]
    nn, k, N = out["bsym"].shape
    root = nn - 1
    v = {}
    exp_slot0 = np.where((np.arange(nn) == root)[:, None], 0, pd[np.maximum(parent, 0)])
    v["handoff_slot"] = (int(np.sum(out["bslot0"] != exp_slot0)), nn * k)
    v["depth_le_N"] = (int(np.sum(pd > N)), nn * k)
    v["gen_depth_eq_path"] = (int(np.sum(out["gen_depth"] != pd)), nn * k)
    used = np.arange(N)[None, None, :] < out["bcount"][..., None]
    bt = out["btime"].astype(np.float64)
    tol = 1e-6
    inside = (bt >= t0[:, None, None] - tol) & (bt <= t[:, None, None] + tol)
    v["time_in_branch"] = (int(np.sum(used & ~inside)), int(used.sum()))
    v["time_nan_iff_unused"] = (int(np.sum(np.isnan(bt) == used)), nn * k * N)
    sym = out["bsym"]
    v["symbol_valid"] = (int(np.sum(used & ((sym < 0) | (sym >= rec.M))) +
                             np.sum(~used & (sym != -1))), nn * k * N)
    tt = ttime.astype(np.float64)
    both = np.isfinite(tt[..., 1:]) & np.isfinite(tt[..., :-1])
    v["queue_order"] = (int(np.sum(both & (tt[..., 1:] < tt[..., :-1]))), int(both.sum()))
    return v


def pair_checks(out, rec, pd, tsym, ttime, torig, rng, acc_E, tree, tv):
    """E (collision past the LCA) and D's prefix assertion, on random tip pairs."""
    n, parent, t = out["n"], out["parent"], out["t"]
    k, N = rec.k, rec.N
    anc = ancestors(parent)
    c1 = rng.integers(0, n, PAIRS_PER_TREE)
    c2 = (c1 + rng.integers(1, n, PAIRS_PER_TREE)) % n
    common = anc[c1] & anc[c2]
    lca = np.argmax(np.where(common, t[None, :], -1.0), axis=1)
    DA = pd[lca]                                                   # (P, k)
    # D: prefix identical in slots < D_A
    slot = np.arange(N)[None, None, :]
    pre = slot < DA[..., None]
    mism = int(np.sum(pre & (tsym[c1] != tsym[c2])))
    ncmp = int(pre.sum())
    # E: first slot past the LCA
    elig = (DA < N) & (pd[c1] > DA) & (pd[c2] > DA)
    P_idx, z_idx = np.nonzero(elig)
    s = DA[P_idx, z_idx]
    y1, y2 = tsym[c1[P_idx], z_idx, s], tsym[c2[P_idx], z_idx, s]
    t1 = ttime[c1[P_idx], z_idx, s].astype(float)
    t2 = ttime[c2[P_idx], z_idx, s].astype(float)
    match = (y1 == y2).astype(float)
    expv = rec.inner(t1, t2)
    on1 = rec.share_at(t1).sum(-1) > 0.1
    on2 = rec.share_at(t2).sum(-1) > 0.1
    cls = on1.astype(int) + on2.astype(int)
    # E_u: one entry per unique (origin of edit 1, origin of edit 2, tape)
    o1, o2 = torig[c1[P_idx], z_idx, s], torig[c2[P_idx], z_idx, s]
    lo = np.minimum(o1, o2).astype(np.int64)
    hi = np.maximum(o1, o2).astype(np.int64)
    nn = parent.size
    ukey = (lo * nn + hi) * k + z_idx
    _, first = np.unique(ukey, return_index=True)
    uniq = np.zeros(ukey.size, bool)
    uniq[first] = True
    for tag, sel in (("E", np.ones(ukey.size, bool)), ("E_u", uniq)):
        bins = {(tag, "lca_depth", int(b)): sel & (s == b) for b in range(N)}
        if tv:
            for c in range(3):
                bins[(tag, "writes_in_on_window", c)] = sel & (cls == c)
        for key, m in bins.items():
            if m.sum() == 0:
                continue
            acc_E.add(key, tree, np.mean(match[m] - expv[m]), match[m].mean(), expv[m].mean(),
                      expv[m].sum(), (1 - expv[m]).sum())
    return mism, ncmp


def tip_stats(pd, tsym, ttime, rec, n):
    """Per-tree tip statistics shared by A, B, B_L and C."""
    N = rec.N
    tip_d, tip_s, tip_t = pd[:n], tsym[:n], ttime[:n].astype(float)
    sig = (tip_s >= 0) & (rec.channel_of[np.maximum(tip_s, 0)] < rec.A)
    out = {}
    for d in range(N + 1):
        out[("depth_pmf", d)] = ((tip_d == d).mean(), int((tip_d == d).sum()), tip_d.size)
    for j in range(N):
        m = tip_d > j
        if m.sum():
            out[("slot_time", j + 1)] = (tip_t[..., j][m].mean(), int(m.sum()), int(m.sum()))
            out[("slot_signal", j + 1)] = (sig[..., j][m].mean(), int(sig[..., j][m].sum()),
                                           int(m.sum()))
    return out, tip_d, sig


def run_config(name, cfg, trees, seed, acc, cacc, rhet=False):
    rng_master = np.random.SeedSequence([seed, {"CONST": 1, "TV": 2, "RHET": 3}[name]])
    ss = rng_master.spawn(len(trees))
    r = None
    if rhet:
        rr = np.random.default_rng(np.random.SeedSequence([seed, 7])).lognormal(0, 0.5, 30)
        r = rr / rr.mean()
    rec = ed.Recorder(k=30, r=r, **cfg)
    N, k = rec.N, rec.k
    tv = name == "TV"
    D_tot: dict = {}
    prefix = [0, 0]
    secs = []
    EB = {}
    if not rhet:
        for d in range(1, N + 1):
            for j in range(1, d + 1):
                EB[(d, j)] = kernel_signal_share(rec, d, j)
    gl = np.minimum(((np.cumsum(rec.xi_L) - rec.xi_L / 2) * N_LGROUPS).astype(int), N_LGROUPS - 1)
    gmass = np.bincount(gl, weights=rec.xi_L, minlength=N_LGROUPS)
    rbins = np.digitize(rec.r, np.quantile(rec.r, [0.2, 0.4, 0.6, 0.8])) if rhet else None
    for ti, (parent, branch) in enumerate(trees):
        n = N_TIPS
        g_main, g_i, g_ii2, g_pairs = [np.random.default_rng(s) for s in ss[ti].spawn(4)]
        t_s = time.perf_counter()
        out = ed.decorate(parent, branch, n, rec, g_main, route="ii")
        secs.append(time.perf_counter() - t_s)
        pd, tsym, ttime, torig = ed.assemble(out)
        for kk, (viol, cmp_) in structural(out, rec, pd, tsym, ttime).items():
            a = D_tot.setdefault(kk, [0, 0]); a[0] += viol; a[1] += cmp_
        tips, tip_d, sig = tip_stats(pd, tsym, ttime, rec, n)

        # A: tip depth pmf
        mu_tip = rec.Lam_T * rec.r                                     # (k,)
        expA = trunc_pmf(mu_tip, N).mean(0)                            # (N+1,)
        nt = n * k
        for d in range(N + 1):
            o = tips[("depth_pmf", d)][0]
            acc.add(("A", "tip_depth", d), ti, o - expA[d], o, expA[d], expA[d] * nt,
                    (1 - expA[d]) * nt)
        # A': internal nodes, 3 time bins
        internal = np.arange(n, parent.size)
        tv_nodes = out["t"][internal]
        mu_nodes = rec.Lam(tv_nodes)                                   # (m, k)
        pm = trunc_pmf(mu_nodes, N)                                    # (m, k, N+1)
        Dn = pd[internal]
        tb = np.minimum((tv_nodes * 3).astype(int), 2)
        for b in range(3):
            m = tb == b
            if not m.any():
                continue
            cnt = m.sum() * k
            for d in range(N + 1):
                o = (Dn[m] == d).mean()
                e_ = pm[m][..., d].mean()
                acc.add(("A'", f"node_time_bin{b}", d), ti, o - e_, o, e_, e_ * cnt,
                        (1 - e_) * cnt)
        if rhet:
            # A_r: tip mean depth per tape vs E[min(Pois(r Lam_T), N)], by r-quantile bin
            md = (np.minimum(np.arange(80), N)[None, :] *
                  poisson.pmf(np.arange(80)[None, :], mu_tip[:, None])).sum(1)
            for b in range(5):
                m = rbins == b
                o = tip_d[:, m].mean()
                e_ = md[m].mean()
                acc.add(("A_r", "r_quintile", b), ti, o - e_, o, e_, 1e9, 1e9)
        else:
            # B: signal share at slot j among tips of final depth d
            for (d, j), e_ in EB.items():
                m = tip_d == d
                cnt = int(m.sum())
                if cnt == 0:
                    continue
                o = sig[..., j - 1][m].mean()
                acc.add(("B", f"d{d}", j), ti, o - e_, o, e_, e_ * cnt, (1 - e_) * cnt)
            # B_L: within-lineage composition, pooled over slots
            ts = tsym[:n]
            lin = (ts >= 0) & (rec.channel_of[np.maximum(ts, 0)] == rec.A)
            idx = ts[lin].astype(int) - rec.offsets[rec.A]
            if idx.size:
                frac = np.bincount(gl[idx], minlength=N_LGROUPS) / idx.size
                for g in range(N_LGROUPS):
                    acc.add(("B_L", "xiL_group", g), ti, frac[g] - gmass[g], frac[g], gmass[g],
                            gmass[g] * idx.size, (1 - gmass[g]) * idx.size)
            # E + D prefix
            mm, nc = pair_checks(out, rec, pd, tsym, ttime, torig, g_pairs, acc, ti, tv)
            prefix[0] += mm; prefix[1] += nc
            # C: route (i) and route (ii) at a second seed, paired per tree
            for label, gen, route in (("i", g_i, "i"), ("ii_seed2", g_ii2, "ii")):
                o2 = ed.decorate(parent, branch, n, rec, gen, route=route)
                pd2, ts2, tt2, _ = ed.assemble(o2)
                tips2, _, _ = tip_stats(pd2, ts2, tt2, rec, n)
                for key, (val, pos, tot) in tips.items():
                    if key not in tips2:
                        continue
                    val2, pos2, tot2 = tips2[key]
                    ep = min(pos, pos2)
                    en = min(tot - pos, tot2 - pos2) if key[0] != "slot_time" else ep
                    cacc[label].add(("C", key[0], key[1]), ti, val - val2, val, val2, ep, en)
    D_tot["prefix_to_lca"] = prefix
    return rec, D_tot, float(np.mean(secs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=24092026)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t_start = time.perf_counter()
    per_cell = 2 if a.quick else 20
    trees = []
    for rho, th in TREE_CELLS:
        trees += ed.load_library_trees(N_TIPS, rho, th)[:per_cell]
    split = np.random.default_rng(a.seed + 1).permutation(len(trees))
    out = {"settings": {"seed": a.seed, "n": N_TIPS, "trees": len(trees),
                        "tree_cells_rho_theta": TREE_CELLS, "pairs_per_tree": PAIRS_PER_TREE,
                        "min_expected": MIN_EXPECTED, "z_pass": Z_PASS,
                        "configs": {k: {kk: list(map(list, v)) if kk == "shares" else list(v)
                                        for kk, v in c.items()} for k, c in CONFIGS.items()},
                        "rhet": "CONST with r_z ~ lognormal(0, 0.5), mean 1 (added)"},
           "checks": [], "D": {}, "seconds_per_tree_n210": {}}
    all_rows = []
    for name, cfg, rhet in (("CONST", CONFIGS["CONST"], False), ("TV", CONFIGS["TV"], False),
                            ("RHET", CONFIGS["CONST"], True)):
        acc = Acc()
        cacc = {"i": Acc(), "ii_seed2": Acc()}
        tc = time.perf_counter()
        rec, D_tot, spt = run_config(name, cfg, trees, a.seed, acc, cacc, rhet)
        out["seconds_per_tree_n210"][name] = round(spt, 4)
        out["settings"][f"recorder_{name}"] = {k: v for k, v in rec.settings().items() if k != "r"}
        out["settings"][f"recorder_{name}"]["q_L"] = float(rec.q_chan[-1])
        out["settings"][f"recorder_{name}"]["M_lineage"] = int(rec.xi_L.size)
        rows = acc.summarise(split)
        if not rhet:
            ci = {tuple(r_["key"]): r_ for r_ in cacc["i"].summarise()}
            cn = {tuple(r_["key"]): r_ for r_ in cacc["ii_seed2"].summarise()}
            for key, r_ in ci.items():
                r_["z_null_split"] = cn[key]["z"] if key in cn else np.nan
                rows.append(r_)
        for r_ in rows:
            r_["config"] = name
        all_rows += rows
        out["D"][name] = {k: {"violations": v[0], "comparisons": v[1]} for k, v in D_tot.items()}
        print(f"{name}: {len(trees)} trees, {time.perf_counter() - tc:.0f} s", flush=True)

    # cost + D at the largest design size
    n_big, rho_big, th_big = BIG
    pb, bb = ed.load_library_trees(n_big, rho_big, th_big)[0]
    for name in ("CONST", "TV"):
        rec = ed.Recorder(k=30, **CONFIGS[name])
        cost = {}
        for route in ("ii", "i"):
            g = np.random.default_rng(np.random.SeedSequence([a.seed, 99, len(route)]))
            t_s = time.perf_counter()
            ob = ed.decorate(pb, bb, n_big, rec, g, route=route)
            cost[route] = time.perf_counter() - t_s
        t_s = time.perf_counter()
        pdb, tsb, ttb, _ = ed.assemble(ob)
        cost["assemble"] = time.perf_counter() - t_s
        for kk, (viol, cmp_) in structural(ob, rec, pdb, tsb, ttb).items():
            out["D"][name].setdefault(kk + f"_n{n_big}", {"violations": viol, "comparisons": cmp_})
        out[f"seconds_n{n_big}_{name}"] = {k: round(v, 3) for k, v in cost.items()}

    # round trip of the storage format
    rec = ed.Recorder(k=30, **CONFIGS["TV"])
    o = ed.decorate(*trees[0], N_TIPS, rec, np.random.default_rng(a.seed))
    path = RES / "edit_library" / "validation_roundtrip.npz"
    ed.save_decorated(path, o, rec, {"seed": a.seed, "tree": "validation tree 0"})
    back = ed.load_decorated(path)
    rt_ok = all(np.array_equal(back[k], o[k], equal_nan=(k == "btime"))
                for k in ("parent", "bcount", "bslot0", "bsym", "btime")) and \
        np.allclose(back["t"], o["t"])
    out["storage_roundtrip_identical"] = bool(rt_ok)

    # verdict
    tested = [r_ for r_ in all_rows if r_["tested"]]
    zs = np.array([abs(r_["z"]) for r_ in tested])
    zn = np.array([r_["z_null_split"] for r_ in tested if np.isfinite(r_["z_null_split"])])
    trips = [r_ for r_ in tested if abs(r_["z"]) > Z_PASS]
    d_viol = sum(v["violations"] for c in out["D"].values() for v in c.values())
    out["checks"] = all_rows
    out["summary"] = {
        "cells_tested": int(len(tested)), "cells_excluded": int(len(all_rows) - len(tested)),
        "worst_abs_z": float(zs.max()), "cells_over_3": int(len(trips)),
        "expected_trips_by_chance": float(len(tested) * 0.0027),
        "null_z_sd": float(zn.std(ddof=1)), "null_worst_abs_z": float(np.abs(zn).max()),
        "null_cells": int(zn.size), "D_violations_total": int(d_viol),
        "roundtrip_identical": bool(rt_ok),
        "verdict": "PASS" if (zs.max() <= Z_PASS and d_viol == 0 and rt_ok) else
                   ("RERUN (single trip)" if len(trips) == 1 and d_viol == 0 and rt_ok
                    else "FAIL")}
    out["seconds_total"] = round(time.perf_counter() - t_start, 1)
    name = "validate_editing_quick.json" if a.quick else "validate_editing.json"
    (RES / name).write_text(json.dumps(out, indent=1, default=float))

    # the table
    print(f"\n{'check':5s} {'config':6s} {'tested':>6s} {'excl':>5s} {'worst|z|':>9s}  verdict")
    groups: dict = {}
    for r_ in all_rows:
        groups.setdefault((r_["key"][0], r_["config"]), []).append(r_)
    for (chk, cfg), rows in sorted(groups.items()):
        t_ = [x for x in rows if x["tested"]]
        w = max((abs(x["z"]) for x in t_), default=float("nan"))
        print(f"{chk:5s} {cfg:6s} {len(t_):6d} {len(rows) - len(t_):5d} {w:9.2f}  "
              f"{'PASS' if w <= Z_PASS else 'TRIP'}")
    for cfg, dd in out["D"].items():
        tv_ = sum(v["violations"] for v in dd.values())
        tc_ = sum(v["comparisons"] for v in dd.values())
        print(f"D     {cfg:6s}  {tv_} violations of {tc_:.3g} comparisons")
    s = out["summary"]
    print(f"\ncells tested {s['cells_tested']} (excluded {s['cells_excluded']}); worst |z| "
          f"{s['worst_abs_z']:.2f}; cells > 3: {s['cells_over_3']} (chance expects "
          f"{s['expected_trips_by_chance']:.2f}); null z sd {s['null_z_sd']:.2f} over "
          f"{s['null_cells']} cells, null worst {s['null_worst_abs_z']:.2f}; round trip "
          f"{'identical' if rt_ok else 'DIFFERS'}  =>  {s['verdict']}")
    print(f"seconds/tree n={N_TIPS}: {out['seconds_per_tree_n210']}; n={n_big}: "
          f"{out[f'seconds_n{n_big}_CONST']}  total {out['seconds_total']} s")


if __name__ == "__main__":
    main()
