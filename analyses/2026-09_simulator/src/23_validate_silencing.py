#!/usr/bin/env python3
r"""
23_validate_silencing.py -- validate the rung-2 loss layer (22_silencing.py) against closed forms.
Approved 2026-10-01 ("PROPOSAL (2026-10-01, rev. 1)", ../CLAUDE.md).

A failed check is a BUG, not a finding.  Every entry is a deviation from a closed form in standard
errors, the se taken ACROSS TREES; the yardstick is calibrated by a half-split of the trees, and
Acc / split_pairs are 18_validate_dropout.py's, so scoring is identical to rung 1's.

FIXTURES (test cases, not science; n = 210, k = 30, N = 6, 200 library trees as in 18).
  OFF       Hbar 0 on rung-1 BASE                      -- separability: must equal rung 1 bit for bit
  CONST     Hbar 0.25, zeta_0 constant
  EARLY     Hbar 0.25, zeta_0 = 3 on [0,.25), 0 on [.25,.6), 0.625 on [.6,1)   (zero window)
  LATE      Hbar 0.25, zeta_0 = 0 on [0,.7), 10/3 on [.7,1)                    (mostly singletons)
  COMPOSED  rung-1 SLOW (Lam_T 2, pause-and-burst, pbar 0.3, CIS panel) + CONST silencing

CHECKS (Hb = Hbar, W_z = W_zeta, x = the pair's last common ancestor)
  G1  per tape: mean_c H_cz vs Hb                                         CONST/EARLY/LATE
  G2  one tip from each side of every internal node, per tape: H1 H2 vs
      1 - 2(1-Hb) + (1-Hb)^{2 - W_z(t_x)}.  POOLED cell in verdict; 10 W_z-bins secondary
  G3  visible events per tree by clade-size bin vs k sum_e (1-Hb)^{W_z(t0)} [1-(1-Hb)^{dW_z}]
  G4  loss-time PIT  (1-e^{-Gam_T[W_z(t*)-W_z(t0)]}) / (1-e^{-Gam_e})  by decile vs 0.1
  G5  COMPOSED: P(Y=1) per tape vs (1-pi)(1-Hb)(1-e^{-mu_z}); pooled pair Y1 Y2 vs
      (1-pi1)(1-pi2)(1-Hb)^{2-W_z(t_x)} (1 - 2e^{-mu} + e^{-mu(2-W(t_x))})
  diagnostics (outside the verdict; they should MISS, showing the checks have teeth):
  G2n  H1 H2 vs the independence prediction Hb^2                           CONST
  G3x  EARLY's visible events scored against CONST's prediction on the same tree
  H   hard: solve exact; Dollo monotone (no node regains a tape); H at the tips recomputed
      independently as "any L on the root-to-tip path"; sum_e n_e V_ez = sum_c H_cz exactly;
      t* inside its branch and outside every zero-hazard window; no L where Gam_e = 0;
      nested sets across Hbar 0.02 < 0.10 < 0.25 < 0.33 on the same draws; replay identical;
      OFF: H = 0, Y identical to rung 1's observe, exported caches identical array for array;
      COMPOSED: Y = (1-X)(1-H)1[D>=1], exported recovered = Y, term ABSENT wherever H.

READING RULE (fixed before the run): PASS = worst |z| <= 3.0 over every tested verdict cell AND
zero H violations.  Tested = expected positives and negatives each >= 100 summed over trees, and
>= 20 trees; excluded cells are reported.  One trip -> re-run at a new seed; two, or a sign
pattern within one check, is a defect.

usage:  23_validate_silencing.py [--seed S] [--quick]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, _HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


vd = _load("validate_dropout", "18_validate_dropout.py")
sl = _load("silencing", "22_silencing.py")
dl, ed = vd.dl, vd.ed

RES = _HERE.parent / "results"
OUTDIR = RES / "silencing_validation"
N_TIPS, K, N = vd.N_TIPS, vd.K, vd.N
HB = 0.25
SIL = {
    "CONST": sl.SilencingParams(HB),
    "EARLY": sl.SilencingParams(HB, (0.0, 0.25, 0.6, 1.0), (3.0, 0.0, 0.625)),
    "LATE":  sl.SilencingParams(HB, (0.0, 0.7, 1.0), (0.0, 10.0 / 3.0)),
}
NEST = (0.02, 0.10, 0.25, 0.33)
SIZE_BINS = ((1, 1), (2, 2), (3, 4), (5, 16), (17, 64), (65, N_TIPS - 1), (N_TIPS, N_TIPS))


def hv(H, key, viol, cmp_):
    a = H.setdefault(key, [0, 0])
    a[0] += int(viol)
    a[1] += int(cmp_)


def ancestors(parent):
    """anc[x, y] = y is on the root-to-x path (inclusive)."""
    nn = parent.size
    anc = np.zeros((nn, nn), bool)
    for x in range(nn - 1, -1, -1):
        if parent[x] >= 0:
            anc[x] = anc[parent[x]]
        anc[x, x] = True
    return anc


def size_bin_label(lo, hi):
    return f"{lo}" if lo == hi else f"{lo}-{hi}"


def score_tree(name, ti, parent, t, n, sp, out, acc, H, rng_pair, extra=None):
    """G1-G4 + hard assertions for one tree's loss layer."""
    shape = sl.Shape(sp.knots, sp.levels)
    Hb, G = sp.Hbar, sp.Gam_T
    L, V, Hn, tstar, Gam_e = out["L"], out["V"], out["H"], out["tstar"], out["Gam_e"]
    t0, t1 = out["t0"], out["t1"]
    nn = parent.size
    Ht = Hn[:n]
    # ---- hard ---------------------------------------------------------------------------
    nonroot = np.arange(nn - 1)
    hv(H, "dollo_monotone", np.sum(Hn[parent[nonroot]] & ~Hn[nonroot]), (nn - 1) * K)
    anc = ancestors(parent)
    H_alt = (anc[:n].astype(np.int32) @ L.astype(np.int32)) > 0
    hv(H, "H_tips_independent", np.sum(H_alt != Ht), Ht.size)
    sz = sl.clade_sizes(parent, n)
    hv(H, "size_weighted_V_equals_H", np.sum((sz[:, None] * V).sum(0) != Ht.sum(0)), K)
    ee, zz = np.nonzero(L)
    ts = tstar[ee, zz]
    hv(H, "tstar_in_branch", np.sum((ts < t0[ee] - 1e-12) | (ts > t1[ee] + 1e-12)), ts.size)
    for a, b in shape.zero_windows():
        hv(H, "tstar_outside_zero_window", np.sum((ts > a + 1e-12) & (ts < b - 1e-12)), ts.size)
    hv(H, "no_L_where_hazard_0", np.sum(L[Gam_e == 0]), int((Gam_e == 0).sum()) * K)
    # ---- G1 -----------------------------------------------------------------------------
    for z in range(K):
        acc.add(("G1", name, z), ti, Ht[:, z].mean() - Hb, Ht[:, z].mean(), Hb, n * Hb,
                n * (1 - Hb))
    # ---- G2 / G2n -------------------------------------------------------------------------
    nodes, c1, c2 = vd.split_pairs(parent, n, rng_pair)
    Wz = shape.W_of(t[nodes])
    h12 = (Ht[c1] & Ht[c2]).astype(float)
    Eh = np.broadcast_to((1 - 2 * (1 - Hb) + (1 - Hb) ** (2 - Wz))[:, None], h12.shape)
    acc.add(("G2", name, "pooled"), ti, np.mean(h12 - Eh), h12.mean(), Eh.mean(), Eh.sum(),
            (1 - Eh).sum())
    wb = np.minimum((Wz * 10).astype(int), 9)
    for b in range(10):
        m = wb == b
        if m.any():
            acc.add(("G2", name, b), ti, np.mean(h12[m] - Eh[m]), h12[m].mean(), Eh[m].mean(),
                    Eh[m].sum(), (1 - Eh[m]).sum())
    if name == "CONST":
        En = np.full(h12.shape, Hb ** 2)
        acc.add(("G2n", name, "pooled"), ti, np.mean(h12 - En), h12.mean(), Hb ** 2, En.sum(),
                (1 - En).sum())
    # ---- G3 / G3x -------------------------------------------------------------------------
    pv = sl.p_visible(parent, t, sp)
    pvx = sl.p_visible(parent, t, SIL["CONST"]) if name == "EARLY" else None
    for lo, hi in SIZE_BINS:
        m = (sz >= lo) & (sz <= hi)
        if not m.any():
            continue
        obs = float(V[m].sum())
        opp = float(m.sum() * K)
        lab = size_bin_label(lo, hi)
        ex = K * float(pv[m].sum())
        acc.add(("G3", name, lab), ti, obs - ex, obs, ex, ex, opp - ex)
        if pvx is not None:
            exx = K * float(pvx[m].sum())
            acc.add(("G3x", name, lab), ti, obs - exx, obs, exx, exx, opp - exx)
    # ---- G4 -------------------------------------------------------------------------------
    if ts.size:
        W0 = shape.W_of(t0[ee])
        u = -np.expm1(-G * (shape.W_of(ts) - W0)) / -np.expm1(-Gam_e[ee])
        dec = np.minimum((u * 10).astype(int), 9)
        for b in range(10):
            f = float(np.mean(dec == b))
            acc.add(("G4", name, b), ti, f - 0.1, f, 0.1, 0.1 * ts.size, 0.9 * ts.size)
    # ---- descriptive (fixture info, not science) ------------------------------------------
    if extra is not None:
        extra["Hbar_realised"].append(float(Ht.mean()))
        extra["heritable_share_expected"].append(sl.heritable_share(parent, t, n, sp))
        extra["visible_events_per_tape"].append(float(V.sum() / K))
        extra["stem_event_frac"].append(float(V[nn - 1].mean()))


def run_tree_fixture(fi, name, sp, trees, seed, acc, H, replay_n=5):
    extra = {"Hbar_realised": [], "heritable_share_expected": [], "visible_events_per_tape": [],
             "stem_event_frac": []}
    hv(H, "solve_exact", abs(1 - np.exp(-sp.Gam_T) - sp.Hbar) > 1e-12, 1)
    for ti, (parent, branch) in enumerate(trees):
        n = N_TIPS
        t = ed.node_times(parent, branch, n)
        sE = np.random.SeedSequence([seed, 10, fi, ti])
        E = sl.loss_draws(parent.size, K, np.random.default_rng(sE))
        out = sl.silence(parent, t, E, sp)
        score_tree(name, ti, parent, t, n, sp, out, acc, H,
                   np.random.default_rng(np.random.SeedSequence([seed, 11, fi, ti])), extra)
        if ti < replay_n:
            E2 = sl.loss_draws(parent.size, K, np.random.default_rng(sE))
            o2 = sl.silence(parent, t, E2, sp)
            v = (np.sum(o2["L"] != out["L"]) + np.sum(o2["H"] != out["H"]) +
                 np.sum(~((o2["tstar"] == out["tstar"]) |
                          (np.isnan(o2["tstar"]) & np.isnan(out["tstar"])))))
            hv(H, "replay_identical", v, 3 * E.size)
        if name == "CONST":
            prev = None
            for hb in NEST:
                Hh = sl.silence(parent, t, E, sl.SilencingParams(hb, sp.knots, sp.levels))["H"]
                if prev is not None:
                    hv(H, "nested_across_Hbar", np.sum(prev & ~Hh), Hh.size)
                prev = Hh
    hr = np.array(extra["Hbar_realised"])
    return {"Hbar_target": sp.Hbar, "Gam_T": sp.Gam_T,
            "Hbar_realised_mean": float(hr.mean()), "Hbar_realised_sd_across_trees":
            float(hr.std(ddof=1)), "heritable_share_expected_mean":
            float(np.mean(extra["heritable_share_expected"])),
            "visible_events_per_tape_per_tree": float(np.mean(extra["visible_events_per_tape"])),
            "stem_event_frac_mean": float(np.mean(extra["stem_event_frac"]))}


def rung1_panel(fi, fx, seed):
    pdraw = dl.panel_draws(K, np.random.default_rng(np.random.SeedSequence([seed, 1, fi])))
    u, r = dl.panel_tapes(pdraw, fx["panel"])
    beta, _ = dl.panel_beta(pdraw, u, fx["panel"], vd.MASK)
    rec = ed.Recorder(k=K, Lam_T=fx["Lam_T"], r=r, shares=((0.0,),), **fx["rate"])
    return rec, beta


def run_off(fi, trees, seed, H):
    """Hbar = 0 on rung-1 BASE: the loss layer must leave rung 1 untouched, bit for bit."""
    fx = vd.FIXTURES["BASE"]
    rec, beta = rung1_panel(fi, fx, seed)
    sp = sl.SilencingParams(0.0)
    p1, p2 = [], []
    for ti, (parent, branch) in enumerate(trees):
        n = N_TIPS
        D, sym, out, _ = dl.decorate_tips(parent, branch, n, rec,
                                          np.random.default_rng([seed, 2, fi, ti]))
        md = dl.mask_draws(n, K, np.random.default_rng([seed, 3, fi, ti]))
        X, _ = dl.technical_mask(md, beta, vd.MASK.s_alpha)
        Y1 = dl.observe(X, D)
        E = sl.loss_draws(parent.size, K, np.random.default_rng([seed, 10, fi, ti]))
        so = sl.silence(parent, out["t"], E, sp)
        Y2 = sl.observe(X, D, so["H"][:n])
        hv(H, "OFF_H_zero", so["H"].sum() + so["L"].sum(), 2 * so["H"].size)
        hv(H, "OFF_Y_identical", np.sum(Y1 != Y2), Y1.size)
        kept = dl.retained(Y1, md, None)
        p1.append({"D": D, "sym": sym, "Y": Y1, "kept": kept, "name": f"OFF_tree{ti}"})
        p2.append({"D": D, "sym": sym, "Y": Y2, "kept": kept, "name": f"OFF_tree{ti}"})
    dl.export_arm("simOFF_rung1", p1, 0, rec.M, OUTDIR)
    dl.export_arm("simOFF_rung2", p2, 0, rec.M, OUTDIR)
    for kind in ("dropout_matrix", "prefix_codes6"):
        a = np.load(OUTDIR / f"{kind}_simOFF_rung1.npz")
        b = np.load(OUTDIR / f"{kind}_simOFF_rung2.npz")
        v = int(sorted(a.files) != sorted(b.files))
        v += sum(int(not np.array_equal(a[f], b[f])) for f in a.files if f in b.files)
        hv(H, f"OFF_export_identical_{kind}", v, len(a.files))
    return {"trees": len(trees)}


def run_composed(fi, trees, seed, acc, H):
    """rung-1 SLOW + CONST silencing: the three factors compose as independent."""
    fx = vd.FIXTURES["SLOW"]
    rec, beta = rung1_panel(fi, fx, seed)
    sp = SIL["CONST"]
    Hb = sp.Hbar
    mu = rec.Lam_T * rec.r
    e0 = np.exp(-mu)
    name = "COMPOSED"
    parts = []
    for ti, (parent, branch) in enumerate(trees):
        n = N_TIPS
        D, sym, out, _ = dl.decorate_tips(parent, branch, n, rec,
                                          np.random.default_rng([seed, 2, fi, ti]))
        md = dl.mask_draws(n, K, np.random.default_rng([seed, 3, fi, ti]))
        X, pi = dl.technical_mask(md, beta, vd.MASK.s_alpha)
        t = out["t"]
        E = sl.loss_draws(parent.size, K, np.random.default_rng([seed, 10, fi, ti]))
        so = sl.silence(parent, t, E, sp)
        Ht = so["H"][:n]
        Y = sl.observe(X, D, Ht)
        hv(H, "COMPOSED_Y_identity", np.sum(Y != ((~X) & (~Ht) & (D >= 1))), Y.size)
        P = (1 - pi) * (1 - Hb) * (1 - e0)[None, :]
        for z in range(K):
            acc.add(("G5", name, z), ti, np.mean(Y[:, z] - P[:, z]), Y[:, z].mean(),
                    P[:, z].mean(), P[:, z].sum(), (1 - P[:, z]).sum())
        nodes, c1, c2 = vd.split_pairs(parent, n,
                                       np.random.default_rng([seed, 11, fi, ti]))
        Wa = rec.W_of(t[nodes])
        Wz = sl.Shape(sp.knots, sp.levels).W_of(t[nodes])
        y12 = (Y[c1] & Y[c2]).astype(float)
        Ey = ((1 - pi[c1]) * (1 - pi[c2]) * ((1 - Hb) ** (2 - Wz))[:, None] *
              (1 - 2 * e0[None, :] + np.exp(-mu[None, :] * (2 - Wa[:, None]))))
        acc.add(("G5", name, "pair_pooled"), ti, np.mean(y12 - Ey), y12.mean(), Ey.mean(),
                Ey.sum(), (1 - Ey).sum())
        kept = dl.retained(Y, md, None)
        parts.append({"D": D, "sym": sym, "Y": Y, "kept": kept, "H": Ht,
                      "name": f"COMPOSED_tree{ti}"})
    dl.export_arm("simCOMPOSED", parts, 0, rec.M, OUTDIR)
    z = np.load(OUTDIR / "dropout_matrix_simCOMPOSED.npz")
    Yall = np.concatenate([p["Y"] for p in parts])
    Hall = np.concatenate([p["H"] for p in parts])
    hv(H, "COMPOSED_export_recovered", np.sum(z["recovered"] != Yall), Yall.size)
    hv(H, "COMPOSED_export_ABSENT_where_H", np.sum(z["term"][Hall] != dl.ABSENT), int(Hall.sum()))
    return {"mu_range": [float(mu.min()), float(mu.max())], "observed_missing":
            float(1 - Yall.mean()), "silenced_tips": float(Hall.mean())}


def in_verdict(key):
    chk, fx, cell = key[0], key[1], key[2]
    if chk in ("G2n", "G3x"):
        return False
    if chk == "G2" and cell != "pooled":
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1102026)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t_start = time.perf_counter()
    per_cell = 2 if a.quick else 20
    trees = []
    for rho, th in vd.TREE_CELLS:
        tt = ed.load_library_trees(N_TIPS, rho, th)
        assert len(tt) >= per_cell
        trees += tt[:per_cell]
    split = np.random.default_rng(a.seed + 1).permutation(len(trees))
    acc, H, infos = vd.Acc(), {}, {}
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tc = time.perf_counter()
    infos["OFF"] = run_off(0, trees, a.seed, H)
    print(f"OFF: {time.perf_counter() - tc:.0f} s", flush=True)
    for fi, (name, sp) in enumerate(SIL.items(), start=1):
        tc = time.perf_counter()
        infos[name] = run_tree_fixture(fi, name, sp, trees, a.seed, acc, H)
        print(f"{name}: {time.perf_counter() - tc:.0f} s  {infos[name]}", flush=True)
    tc = time.perf_counter()
    infos["COMPOSED"] = run_composed(9, trees, a.seed, acc, H)
    print(f"COMPOSED: {time.perf_counter() - tc:.0f} s  {infos['COMPOSED']}", flush=True)

    rows = acc.summarise(split)
    for r_ in rows:
        r_["in_verdict"] = in_verdict(r_["key"])
    tested = [r_ for r_ in rows if r_["tested"] and r_["in_verdict"]]
    zs = np.array([abs(r_["z"]) for r_ in tested])
    zn = np.array([r_["z_null_split"] for r_ in tested if np.isfinite(r_["z_null_split"])])
    trips = [r_ for r_ in tested if abs(r_["z"]) > vd.Z_PASS]
    h_viol = int(sum(v[0] for v in H.values()))
    summary = {"cells_tested": len(tested),
               "cells_excluded": int(sum(1 for r_ in rows if r_["in_verdict"] and not r_["tested"])),
               "worst_abs_z": float(zs.max()), "cells_over_3": len(trips),
               "trips": [r_["key"] for r_ in trips],
               "expected_trips_by_chance": float(len(tested) * 0.0027),
               "null_z_sd": float(zn.std(ddof=1)), "null_worst_abs_z": float(np.abs(zn).max()),
               "null_cells": int(zn.size), "H_violations_total": h_viol,
               "verdict": "PASS" if (zs.max() <= vd.Z_PASS and h_viol == 0) else
                          ("RERUN (single trip)" if len(trips) == 1 and h_viol == 0 else "FAIL")}
    out = {"settings": {"seed": a.seed, "n": N_TIPS, "k": K, "N": N, "trees": len(trees),
                        "tree_cells_rho_theta": vd.TREE_CELLS,
                        "silencing": {k: json.loads(sl.settings_json(v)) for k, v in SIL.items()},
                        "composed_rung1": "SLOW", "off_rung1": "BASE", "nest_Hbar": NEST,
                        "size_bins": SIZE_BINS, "min_expected": vd.MIN_EXPECTED,
                        "z_pass": vd.Z_PASS},
           "fixture_info": infos, "checks": rows, "H": H, "summary": summary,
           "seconds_total": round(time.perf_counter() - t_start, 1)}
    fname = "validate_silencing_quick.json" if a.quick else "validate_silencing.json"
    (RES / fname).write_text(json.dumps(out, indent=1, default=float))

    print(f"\n{'check':4s} {'fixture':8s} {'cell':6s} {'tested':>6s} {'excl':>5s} {'worst|z|':>9s}")
    groups: dict = {}
    for r_ in rows:
        part = "pooled" if "pooled" in r_["key"][2] else ("bins" if r_["key"][0] == "G2" else "")
        groups.setdefault((r_["key"][0], r_["key"][1], part), []).append(r_)
    for (chk, fx, part), rr in sorted(groups.items()):
        t_ = [x for x in rr if x["tested"]]
        w = max((abs(x["z"]) for x in t_), default=float("nan"))
        note = "" if all(x["in_verdict"] for x in rr) else "  (not in verdict)"
        print(f"{chk:4s} {fx:8s} {part:6s} {len(t_):6d} {len(rr) - len(t_):5d} {w:9.2f}{note}")
    for kk, v in H.items():
        print(f"H    {kk:34s} {v[0]} violations of {v[1]:.3g}")
    s = summary
    print(f"\ncells tested {s['cells_tested']} (excluded {s['cells_excluded']}); worst |z| "
          f"{s['worst_abs_z']:.2f}; cells > 3: {s['cells_over_3']} (chance expects "
          f"{s['expected_trips_by_chance']:.2f}); null z sd {s['null_z_sd']:.2f} over "
          f"{s['null_cells']} cells; H violations {h_viol}  =>  {s['verdict']}")
    print(f"total {out['seconds_total']} s")


if __name__ == "__main__":
    main()
