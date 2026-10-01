#!/usr/bin/env python3
r"""
18_validate_dropout.py -- validate the rung-1 dropout layer (17_dropout.py) against closed forms.
Approved 2026-09-30 ("PROPOSAL (2026-09-29)", ../CLAUDE.md, as amended).

A failed check is a BUG, not a finding.  Every entry is a deviation from a closed form in standard
errors, the se taken ACROSS TREES (cells of one tree share ancestry), except F0, whose unit is a
panel.  Yardstick calibrated by a half-split of the trees, as in 15_validate_editing.py.

FIXTURES (test cases, not science; N = 6, k = 30, n = 210, 200 library trees, Lam_pre = 0,
lineage channel only).  One panel per fixture, shared by all 200 trees (as in an experiment).
  BASE    Lam_T 5.5, constant rate, phi 0, s_r 0;         pbar 0.3, s_beta 1, s_alpha 0.7
  CIS     BASE + phi 0.15, r_cl 0.6, s_r 0.3, Delta_beta 2
  SLOW    CIS's panel at Lam_T 2.0 with Session 18's pause-and-burst rate -- censoring is common,
          so F4o has something to detect, and W(t_a) != t_a is exercised
  FILTER  BASE + ClonalBC retention psi = (0: .3, 12: .5, 20: .8, 25: .95); QC R_min = 12 at export

CHECKS
  F0   the panel law over 20,000 panels (unit = panel): E[u] = phi; E[r] = 1;
       E[r u] = phi r_cl rbar_0; E[beta - mu_beta] = Delta_beta phi;
       E[(beta - mu_beta) ln r] = Delta_beta phi (ln rbar_0 + ln r_cl - s_r^2/2)
  F1   technical missing vs sigma(alpha_c + beta_z), 10 bins of that probability
  F2   observed missing per tape vs 1 - (1 - pi_cz)(1 - e^{-mu_z}),  mu_z = r_z Lam_T
  F3   depth given observed, by locus class, vs P(min(Pois(mu_z), N) = d) / (1 - e^{-mu_z})
  F4t  one tip from each side of every internal node, per tape: X1 X2 vs pi1 pi2 -- the mask
       must not read the tree.  10 bins of W(t_a)
  F4o  the same pairs: O1 O2 vs 1 - (1-pi1)(1-e) - (1-pi2)(1-e)
                                 + (1-pi1)(1-pi2)(1 - 2e + e^{-mu (2 - W(t_a))}),  e = e^{-mu}
  F4n  (DIAGNOSTIC, outside the verdict) O1 O2 vs the independence prediction P(O1) P(O2) --
       shows F4o has teeth where censoring is common
  F5   (FILTER) ClonalBC retention vs psi(R_c), per step
  H    hard assertions: the mu_beta solve hits pbar; Y == ~X & (D >= 1); Y = 0 wherever D = 0;
       decoration and mask replay bit-identical; the exported caches -- term / depth / recovered /
       clone / R_c, prefix rows = (R_c >= R_min) & kept, codes -1 exactly off the observed prefix,
       and code <-> prefix a bijection at every depth

⚠ SCORING AMENDED 2026-09-30 (Justin), after the seed-30092026 FAIL: see in_verdict().  F4 is
scored on one pooled cell per fixture; F4o counts only in SLOW and CIS; F4d (both relatives
unedited vs e^{-mu (2 - W)}) is added as a diagnostic.  The generator and closed forms are unchanged.

READING RULE (fixed before the run): PASS = worst |z| <= 3.0 over every tested cell (F4n excluded)
AND zero H violations.  A cell is tested if its expected positive and negative counts, summed over
trees, are each >= 100 and it has >= 20 trees; excluded cells are reported.  One trip -> re-run at
a new seed; two, or a sign pattern within a check, is a defect.

usage:  18_validate_dropout.py [--seed S] [--quick]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
import time

import numpy as np
from scipy.stats import poisson

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("dropout", _HERE / "17_dropout.py")
dl = importlib.util.module_from_spec(_spec)
sys.modules["dropout"] = dl
_spec.loader.exec_module(dl)
ed = dl.ed

RES = _HERE.parent / "results"
N_TIPS, K, N = 210, 30, 6
TREE_CELLS = [(rho, th) for rho in (0.0005, 0.002, 0.02, 0.1, 0.25) for th in (0.3, 0.7)]
MIN_EXPECTED, Z_PASS, N_PANELS = 100.0, 3.0, 20_000
TV = dict(rate_knots=(0.0, 0.25, 0.35, 0.6, 1.0), rate_levels=(1.2, 0.0, 1.6, 0.8))
CONST = dict(rate_knots=(0.0, 1.0), rate_levels=(1.0,))
PSI = ((0, 0.3), (12, 0.5), (20, 0.8), (25, 0.95))
R_MIN_FILTER = 12
MASK = dl.MaskParams(pbar=0.3, s_beta=1.0, s_alpha=0.7)
P_CIS = dl.PanelParams(phi=0.15, r_cl=0.6, s_r=0.3, delta_beta=2.0)
FIXTURES = {
    "BASE":   dict(Lam_T=5.5, rate=CONST, panel=dl.PanelParams(), psi=None),
    "CIS":    dict(Lam_T=5.5, rate=CONST, panel=P_CIS, psi=None),
    "SLOW":   dict(Lam_T=2.0, rate=TV, panel=P_CIS, psi=None),
    "FILTER": dict(Lam_T=5.5, rate=CONST, panel=dl.PanelParams(), psi=PSI),
}


class Acc:
    """Per-cell, per-unit record: (residual, observed, expected, expected+, expected-)."""

    def __init__(self):
        self.cells: dict = {}

    def add(self, key, unit, resid, obs, exp, epos, eneg):
        self.cells.setdefault(key, {})[unit] = (float(resid), float(obs), float(exp),
                                                float(epos), float(eneg))

    def summarise(self, split_perm=None, min_units=20):
        rows = []
        for key, per in self.cells.items():
            units = sorted(per)
            v = np.array([per[u] for u in units])
            epos, eneg = v[:, 3].sum(), v[:, 4].sum()
            tested = epos >= MIN_EXPECTED and eneg >= MIN_EXPECTED and len(units) >= min_units
            r = v[:, 0]
            sd = r.std(ddof=1) if r.size > 1 else 0.0
            z = float(r.mean() / (sd / np.sqrt(r.size))) if tested and sd > 0 else \
                (0.0 if tested else np.nan)
            znull = np.nan
            if tested and split_perm is not None:
                pos = {t: i for i, t in enumerate(split_perm)}
                h = np.array([pos.get(t, -1) % 2 for t in units])
                a, b = r[h == 0], r[h == 1]
                if a.size > 5 and b.size > 5:
                    znull = float((a.mean() - b.mean()) /
                                  np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size))
            rows.append({"key": list(map(str, key)), "n_units": len(units), "tested": bool(tested),
                         "obs": float(v[:, 1].mean()), "exp": float(v[:, 2].mean()), "z": z,
                         "z_null_split": znull, "expected_pos": float(epos),
                         "expected_neg": float(eneg)})
        return rows


def trunc_pmf(mu, N):
    p = poisson.pmf(np.arange(N), np.asarray(mu, float)[..., None])
    return np.concatenate([p, 1.0 - p.sum(-1, keepdims=True)], axis=-1)


def split_pairs(parent, n, rng):
    """For every internal node: its id and one random tip from each child's subtree."""
    nn = parent.size
    anc = np.zeros((nn, nn), bool)
    for x in range(nn - 1, -1, -1):
        if parent[x] >= 0:
            anc[x] = anc[parent[x]]
        anc[x, x] = True
    kids = [[] for _ in range(nn)]
    for x in range(nn - 1):
        kids[parent[x]].append(x)
    nodes, c1, c2 = [], [], []
    for x in range(n, nn):
        assert len(kids[x]) == 2, f"node {x} has {len(kids[x])} children"
        a = np.flatnonzero(anc[:n, kids[x][0]])
        b = np.flatnonzero(anc[:n, kids[x][1]])
        nodes.append(x); c1.append(rng.choice(a)); c2.append(rng.choice(b))
    return np.array(nodes), np.array(c1), np.array(c2)


def check_F0(seed, acc):
    """The panel law, unit = panel."""
    for tag, pp in (("CIS", P_CIS), ("STRONG", dl.PanelParams(0.3, 0.4, 0.3, 2.0))):
        rbar0 = 1 / (1 - pp.phi + pp.phi * pp.r_cl)
        exp_ = {"E[u]": pp.phi, "E[r]": 1.0, "E[r u]": pp.phi * pp.r_cl * rbar0,
                "E[beta-mu]": pp.delta_beta * pp.phi,
                "E[(beta-mu) ln r]": pp.delta_beta * pp.phi *
                (np.log(rbar0) + np.log(pp.r_cl) - pp.s_r ** 2 / 2)}
        ss = np.random.SeedSequence([seed, 50, len(tag)]).spawn(N_PANELS)
        for i, s in enumerate(ss):
            dr = dl.panel_draws(K, np.random.default_rng(s))
            u, r = dl.panel_tapes(dr, pp)
            b = dl.beta_offset(dr, u, pp, MASK)
            obs = {"E[u]": u.mean(), "E[r]": r.mean(), "E[r u]": (r * u).mean(),
                   "E[beta-mu]": b.mean(), "E[(beta-mu) ln r]": (b * np.log(r)).mean()}
            for key, o in obs.items():
                acc.add(("F0", tag, key), i, o - exp_[key], o, exp_[key], 1e9, 1e9)


def run_fixture(ci, name, fx, trees, seed, acc, H, exports, replay_n=5):
    pdraw = dl.panel_draws(K, np.random.default_rng(np.random.SeedSequence([seed, 1, ci])))
    u, r = dl.panel_tapes(pdraw, fx["panel"])
    beta, mu_beta = dl.panel_beta(pdraw, u, fx["panel"], MASK)
    base = beta - mu_beta
    hit = abs(dl.mean_missing(mu_beta, base, MASK.s_alpha) - MASK.pbar)
    H.setdefault("mu_beta_solve", [0, 0])
    H["mu_beta_solve"][0] += int(hit > 1e-9); H["mu_beta_solve"][1] += 1
    rec = ed.Recorder(k=K, Lam_T=fx["Lam_T"], r=r, shares=((0.0,),), **fx["rate"])
    mu = rec.Lam_T * rec.r                                       # (k,) total expected edits
    e0 = np.exp(-mu)
    q3 = trunc_pmf(mu, N)[:, 1:] / (1 - e0)[:, None]              # (k, N): P(D=d | D>=1)
    info = {"closed_tapes": int(u.sum()), "mu_beta": float(mu_beta),
            "mu_range": [float(mu.min()), float(mu.max())],
            "P_unedited_range": [float(e0.min()), float(e0.max())]}
    parts, secs = [], []
    for ti, (parent, branch) in enumerate(trees):
        n = N_TIPS
        sdec = np.random.SeedSequence([seed, 2, ci, ti])
        smask = np.random.SeedSequence([seed, 3, ci, ti])
        spair = np.random.SeedSequence([seed, 4, ci, ti])
        t_s = time.perf_counter()
        D, sym, out, pd = dl.decorate_tips(parent, branch, n, rec, np.random.default_rng(sdec))
        md = dl.mask_draws(n, K, np.random.default_rng(smask))
        X, pi = dl.technical_mask(md, beta, MASK.s_alpha)
        Y = dl.observe(X, D)
        O = ~Y
        secs.append(time.perf_counter() - t_s)
        # H: observed view
        for key, viol, cmp_ in (("Y_identity", int(np.sum(Y != ((~X) & (D >= 1)))), Y.size),
                                ("Y0_where_D0", int(np.sum(Y & (D == 0))), int((D == 0).sum()))):
            a = H.setdefault(key, [0, 0]); a[0] += viol; a[1] += cmp_
        if ti < replay_n:
            D2, sym2, _, _ = dl.decorate_tips(parent, branch, n, rec, np.random.default_rng(sdec))
            X2, _ = dl.technical_mask(dl.mask_draws(n, K, np.random.default_rng(smask)), beta,
                                      MASK.s_alpha)
            v = int(np.sum(D2 != D) + np.sum(sym2 != sym) + np.sum(X2 != X))
            a = H.setdefault("replay_identical", [0, 0]); a[0] += v; a[1] += D.size * (2 + N)
        # F1
        b1 = np.minimum((pi * 10).astype(int), 9)
        for b in range(10):
            m = b1 == b
            if m.any():
                acc.add(("F1", name, b), ti, np.mean(X[m] - pi[m]), X[m].mean(), pi[m].mean(),
                        pi[m].sum(), (1 - pi[m]).sum())
        # F2
        P2 = 1 - (1 - pi) * (1 - e0)[None, :]
        for z in range(K):
            acc.add(("F2", name, z), ti, np.mean(O[:, z] - P2[:, z]), O[:, z].mean(),
                    P2[:, z].mean(), P2[:, z].sum(), (1 - P2[:, z]).sum())
        # F3
        for cls in (0, 1):
            zc = u == cls
            m = Y & zc[None, :]
            if not m.any():
                continue
            cnt = int(m.sum())
            zz = np.nonzero(m)[1]
            Dm = D[m]
            for d in range(1, N + 1):
                qe = q3[zz, d - 1]
                acc.add(("F3", name, f"class{cls}_d{d}"), ti, np.mean((Dm == d) - qe),
                        (Dm == d).mean(), qe.mean(), qe.sum(), (1 - qe).sum())
        # F4t / F4o / F4n
        nodes, c1, c2 = split_pairs(out["parent"], n, np.random.default_rng(spair))
        Wa = rec.W_of(out["t"][nodes])                             # (m,)
        wb = np.minimum((Wa * 10).astype(int), 9)
        p1, p2 = pi[c1], pi[c2]                                    # (m, k)
        x12 = (X[c1] & X[c2]).astype(float)
        o12 = (O[c1] & O[c2]).astype(float)
        ee = e0[None, :]
        both_edited = 1 - 2 * ee + np.exp(-mu[None, :] * (2 - Wa[:, None]))
        Eo = 1 - (1 - p1) * (1 - ee) - (1 - p2) * (1 - ee) + (1 - p1) * (1 - p2) * both_edited
        En = P2[c1] * P2[c2]
        Et = p1 * p2
        d12 = ((D[c1] == 0) & (D[c2] == 0)).astype(float)
        Ed = np.broadcast_to(np.exp(-mu[None, :] * (2 - Wa[:, None])), d12.shape)
        # ⚠ SCORING CHANGE 2 (Justin, 2026-09-30): each F4 check is scored on ONE pooled cell per
        # fixture (all nodes x tapes of a tree); the 10 W-bins are secondary.  Bins share their
        # trees, so a sign pattern across them is not independent evidence (README Session 20).
        for tag, obs, ex in (("F4t", x12, Et), ("F4o", o12, Eo), ("F4n", o12, En),
                             ("F4d", d12, Ed)):
            acc.add((tag, name, "pooled"), ti, np.mean(obs - ex), obs.mean(), ex.mean(),
                    ex.sum(), (1 - ex).sum())
        for b in range(10):
            m = wb == b
            if not m.any():
                continue
            for tag, obs, ex in (("F4t", x12, Et), ("F4o", o12, Eo), ("F4n", o12, En)):
                acc.add((tag, name, b), ti, np.mean(obs[m] - ex[m]), obs[m].mean(), ex[m].mean(),
                        ex[m].sum(), (1 - ex[m]).sum())
        # F5 + filters
        kept = dl.retained(Y, md, fx["psi"])
        if fx["psi"] is not None:
            R = Y.sum(1)
            pv = dl.psi_of(R, fx["psi"])
            lo = [s[0] for s in fx["psi"]] + [K + 1]
            for j in range(len(fx["psi"])):
                m = (R >= lo[j]) & (R < lo[j + 1])
                if m.any():
                    acc.add(("F5", name, f"R{lo[j]}-{lo[j + 1] - 1}"), ti,
                            np.mean(kept[m] - pv[m]), kept[m].mean(), pv[m].mean(),
                            pv[m].sum(), (1 - pv[m]).sum())
        parts.append({"D": D, "sym": sym, "Y": Y, "kept": kept, "name": f"{name}_tree{ti}"})
    if name in exports:
        R_min = R_MIN_FILTER if fx["psi"] is not None else 0
        H.update(check_export(name, parts, R_min, rec.M))
    info["seconds_per_tree"] = float(np.mean(secs))
    return info


def check_export(name, parts, R_min, M):
    """Write the caches, read them back, and check every encoding rule."""
    outdir = RES / "dropout_validation"
    dl.export_arm(f"sim{name}", parts, R_min, M, outdir)
    z = np.load(outdir / f"dropout_matrix_sim{name}.npz")
    codes = np.load(outdir / f"prefix_codes6_sim{name}.npz")["codes"]
    D = np.concatenate([p["D"] for p in parts]).astype(int)
    Y = np.concatenate([p["Y"] for p in parts])
    S = np.concatenate([p["sym"] for p in parts]).astype(np.int64)
    kept = np.concatenate([p["kept"] for p in parts])
    cid = np.concatenate([np.full(p["D"].shape[0], i) for i, p in enumerate(parts)])
    h = {}
    h[f"{name}_recovered"] = (int(np.sum(z["recovered"] != Y)), Y.size)
    term_exp = np.where(~Y, dl.ABSENT, np.where(D >= N, dl.COMPLETE, dl.UNEDITED))
    h[f"{name}_term"] = (int(np.sum(z["term"] != term_exp)), Y.size)
    h[f"{name}_depth"] = (int(np.sum(z["depth"] != np.where(Y, D, 0))), Y.size)
    h[f"{name}_Rc"] = (int(np.sum(z["recovered"].sum(1) != Y.sum(1))), Y.shape[0])
    h[f"{name}_clone"] = (int(np.sum(z["clone"] != np.where(kept, cid, -1))), Y.shape[0])
    rows = (Y.sum(1) >= R_min) & kept
    h[f"{name}_prefix_rows"] = (int(codes.shape[0] != rows.sum()), 1)
    Dr, Yr, Sr = D[rows], Y[rows], S[rows]
    viol_det, viol_bij, cmp_det, cmp_bij = 0, 0, 0, 0
    for d in range(N):
        det = Yr & (Dr > d)
        viol_det += int(np.sum((codes[..., d] >= 0) != det)); cmp_det += det.size
        c = codes[..., d][det]
        pre = Sr[..., :d + 1][det]
        if c.size:
            _, inv = np.unique(pre, axis=0, return_inverse=True)
            nc, npre = np.unique(c).size, inv.max() + 1
            npair = np.unique(c.astype(np.int64) * (npre + 1) + inv.ravel()).size
            viol_bij += int(not (nc == npre == npair)); cmp_bij += 1
    h[f"{name}_codes_determined"] = (viol_det, cmp_det)
    h[f"{name}_codes_bijection"] = (viol_bij, cmp_bij)
    return {k: list(v) for k, v in h.items()}


CENSORING_FIXTURES = ("SLOW", "CIS")


def in_verdict(key):
    """Scoring as amended by Justin 2026-09-30, after the seed-30092026 FAIL (README Session 20):
      - F4 checks count through their POOLED cell only; W-bins are secondary;
      - F4o counts only where censoring is not negligible (SLOW, CIS) -- elsewhere it duplicates
        F4t (z correlation 0.98-1.00);
      - F4n (independence prediction) and F4d (both unedited) are diagnostics."""
    chk, fx, cell = key[0], key[1], key[2]
    if chk in ("F4n", "F4d"):
        return False
    if chk in ("F4t", "F4o") and cell != "pooled":
        return False
    if chk == "F4o" and fx not in CENSORING_FIXTURES:
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=29092026)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t_start = time.perf_counter()
    global N_PANELS
    per_cell = 2 if a.quick else 20
    if a.quick:
        N_PANELS = 500
    trees = []
    for rho, th in TREE_CELLS:
        tt = ed.load_library_trees(N_TIPS, rho, th)
        assert len(tt) >= per_cell
        trees += tt[:per_cell]
    split = np.random.default_rng(a.seed + 1).permutation(len(trees))
    acc0, acc, H = Acc(), Acc(), {}
    check_F0(a.seed, acc0)
    infos = {}
    for ci, (name, fx) in enumerate(FIXTURES.items()):
        tc = time.perf_counter()
        infos[name] = run_fixture(ci, name, fx, trees, a.seed, acc, H,
                                  exports=("BASE", "FILTER"))
        print(f"{name}: {len(trees)} trees, {time.perf_counter() - tc:.0f} s  {infos[name]}",
              flush=True)
    rows = acc0.summarise(None, min_units=20) + acc.summarise(split)
    for r_ in rows:
        r_["in_verdict"] = in_verdict(r_["key"])
    tested = [r_ for r_ in rows if r_["tested"] and r_["in_verdict"]]
    zs = np.array([abs(r_["z"]) for r_ in tested])
    zn = np.array([r_["z_null_split"] for r_ in tested if np.isfinite(r_["z_null_split"])])
    trips = [r_ for r_ in tested if abs(r_["z"]) > Z_PASS]
    h_viol = int(sum(v[0] for v in H.values()))
    summary = {"cells_tested": len(tested),
               "cells_excluded": int(sum(1 for r_ in rows if r_["in_verdict"] and not r_["tested"])),
               "worst_abs_z": float(zs.max()), "cells_over_3": len(trips),
               "expected_trips_by_chance": float(len(tested) * 0.0027),
               "null_z_sd": float(zn.std(ddof=1)), "null_worst_abs_z": float(np.abs(zn).max()),
               "null_cells": int(zn.size), "H_violations_total": h_viol,
               "verdict": "PASS" if (zs.max() <= Z_PASS and h_viol == 0) else
                          ("RERUN (single trip)" if len(trips) == 1 and h_viol == 0 else "FAIL")}
    out = {"settings": {"seed": a.seed, "n": N_TIPS, "k": K, "N": N, "trees": len(trees),
                        "tree_cells_rho_theta": TREE_CELLS, "n_panels_F0": N_PANELS,
                        "mask": dl.asdict(MASK), "psi": PSI, "R_min_filter": R_MIN_FILTER,
                        "fixtures": {k: {"Lam_T": v["Lam_T"], "rate": v["rate"],
                                         "panel": dl.asdict(v["panel"]), "psi": v["psi"]}
                                     for k, v in FIXTURES.items()},
                        "min_expected": MIN_EXPECTED, "z_pass": Z_PASS},
           "fixture_info": infos, "checks": rows, "H": H, "summary": summary,
           "seconds_total": round(time.perf_counter() - t_start, 1)}
    fname = "validate_dropout_quick.json" if a.quick else "validate_dropout.json"
    (RES / fname).write_text(json.dumps(out, indent=1, default=float))

    print(f"\n{'check':4s} {'fixture':7s} {'cell':6s} {'tested':>6s} {'excl':>5s} {'worst|z|':>9s}")
    groups: dict = {}
    for r_ in rows:
        part = "pooled" if r_["key"][2] == "pooled" else ("bins" if r_["key"][0].startswith("F4")
                                                          else "")
        groups.setdefault((r_["key"][0], r_["key"][1], part), []).append(r_)
    for (chk, fx, part), rr in sorted(groups.items()):
        t_ = [x for x in rr if x["tested"]]
        w = max((abs(x["z"]) for x in t_), default=float("nan"))
        note = "" if all(x["in_verdict"] for x in rr) else "  (not in verdict)"
        print(f"{chk:4s} {fx:7s} {part:6s} {len(t_):6d} {len(rr) - len(t_):5d} {w:9.2f}{note}")
    for kk, v in H.items():
        print(f"H    {kk:28s} {v[0]} violations of {v[1]:.3g}")
    s = summary
    print(f"\ncells tested {s['cells_tested']} (excluded {s['cells_excluded']}); worst |z| "
          f"{s['worst_abs_z']:.2f}; cells > 3: {s['cells_over_3']} (chance expects "
          f"{s['expected_trips_by_chance']:.2f}); null z sd {s['null_z_sd']:.2f} over "
          f"{s['null_cells']} cells; H violations {h_viol}  =>  {s['verdict']}")
    print(f"total {out['seconds_total']} s")


if __name__ == "__main__":
    main()
