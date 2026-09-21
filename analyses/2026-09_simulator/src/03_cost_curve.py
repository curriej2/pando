#!/usr/bin/env python3
r"""
03_cost_curve.py -- is pure forward simulation affordable, and if not, where does it break?

THE QUESTION.  If forward-only is affordable through Subclone then S3.3.5's
route (c) and the owed branching-time density verification are dropped and the
BDS density never enters the project.  If it breaks at some clone size n-dagger
we learn which arms need the density route and the equivalence check becomes
load-bearing.  Either answer closes a live item.

THE COST MODEL BEING TESTED.  With n = observed cells in a clone and rho = the
capture fraction (dimensionless), the envelope prediction is

        E(n)  ~  kappa * n^2 / rho        lineage-events per accepted tree

where one lineage-event is one birth or death draw, and kappa is an O(1)
constant absorbing the extinction share and the birth:death ratio.  The two
factors are recorded SEPARATELY -- attempts per acceptance (predicted ~2.7n) and
peak live lineages (predicted ~n/rho) -- because if only their product were
recorded a wrong exponent would be invisible.  NULL / no-effect value:
kappa = 2.7 and exponent gamma = 2.  Wall-clock adds tau = seconds per
lineage-event, which is an implementation constant no theory supplies.

READING RULE.
  * measured E(n) within ~3x of events_pred() AND fitted gamma in [1.7, 2.3]
    => the cost model holds and extrapolating to unmeasured n is safe.
  * gamma > 2.3 => extrapolation is VOID; do not project Subclone from Mouse3.
  * per arm at R replicates: under 1 h => forward-only, drop the density.
    Over 7 d (the partition walltime cap) => density required for that arm.
    In between => a per-arm call.

⚠ Every number here is a COST, not a result about trees.  Correctness is
02_validate_tree.py's job and this script assumes it passed.

⚠ The per-arm projection uses clone sizes WITHOUT the paper's per-cell tape
filter (>=100 tapes Initial/Subclone, >=20 mice; D.4c), which needs the edit
tables.  Unfiltered clones are larger, so the projection is an OVER-estimate --
conservative in the direction that matters for a feasibility call.

usage:  scripts/submit.sh analyses/2026-09_simulator/src/03_cost_curve.py \
            --n 4 32 210 --reps 5 --budget-s 600
"""
from __future__ import annotations

import argparse
import collections
import csv
import importlib.util
import json
import pathlib
import resource
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
_ANA = _HERE.parent
_spec = importlib.util.spec_from_file_location("bdtree", _HERE / "01_bdtree.py")
bd = importlib.util.module_from_spec(_spec)
sys.modules["bdtree"] = bd          # @dataclass resolves annotations via sys.modules
_spec.loader.exec_module(bd)

CLONES = pathlib.Path("/data1/choij10/justin/pando/data/cancer_metastasis/"
                      "clonalbc_percell_hamming1_corrected.csv")
KAPPA_NULL = 2.7          # exact-hit rejection factor: P(k = n) ~ e^{-1}/n => ~2.7n attempts


def events_pred(n, rho, b, delta, T):
    r"""Predicted lineage-events per accepted tree -- the null the measurement is read against.

        attempts       ~  2.7 n / (1 - alpha)          exact-hit rejection, plus extinct replicates
        events/attempt ~  N(T) (1+theta)/(1-theta)     births B = N/(1-theta), deaths = theta B

    ⚠⚠ CORRECTION (2026-09-20, before any run).  The proposal quoted a FLAT null of
    2.7 n^2/rho and a 3x reading band.  That is the theta = 0 case only: the turnover
    factor (1+theta)/(1-theta)^2 is 1.0 / 2.65 / 28.0 at theta = 0 / 0.3 / 0.75, so a
    flat null would have reported a 28x "failure" at high turnover that is not a
    failure at all.  The n-EXPONENT gamma is untouched by this (the factor is
    n-independent), which is why gamma and not the ratio is the load-bearing test.
    """
    theta = delta / b
    alpha, _ = bd.bd_alpha_beta(b, delta, T)
    return (KAPPA_NULL * n / max(1 - alpha, 1e-9)) * (n / rho) * (1 + theta) / (1 - theta)
SEC_PER_HOUR = 3600.0


def peak_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # KB -> GB on Linux


def measure_cell(n, turnover, rho, reps, budget_s, seed):
    """One (n, turnover) grid point.  Returns per-accepted-tree cost measurements."""
    b, delta, T = bd.rate_params(n, rho, turnover)
    got, rejects = [], []
    t0 = time.perf_counter()
    for i in range(reps):
        if time.perf_counter() - t0 > budget_s:
            break
        tr, st = bd.simulate_tree(b, delta, T, rho, n, seed=seed + 1000 * i,
                                  tol=0, keep_rejects=True)
        if tr is None:
            continue
        got.append(st)
        rejects.extend(st.rejects)

    if not got:
        return {"n": n, "turnover": turnover, "accepted": 0,
                "note": f"no acceptance inside {budget_s:.0f} s"}

    att = np.array([s.attempts for s in got], float)
    ev = np.array([s.events_scan + s.events_build for s in got], float)
    sec = np.array([s.seconds for s in got], float)
    live = np.array([s.peak_live for s in got], float)

    # tolerance-window discount, free from the recorded reject tip counts
    k_all = np.array(rejects + [n] * len(got), float)
    p0 = max(float(np.mean(k_all == n)), 1e-12)
    disc = {str(w): float(np.mean(np.abs(k_all - n) <= w) / p0) for w in (0, 1, 2, 5)}

    return {
        "n": n, "turnover": turnover, "b": b, "delta": delta, "T": T, "rho": rho,
        "accepted": len(got),
        "attempts_per_accept": float(att.mean()),
        "attempts_pred": 2.7 * n / max(1 - bd.bd_alpha_beta(b, delta, T)[0], 1e-9),
        "peak_live": float(live.mean()),
        "peak_live_pred_n_over_rho": n / rho,
        "events_per_accept": float(ev.mean()),
        "events_pred": events_pred(n, rho, b, delta, T),
        "events_ratio_obs_over_pred": float(ev.mean() / events_pred(n, rho, b, delta, T)),
        "sec_per_accept": float(sec.mean()),
        "tau_sec_per_event": float(sec.sum() / ev.sum()),
        "peak_rss_gb": float(peak_rss_gb()),
        "tol_speedup_vs_tol0": disc,
    }


def clone_sizes(path=CLONES):
    """(arm, ClonalBC) -> cell count.  Unfiltered; see the header caveat."""
    arms = collections.defaultdict(collections.Counter)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            # ⚠⚠ FIX 2026-09-20.  This stripped only a trailing _<digit>, which left
            # M1_LL, M1_LN, M2_LV ... -- ORGANS, not arms.  Park clones span organs by
            # construction (the Metient migration subset is clones appearing in >=2
            # organs, D.4c), so grouping by Sample SPLITS clones and understates their
            # size.  The arm is the first token.
            arm = row["Sample"].split("_")[0]
            bc = row["ClonalBC"]
            if bc and bc != "None":
                arms[arm][bc] += 1
    return {a: np.array(sorted(c.values(), reverse=True)) for a, c in arms.items()}


def project(cells, sizes, reps, rho):
    """Extrapolate measured cost to whole arms.

    ⚠⚠ FIX 2026-09-20.  The first version fitted a power law to EVENTS and multiplied
    by a MEAN tau (seconds per lineage-event).  tau is not a constant: it ranged
    6.9e-6 s/event at n=4 to 4.4e-8 at n=3,387 in the rho=0.5 pilot, a 155x spread,
    because at small n (or large rho) each accepted tree is many cheap ATTEMPTS and
    per-attempt Python overhead dominates, while at large n the vectorised per-event
    cost does.  A single mean tau is set by the small-n points and then applied at
    n = 27,224, inflating that projection ~36x -- enough to flip a verdict.

    The cost genuinely has two terms, so fit two:

        seconds per accepted tree  =  c1 * attempts  +  c2 * lineage-events

    c1 is per-attempt overhead (s), c2 is the marginal per-event cost (s).  Fitted by
    non-negative least squares across the grid, then evaluated per clone at the
    PREDICTED attempts and events for that n.  The events power law is still reported,
    because gamma is the check on the cost MODEL, but it no longer drives the hours.
    """
    ok = [c for c in cells if c.get("accepted")]
    if len(ok) < 2:
        return {"note": "need >= 2 measured grid points to project"}
    x = np.log(np.array([c["n"] for c in ok], float))
    y = np.log(np.array([c["events_per_accept"] for c in ok], float))
    gamma, inter = np.polyfit(x, y, 1)

    A = np.array([[c["attempts_per_accept"], c["events_per_accept"]] for c in ok], float)
    sec = np.array([c["sec_per_accept"] for c in ok], float)
    (c1, c2), *_ = np.linalg.lstsq(A, sec, rcond=None)
    c1, c2 = max(float(c1), 0.0), max(float(c2), 0.0)
    resid = float(np.max(np.abs(A @ [c1, c2] - sec) / np.maximum(sec, 1e-12)))
    turn = ok[0]["turnover"]

    out = {"fitted_gamma": float(gamma), "gamma_null": 2.0,
           "cost_c1_sec_per_attempt": c1, "cost_c2_sec_per_event": c2,
           "cost_fit_worst_rel_resid": resid,
           "extrapolation_safe": bool(1.7 <= gamma <= 2.3), "arms": {}}
    for arm, s in sorted(sizes.items()):
        s = s[s >= 2]
        b_, d_, T_ = bd.rate_params(1, 1.0, turn), None, None       # placeholder, replaced below
        hrs = 0.0
        for n_c in s:
            bb, dd, TT = bd.rate_params(int(n_c), rho, turn)
            at = 2.7 * n_c / max(1 - bd.bd_alpha_beta(bb, dd, TT)[0], 1e-9)
            ev = events_pred(int(n_c), rho, bb, dd, TT)
            hrs += (c1 * at + c2 * ev)
        hours = float(hrs * reps / SEC_PER_HOUR)
        out["arms"][arm] = {
            "clones": int(s.size), "median_cells": int(np.median(s)), "max_cells": int(s.max()),
            "core_hours_at_R": round(hours, 2), "R": reps,
            "verdict": ("forward-only" if hours < 1 else
                        "density required" if hours > 24 * 7 else "per-arm call"),
        }
    return out


def figure(cells, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for turn in sorted({c["turnover"] for c in cells}):
        pts = sorted([c for c in cells if c["turnover"] == turn and c.get("accepted")],
                     key=lambda c: c["n"])
        if pts:
            ax.plot([c["n"] for c in pts], [c["events_per_accept"] for c in pts],
                    "o-", label=f"turnover {turn:g}")
    if cells:
        rho = cells[0]["rho"]
        for turn in sorted({c["turnover"] for c in cells}):
            pts = sorted([c for c in cells if c["turnover"] == turn], key=lambda c: c["n"])
            ax.plot([c["n"] for c in pts], [c["events_pred"] for c in pts], "k--", lw=0.8,
                    label="prediction" if turn == min(t["turnover"] for t in cells) else None)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("clone size $n$ (cells)")
    ax.set_ylabel("lineage-events per accepted tree")
    ax.set_title(f"forward-simulation cost, " + r"$\rho$" + f" = {cells[0]['rho']:g}" if cells else "")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[4, 32, 210])
    ap.add_argument("--turnover", type=float, nargs="+", default=[0.0, 0.3, 0.75])
    ap.add_argument("--rho", type=float, default=8e-4, help="Park-like capture fraction")
    ap.add_argument("--reps", type=int, default=5, help="accepted trees per grid point")
    ap.add_argument("--project-R", type=int, default=100, help="replicates per clone in the projection")
    ap.add_argument("--budget-s", type=float, default=600.0, help="wall-clock cap per grid point")
    ap.add_argument("--seed", type=int, default=20260920)
    ap.add_argument("--tag", default="pilot")
    a = ap.parse_args()

    cells = []
    for turn in a.turnover:
        for n in a.n:
            c = measure_cell(n, turn, a.rho, a.reps, a.budget_s, a.seed + n)
            cells.append(c)
            if c.get("accepted"):
                print(f"n={n:<6} turn={turn:<5} acc={c['accepted']:<3} "
                      f"att/acc={c['attempts_per_accept']:.3g} (pred {c['attempts_pred']:.3g})  "
                      f"events/acc={c['events_per_accept']:.3g} "
                      f"(x{c['events_ratio_obs_over_pred']:.2f} pred)  "
                      f"{c['sec_per_accept']:.3g} s  RSS {c['peak_rss_gb']:.2f} GB", flush=True)
            else:
                print(f"n={n:<6} turn={turn:<5} {c.get('note')}", flush=True)

    rep = {"params": vars(a), "cells": cells}
    try:
        rep["projection"] = project([c for c in cells if c["turnover"] == 0.3],
                                    clone_sizes(), a.project_R, a.rho)
    except Exception as exc:                                   # clone table absent or renamed
        rep["projection"] = {"note": f"projection skipped: {exc}"}

    res = _ANA / "results" / f"cost_{a.tag}.json"
    res.parent.mkdir(parents=True, exist_ok=True)
    res.write_text(json.dumps(rep, indent=2, default=str))
    fig = figure([c for c in cells if c.get("accepted")], _ANA / "figures" / f"cost_{a.tag}.png")

    p = rep["projection"]
    if "arms" in p:
        print(f"\nfitted gamma = {p['fitted_gamma']:.2f} (null 2.0), "
              f"extrapolation {'SAFE' if p['extrapolation_safe'] else 'VOID'}\n"
              f"cost model: {p['cost_c1_sec_per_attempt']:.3e} s/attempt + "
              f"{p['cost_c2_sec_per_event']:.3e} s/event, "
              f"worst relative residual {p['cost_fit_worst_rel_resid']:.1%}")
        for arm, v in p["arms"].items():
            print(f"  {arm:<12} {v['clones']:>5} clones  median {v['median_cells']:>5}  "
                  f"max {v['max_cells']:>6}  {v['core_hours_at_R']:>12.2f} core-h at R={v['R']}"
                  f"   {v['verdict']}")
    print(f"\nwrote {res}\nwrote {fig}")


if __name__ == "__main__":
    sys.exit(main())
