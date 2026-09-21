#!/usr/bin/env python3
r"""
02_validate_tree.py -- does 01_bdtree.py simulate the process it claims to?

Five checks, each against something ALREADY VERIFIED independently in
notes/sciphy_notes.md S3.3, so this script tests the code and not the algebra.

  A  count law      P(N(T)=n) vs the geometric (1-a)(1-b)b^{n-1}; E[N(T)] vs
                    e^{(b-delta)T}; P(extinct) vs alpha.        [tests the Gillespie loop]
  B  structure      every accepted n-tip tree has exactly n-1 branching times,
                    2n-1 nodes, 2 children per internal node, times inside (0,T).
                    A missed degree-2 suppression shows up here and nowhere else.
                                                                [tests the pruning]
  C  root split     P(smaller root clade = k) vs 2/(n-1), or 1/(n-1) at k = n/2.
                    S3.3.3 verified this is parameter-free, so a skew here means
                    pruning broke exchangeability.              [tests the pruning]
  D  histories      labelled-history frequencies at n=4 (18 of them) and n=5
                    (180) against uniform, with tip labels randomly permuted.
                    Catches a systematically wrong topology that C can miss.
  E  naive vs fast  the two implementations compared in DISTRIBUTION (root split
                    and branching-time ECDF).  They share no code path below
                    reconstruct(), so agreement exercises the two-pass design.

REPORTING.  Deviations are quoted in standard errors, never as p-values -- with
10^5 pooled quantities a KS test rejects almost anything, and this project has
already recorded that significance saturates.  Flag threshold: 3 se.

⚠ This script decides whether the simulator is CORRECT.  It says nothing about
whether forward simulation is AFFORDABLE -- that is 03_cost_curve.py.

usage:  scripts/submit.sh analyses/2026-09_simulator/src/02_validate_tree.py \
            [--reps 20000] [--out results/validate.json]
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("bdtree", _HERE / "01_bdtree.py")
bd = importlib.util.module_from_spec(_spec)
sys.modules["bdtree"] = bd          # @dataclass resolves annotations via sys.modules
_spec.loader.exec_module(bd)

# ⚠⚠ CORRECTION 2026-09-20, after the first pair of runs.  The reading rule
# originally flagged at a FLAT 3.0 se.  Every check reports the MAXIMUM |z| over
# its cells, so the threshold must grow with the number of cells: the null p99 of
# max|z| is 2.56 at m=2, 3.03 at m=4, 3.44 at m=18 and 4.42 at m=180 (N=4000).
# A flat 3.0 therefore called rho8e4 a FAILURE at 3.14 on a 180-cell check whose
# statistic sat at the 62nd percentile of its own null.  Critical values are now
# calibrated per check by multinomial Monte Carlo at the actual (m, N).
# ⚠ The closed form sqrt(2 ln 2m) is NOT an adequate substitute -- it reads 2.04
# where the null mean is 1.38 (m=4) and 3.43 where it is 2.99 (m=180), because
# cell counts at these N are skewed, not normal.  Calibrate, do not approximate.
FLAG_Q = 99.0       # percentile of the null of max|z| at which a check is flagged
_CRIT_CACHE: dict = {}


def crit_value(m, n_draws, q=FLAG_Q, B=40000, seed=12345):
    """Null q-th percentile of max|z| over m equiprobable cells with n_draws draws."""
    key = (int(m), int(n_draws), float(q))
    if key not in _CRIT_CACHE:
        rng = np.random.default_rng(seed)
        p = 1.0 / max(int(m), 2)
        se = np.sqrt(p * (1 - p) / n_draws)
        c = rng.multinomial(n_draws, [p] * max(int(m), 2), size=B)
        _CRIT_CACHE[key] = float(np.percentile((np.abs(c / n_draws - p) / se).max(axis=1), q))
    return _CRIT_CACHE[key]


def se_dev(obs_count, n_total, p_expected):
    """Deviation of an observed frequency from its expectation, in se units."""
    if n_total == 0 or p_expected <= 0 or p_expected >= 1:
        return float("nan")
    se = np.sqrt(p_expected * (1 - p_expected) / n_total)
    return (obs_count / n_total - p_expected) / se


def clade_sizes(tree):
    """Number of tips below every node of a ReconTree.

    ⚠⚠ FIXED 2026-09-20.  The first version iterated node ids DESCENDING, which
    is the right order for the full simulated tree (children get larger ids than
    parents) and exactly the WRONG order for a ReconTree, where slots are handed
    out from the root downward so a child has a SMALLER id than its parent.  A
    parent was therefore accumulated before its children and ended up counting
    only its direct leaf children.  Ordering by node TIME instead is correct
    under either convention: tips sit at T, and a coalescence is always older
    than the coalescences below it, so descending time visits children first.
    """
    k = np.zeros(tree.parent.size, np.int64)
    k[: tree.n_tips] = 1
    for x in np.argsort(-tree.time, kind="stable"):
        p = int(tree.parent[x])
        if p != -1:
            k[p] += k[int(x)]
    return k


def root_split(tree):
    """Size of the SMALLER of the two clades below the root."""
    root = tree.parent.size - 1
    kids = [x for x in range(tree.parent.size) if tree.parent[x] == root]
    k = clade_sizes(tree)
    return int(min(k[kids[0]], k[kids[1]]))


def history_key(tree, rng):
    """Canonical labelled history, with tip labels permuted uniformly at random.

    Tip ids from the simulator come from creation order, which is not a
    meaningful labelling; permuting makes uniformity over labelled histories the
    correct statement of 'every labelled history is equally likely' (S3.3.3).
    """
    lab = rng.permutation(tree.n_tips)
    clade: dict[int, frozenset] = {x: frozenset([int(lab[x])]) for x in range(tree.n_tips)}
    kids: dict[int, list[int]] = {}
    for x in range(tree.parent.size - 1):
        kids.setdefault(int(tree.parent[x]), []).append(x)
    order = sorted(range(tree.n_tips, tree.parent.size), key=lambda x: -tree.time[x])
    key = []
    for x in order:                      # most recent coalescence first
        a, b_ = (clade[c] for c in kids[x])
        key.append(tuple(sorted((tuple(sorted(a)), tuple(sorted(b_))))))
        clade[x] = a | b_
    return tuple(key)


def structural_faults(tree, n_expect, T):
    f = []
    if tree.n_tips != n_expect:
        f.append(f"n_tips {tree.n_tips} != {n_expect}")
    if tree.parent.size != 2 * n_expect - 1:
        f.append(f"{tree.parent.size} nodes != {2 * n_expect - 1}")
    bt = tree.branching_times
    if bt.size != n_expect - 1:
        f.append(f"{bt.size} branching times != {n_expect - 1}")
    if not np.all((bt > 0) & (bt < T)):
        f.append("a branching time falls outside (0, T)")
    counts = np.bincount(tree.parent[tree.parent >= 0], minlength=tree.parent.size)
    internal = counts[n_expect:]
    if not np.all(internal == 2):
        f.append(f"internal node with {sorted(set(internal.tolist()))} children, expected 2")
    if tree.parent[-1] != -1:
        f.append("root has a parent")
    return f


# --------------------------------------------------------------------------


def check_A(b, delta, T, reps, seed):
    rng = np.random.default_rng(seed)
    n_T = np.array([bd.simulate_counts(b, delta, T, rng).n_T for _ in range(reps)])
    alpha, _ = bd.bd_alpha_beta(b, delta, T)
    out = {"reps": reps, "rows": []}

    out["rows"].append({"what": "P(extinct by T)", "obs": float((n_T == 0).mean()),
                        "exp": float(alpha),
                        "se_dev": float(se_dev(int((n_T == 0).sum()), reps, alpha))})
    m_obs, m_exp = float(n_T.mean()), bd.bd_mean(b, delta, T)
    se = n_T.std(ddof=1) / np.sqrt(reps)
    out["rows"].append({"what": "E[N(T)]", "obs": m_obs, "exp": m_exp,
                        "se_dev": float((m_obs - m_exp) / se) if se > 0 else float("nan")})
    for k in (1, 2, 3, 5, 8):
        p = float(bd.bd_pmf(k, b, delta, T))
        out["rows"].append({"what": f"P(N(T)={k})", "obs": float((n_T == k).mean()), "exp": p,
                            "se_dev": float(se_dev(int((n_T == k).sum()), reps, p))})
    out["worst_se"] = float(np.nanmax(np.abs([r["se_dev"] for r in out["rows"]])))
    out["crit99"] = crit_value(len(out["rows"]), reps)
    return out


def collect_fast(n, turnover, rho, reps, seed, tol=0):
    b, delta, T = bd.rate_params(n, rho, turnover)
    trees, stats = [], []
    for i in range(reps):
        tr, st = bd.simulate_tree(b, delta, T, rho, n, seed=seed + i, tol=tol)
        if tr is not None:
            trees.append(tr)
            stats.append(st)
    return trees, stats, (b, delta, T)


def check_BCD(trees, n, T, seed):
    rng = np.random.default_rng(seed)
    faults = []
    for tr in trees:
        faults.extend(structural_faults(tr, n, T))

    ks = np.array([root_split(tr) for tr in trees])
    rows = []
    for k in range(1, n // 2 + 1):
        p = (1.0 if 2 * k == n else 2.0) / (n - 1)
        rows.append({"k": k, "obs": float((ks == k).mean()), "exp": p,
                     "se_dev": float(se_dev(int((ks == k).sum()), ks.size, p))})

    n_exp = math.factorial(n) * math.factorial(n - 1) // 2 ** (n - 1)
    if n_exp > len(trees) / 20:          # too many cells to say anything; C still applies
        hist, devs = {}, []
    else:
        hist = {}
        for tr in trees:
            key = history_key(tr, rng)
            hist[key] = hist.get(key, 0) + 1
        p = 1.0 / n_exp
        devs = [se_dev(hist.get(k, 0), len(trees), p) for k in hist]

    return {
        "structure_faults": faults[:20],
        "n_structure_faults": len(faults),
        "root_split": rows,
        "root_split_worst_se": float(np.nanmax(np.abs([r["se_dev"] for r in rows]))),
        "pruning_load_nodes_per_branchpoint": float(np.mean(
            [(t.n_ancestors - 1) / (n - 1) for t in trees])) if trees else float("nan"),
        "histories_seen": len(hist), "histories_expected": int(n_exp),
        "history_tested": bool(devs),
        "history_worst_se": float(np.nanmax(np.abs(devs))) if devs else float("nan"),
        "history_devs": [round(float(d), 3) for d in devs],
        "history_crit99": crit_value(n_exp, len(trees)) if devs else float("nan"),
        "root_split_crit99": crit_value(max(n // 2, 2), len(trees)),
    }


def check_E(n, turnover, rho, reps, seed, trees_fast, budget_s=600.0):
    """Same conditioning, other implementation: reject naive draws that miss n."""
    b, delta, T = bd.rate_params(n, rho, turnover)
    rng = np.random.default_rng(seed)
    naive, tries, t0 = [], 0, time.perf_counter()
    while len(naive) < reps and tries < reps * 4000:
        if time.perf_counter() - t0 > budget_s:
            break
        tries += 1
        tr, _ = bd.simulate_tree_naive(b, delta, T, rho, rng)
        if tr is not None and tr.n_tips == n:
            naive.append(tr)
    if len(naive) < 50:
        return {"note": f"only {len(naive)} naive trees in {tries} attempts -- comparison skipped"}

    kf = np.array([root_split(t) for t in trees_fast])
    kn = np.array([root_split(t) for t in naive])
    rows = []
    for k in range(1, n // 2 + 1):
        pf, pn = (kf == k).mean(), (kn == k).mean()
        se = np.sqrt(pf * (1 - pf) / kf.size + pn * (1 - pn) / kn.size)
        rows.append({"k": k, "fast": float(pf), "naive": float(pn),
                     "se_dev": float((pf - pn) / se) if se > 0 else 0.0})

    bf = np.sort(np.concatenate([t.branching_times for t in trees_fast]))
    bn = np.sort(np.concatenate([t.branching_times for t in naive]))
    grid = np.linspace(min(bf[0], bn[0]), max(bf[-1], bn[-1]), 200)
    sup = float(np.max(np.abs(np.searchsorted(bf, grid) / bf.size
                              - np.searchsorted(bn, grid) / bn.size)))
    return {"n_naive": len(naive), "naive_attempts": tries,
            "naive_seconds": round(time.perf_counter() - t0, 1), "root_split": rows,
            "root_split_worst_se": float(np.nanmax(np.abs([r["se_dev"] for r in rows]))),
            "branching_time_sup_gap": sup,
            "sup_gap_ref_1.36sqrt": float(1.36 * np.sqrt(1 / len(trees_fast) + 1 / len(naive)))}


def check_F(n, rho, turnover, reps, budget_s, seed):
    """Check B alone, at LARGE n.  No statistics -- structure is asserted per tree.

    ⚠ C and D cannot scale: they need many trees per cell, and the cell count runs
    away (1.6e6 labelled histories already at n=8).  B does scale, and it is the
    check that would catch a slot-assignment, stack or ancestor-walk failure at the
    sizes this simulator is actually FOR.  Park's largest clone is 10,997 cells;
    everything validated so far was n <= 8.
    """
    b, delta, T = bd.rate_params(n, rho, turnover)
    faults, got, t0 = [], [], time.perf_counter()
    for i in range(reps):
        if time.perf_counter() - t0 > budget_s:
            break
        tr, st = bd.simulate_tree(b, delta, T, rho, n, seed=seed + 7919 * i)
        if tr is None:
            continue
        faults.extend(structural_faults(tr, n, T))
        got.append((st.attempts, st.peak_live, tr.n_ancestors, st.seconds))
    if not got:
        return {"n": n, "rho": rho, "turnover": turnover, "trees": 0,
                "note": f"no acceptance inside {budget_s:.0f} s"}
    a = np.array(got, float)
    return {"n": n, "rho": rho, "turnover": turnover, "trees": len(got),
            "n_structure_faults": len(faults), "structure_faults": faults[:10],
            "mean_attempts": float(a[:, 0].mean()), "mean_peak_live": float(a[:, 1].mean()),
            "pruning_load_nodes_per_branchpoint": float((a[:, 2].mean() - 1) / (n - 1)),
            "sec_per_accept": float(a[:, 3].mean())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20000, help="replicates for check A")
    ap.add_argument("--tree-reps", type=int, default=3000, help="accepted trees for B/C/D")
    ap.add_argument("--rho", type=float, default=0.5, help="capture fraction for the tree checks")
    ap.add_argument("--turnover", type=float, nargs="+", default=[0.0, 0.3, 0.75],
                    help="turnovers swept by B/C/D. S3.3.3 says topology is turnover-INDEPENDENT; "
                         "one value tests the theory, not that the simulator honours it.")
    ap.add_argument("--e-turnover", type=float, default=0.3,
                    help="check E runs at this turnover only -- it tests implementation "
                         "agreement, which is not a turnover question")
    ap.add_argument("--skip-bcd", action="store_true")
    ap.add_argument("--f-grid", default="", help='large-n structural points, COMMA separated: "n:rho,n:rho". '
                         'No spaces -- submit.sh passes args through $* inside --wrap, '
                         'which word-splits a quoted multi-token value.')
    ap.add_argument("--f-reps", type=int, default=3)
    ap.add_argument("--f-budget", type=float, default=900.0)
    ap.add_argument("--seed", type=int, default=20260920)
    ap.add_argument("--skip-naive", action="store_true",
                    help="skip check E. The naive oracle is a per-event Python loop, so at "
                         "small rho it needs ~n/rho iterations per attempt and is not worth "
                         "running there; E is informative at rho ~ 0.5.")
    ap.add_argument("--tag", default="rho50")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.out is None:
        a.out = str(_HERE.parent / "results" / f"validate_{a.tag}.json")
    report = {"params": vars(a), "checks": {}}

    print("A  count law", flush=True)
    for turnover in (0.0, 0.3, 0.75):
        b, delta, T = bd.rate_params(64, 0.5, turnover)
        r = check_A(b, delta, T, a.reps, a.seed)
        r["turnover"] = turnover
        report["checks"].setdefault("A_count_law", []).append(r)
        print(f"   turnover {turnover:<5} worst |dev| = {r['worst_se']:.2f} se", flush=True)

    for turn in ([] if a.skip_bcd else a.turnover):
        for n in (4, 5, 8):
            print(f"B/C/D  n = {n}  turnover = {turn}", flush=True)
            trees, stats, (b, delta, T) = collect_fast(n, turn, a.rho, a.tree_reps, a.seed + n)
            r = check_BCD(trees, n, T, a.seed)
            r.update(n=n, turnover=turn, n_trees=len(trees),
                     mean_attempts=float(np.mean([s.attempts for s in stats])) if stats else 0.0)
            report["checks"].setdefault("BCD", []).append(r)
            print(f"   trees {len(trees)}  structure faults {r['n_structure_faults']}  "
                  f"root-split {r['root_split_worst_se']:.2f}/{r['root_split_crit99']:.2f}  "
                  f"histories {r['histories_seen']}/{r['histories_expected']} "
                  f"{r['history_worst_se']:.2f}/{r['history_crit99']:.2f}  "
                  f"prune-load {r['pruning_load_nodes_per_branchpoint']:.1f}", flush=True)

            if n <= 5 and not a.skip_naive and turn == a.e_turnover:
                e = check_E(n, turn, a.rho, min(a.tree_reps, 2000), a.seed + 100 + n, trees,
                            budget_s=300.0)
                e.update(n=n, turnover=turn)
                report["checks"].setdefault("E_naive_vs_fast", []).append(e)
                if "note" in e:
                    print(f"   naive-vs-fast: {e['note']}", flush=True)
                else:
                    print(f"   naive-vs-fast: root-split worst {e['root_split_worst_se']:.2f} se, "
                          f"branching-time sup gap {e['branching_time_sup_gap']:.4f} "
                          f"(ref {e['sup_gap_ref_1.36sqrt']:.4f})", flush=True)

    for point in a.f_grid.replace(",", " ").split():
        fn, frho = point.split(":")
        print(f"F  structure at n = {fn}, rho = {frho}", flush=True)
        f = check_F(int(fn), float(frho), a.e_turnover, a.f_reps, a.f_budget, a.seed)
        report["checks"].setdefault("F_large_n", []).append(f)
        if f["trees"]:
            print(f"   {f['trees']} trees  faults {f['n_structure_faults']}  "
                  f"attempts {f['mean_attempts']:.0f}  peak live {f['mean_peak_live']:.3g}  "
                  f"prune-load {f['pruning_load_nodes_per_branchpoint']:.1f}  "
                  f"{f['sec_per_accept']:.1f} s/tree", flush=True)
        else:
            print(f"   {f['note']}", flush=True)

    ratios = []          # each check's statistic as a FRACTION of its own critical value
    for r in report["checks"].get("A_count_law", []):
        ratios.append(r["worst_se"] / r["crit99"])
    for r in report["checks"].get("BCD", []):
        ratios.append(r["root_split_worst_se"] / r["root_split_crit99"])
        if r.get("history_tested"):
            ratios.append(r["history_worst_se"] / r["history_crit99"])
        if r["n_structure_faults"]:
            ratios.append(1e9)
    for r in report["checks"].get("F_large_n", []):
        if r.get("n_structure_faults"):
            ratios.append(1e9)
    report["worst_ratio_to_crit"] = float(np.nanmax(ratios)) if ratios else 0.0
    report["verdict"] = "PASS" if report["worst_ratio_to_crit"] < 1.0 else "FAIL"

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\n{report['verdict']}  worst statistic sits at "
          f"{report['worst_ratio_to_crit']:.2f} x its own p{FLAG_Q:g} critical value"
          f"\nwrote {out}")


if __name__ == "__main__":
    sys.exit(main())
