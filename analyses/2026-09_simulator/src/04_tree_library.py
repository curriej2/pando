#!/usr/bin/env python3
r"""
04_tree_library.py -- build the library of simulated trees.

WHAT THIS IS FOR.  Trees are the expensive layer and every later sweep over
editing and dropout parameters re-decorates the SAME trees, which makes those
comparisons paired and removes tree-to-tree variance from the contrast.  The
library is an INPUT; nothing is concluded here.

⭐ SIZE GRID, NOT PARK'S CLONE LIST (Justin, 2026-09-21).  We do not need one
tree per Park clone -- Initial alone has 2,544 of them, nearly all tiny.  We
need trees SPANNING the observed range of clone sizes.  An arm is then
reassembled later by drawing, for each real clone of size n, the library tree
at the nearest stocked size.
⚠ That substitution is only sound if the grid is dense enough in log space.
This grid steps by ~1.78x and pins the five observed arm maxima exactly
(133 Initial, 210 M3, 1,607 M1, 3,387 M2, 11,081 Subclone), so no real clone is
ever more than ~33% from a stocked size.  ⚠ The statistics this library feeds
(the ladder, the variogram, unanimity) are STRONGLY clone-size dependent -- the
project already records that best-k inverts with clone size and that k must be
quoted relative to n_C -- so reassembling the arm's true size DISTRIBUTION by
resampling is not optional, it is the point of stocking the range.

THETA GRID.  theta = delta/b in {0.3,...,0.7}, centred on the only
edit-independent handle in the data: the pool went ~8,000 -> ~160,000 in 10 days
(paper p.4) = a realised net doubling of 55.5 h, against NCI-H1299's intrinsic
22-30 h, giving theta = 1 - r/b = 0.46-0.60.  ⚠ Confounded -- the deficit could
be slower division rather than more death -- so theta is gridded, not fixed.

RHO GRID, per regime, NOT global.
  0.25   Initial/Pre-TX: 37,810 edit-table cells from a stated ~160,000 pool.
  0.10   Subclone lower bound: 8 colonies from single founders, 295-11,081 cells
         each, against ~35,000 expected after 35 d at the realised doubling.
  0.01 } the mice: UNKNOWN.  Growth was tracked by IVIS bioluminescence, which
  0.002}  is relative flux, not a cell count, and no tumour mass or dissociation
          yield is reported.  Gridded until Park supplies a number.
⚠⚠ rho = 8e-4 -- used throughout this analysis before 2026-09-21 -- is SciPhy's
fixed value for HEK293T (notes lines 899, 1742), NOT a Park measurement.  The
mouse rhos above bracket it; Initial and Subclone are 2-3 orders of magnitude off it.
The mouse rhos are capped at the mouse size range (max clone 3,387), since
simulating an 11,081-cell clone at mouse rho describes no arm that exists.

OUTPUT.  One .npz per (n, rho, theta, chunk): rectangular (R, 2n-1) arrays,
`parent` int32 and `branch` float32.  ⚠ Branch LENGTHS, not node times: times
are O(T) and float32 differencing costs ~9% on the shortest measured branch
(1.4e-6 T), while a branch length itself is exact to 3.6e-8 in float32.  Plus
seed provenance (seed, accepting attempt index, per-tree hash) so any tree can
be replayed bit-identically, and edit-free summaries (tree length, B1 balance,
LTT curve) computed at write time so nothing re-walks a tree later.
⚠ Prefix-clade-size-by-depth is NOT computed here -- prefix clades are defined by
shared EDIT prefixes, so that statistic needs the editing layer.

usage:  04_tree_library.py --task K --ntasks N [--reps 20] [--outdir ...]
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("bdtree", _HERE / "01_bdtree.py")
bd = importlib.util.module_from_spec(_spec)
sys.modules["bdtree"] = bd
_spec.loader.exec_module(bd)

# ⚠⚠ The cost model was fitted at turnover 0.3 ONLY (03_cost_curve.py's project()
# filters on it), then applied here up to turnover 0.7.  Measured on this run, the
# actual/predicted ratio rises monotonically with turnover -- median 0.79 at 0.3,
# 1.27 at 0.7, tail to 4.12x -- because total nodes ~ 2 N(T)/(1-theta) is 2.9x N(T)
# at theta=0.3 but 6.7x at theta=0.7, and those extra nodes go through the PYTHON
# genealogy loop (~1 us/event) rather than the vectorised pass-1 (~38 ns/event).
# Aggregate was fine (396 core-h actual vs 339 predicted, 1.17x); the TAIL is what
# hit the 12 h walltime on 4 of 120 tasks.  ⇒ size high-turnover work generously.
C1, C2 = 2.945e-4, 3.835e-8          # fitted cost model, 03_cost_curve.py, at turnover 0.3
ARM_MAX = [133, 210, 1607, 3387, 11081]
THETAS = [0.3, 0.4, 0.5, 0.6, 0.7]
RHOS = {0.25: 11081, 0.10: 11081, 0.01: 3387, 0.002: 3387}    # rho -> max size stocked
LTT_POINTS = 50
CHUNK_SEC = 5400.0                   # split a cell's replicates so no chunk exceeds ~1.5 h


def size_grid():
    g = np.unique(np.round(np.geomspace(2, 11081, 16)).astype(int)).tolist()
    return sorted(set(g + ARM_MAX))


def pred_sec(n, rho, theta):
    b, d, T = bd.rate_params(int(n), rho, theta)
    at = 2.7 * n / max(1 - bd.bd_alpha_beta(b, d, T)[0], 1e-9)
    ev = at * (n / rho) * (1 + theta) / (1 - theta)
    return C1 * at + C2 * ev


def build_jobs(reps, sizes=None, only=None):
    """only = explicit [(n, rho, theta), ...], bypassing the grid -- used to backfill
    cells lost when a task hit the walltime (see the cost-model note below)."""
    if only:
        out = []
        for n, rho, th in only:
            sec = pred_sec(n, rho, th)
            per = max(1, min(reps, int(CHUNK_SEC // max(sec, 1e-9))))
            for r0 in range(0, reps, per):
                k = min(per, reps - r0)
                out.append((int(n), float(rho), float(th), r0, k, sec * k))
        out.sort(key=lambda j: -j[5])
        return out
    return _build_grid_jobs(reps, sizes)


def _build_grid_jobs(reps, sizes=None):
    """(n, rho, theta, rep0, nrep) chunks, each costing <= CHUNK_SEC when possible."""
    jobs = []
    for rho, nmax in RHOS.items():
        for n in (sizes if sizes is not None else size_grid()):
            if n > nmax:
                continue
            for th in THETAS:
                s = pred_sec(n, rho, th)
                per = max(1, min(reps, int(CHUNK_SEC // max(s, 1e-9))))
                for r0 in range(0, reps, per):
                    k = min(per, reps - r0)
                    jobs.append((int(n), float(rho), float(th), r0, k, s * k))
    jobs.sort(key=lambda j: -j[5])               # heaviest first, for round-robin balance
    return jobs


def b1_balance(parent, n_tips):
    """SciPhy's B1: sum over non-root internal nodes of 1/(max edges to a leaf below)."""
    n_nodes = parent.size
    height = np.zeros(n_nodes, np.int64)
    kids = [[] for _ in range(n_nodes)]
    for x in range(n_nodes - 1):
        kids[int(parent[x])].append(x)
    for x in range(n_nodes):                     # ascending: children before parents
        if kids[x]:
            height[x] = 1 + max(height[c] for c in kids[x])
    return float(sum(1.0 / height[x] for x in range(n_tips, n_nodes - 1) if height[x] > 0))


def node_times(parent, branch):
    """Absolute node times from stored branch lengths, root at 0.

    In a ReconTree a parent always has a LARGER id than its child, so descending
    ids visits parents first and one pass suffices.  Provided for readers of the
    library; nothing here calls it.
    """
    t = np.zeros(parent.size, np.float64)
    for x in range(parent.size - 2, -1, -1):
        t[x] = t[int(parent[x])] + branch[x]
    return t


def run_cell(n, rho, theta, rep0, nrep, seed0, outdir):
    b, delta, T = bd.rate_params(n, rho, theta)
    n_nodes = 2 * n - 1
    parent = np.zeros((nrep, n_nodes), np.int32)
    branch = np.zeros((nrep, n_nodes), np.float32)
    seeds = np.zeros(nrep, np.int64); att = np.zeros(nrep, np.int64)
    hashes = []; tlen = np.zeros(nrep); b1 = np.zeros(nrep)
    ltt = np.zeros((nrep, LTT_POINTS), np.int32)
    tgrid = np.linspace(0, T, LTT_POINTS)
    faults = []
    t0 = time.perf_counter()

    for i in range(nrep):
        s = seed0 + 1_000_003 * (rep0 + i)
        tr, st = bd.simulate_tree(b, delta, T, rho, n, seed=s)
        if tr is None:
            faults.append(f"rep {rep0+i}: no acceptance"); continue
        # check B, on every stored tree -- free, and the only check that scales
        bt = tr.branching_times
        if tr.parent.size != n_nodes or bt.size != n - 1 or not np.all((bt > 0) & (bt < T)):
            faults.append(f"rep {rep0+i}: structure")
        cnt = np.bincount(tr.parent[tr.parent >= 0], minlength=n_nodes)
        if not np.all(cnt[n:] == 2):
            faults.append(f"rep {rep0+i}: internal degree")

        bl = np.zeros(n_nodes, np.float64)
        nz = tr.parent >= 0
        bl[nz] = tr.time[nz] - tr.time[tr.parent[nz]]
        parent[i] = tr.parent.astype(np.int32)
        branch[i] = bl.astype(np.float32)
        seeds[i] = s; att[i] = st.attempts - 1
        hashes.append(hashlib.sha1(tr.parent.tobytes() + tr.time.tobytes()).hexdigest()[:16])
        tlen[i] = bl.sum(); b1[i] = b1_balance(tr.parent, n)
        ltt[i] = 1 + np.searchsorted(np.sort(tr.time[n:]), tgrid, side="right")

    out = outdir / f"trees_n{n}_rho{rho:g}_th{theta:g}_r{rep0}.npz"   # n is unique per file
    np.savez_compressed(
        out, parent=parent, branch=branch, seed=seeds, attempt=att,
        tree_length=tlen, b1=b1, ltt=ltt, ltt_time=tgrid,
        meta=json.dumps({"n": n, "rho": rho, "theta": theta, "b": b, "delta": delta, "T": T,
                         "rep0": rep0, "nrep": nrep, "hashes": hashes,
                         "faults": faults, "seconds": time.perf_counter() - t0,
                         "code": "01_bdtree.py two-pass; replay = SeedSequence([seed, attempt])"}))
    return out, faults, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--ntasks", type=int, required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260921)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--sizes", default=None,
                    help="comma-separated sizes REPLACING the default grid. Used for the "
                         "dense small-n pass: 97.1%% of Park clones are <=64 cells, and "
                         "stocking every integer there costs 1.6 core-h, making the "
                         "nearest-size substitution EXACT rather than <=34%% off.")
    ap.add_argument("--thetas", default=None, help="comma list overriding THETAS")
    ap.add_argument("--rhos", default=None,
                    help='comma list of "rho:maxsize" overriding RHOS, e.g. "0.05:11081,0.02:3387"')
    ap.add_argument("--only", default=None,
                    help='explicit cells "n,rho,theta;n,rho,theta;..." bypassing the grid')
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    outdir = pathlib.Path(a.outdir or (_HERE.parent / "results" / "tree_library"))
    outdir.mkdir(parents=True, exist_ok=True)
    if a.thetas:
        THETAS[:] = [float(x) for x in a.thetas.split(",")]
    if a.rhos:
        RHOS.clear()
        RHOS.update({float(c.split(":")[0]): int(c.split(":")[1]) for c in a.rhos.split(",")})
    sizes = [int(x) for x in a.sizes.split(",")] if a.sizes else None
    only = ([tuple(float(v) for v in c.split(",")) for c in a.only.split(";")]
            if a.only else None)
    jobs = build_jobs(a.reps, sizes, only)
    mine = jobs[a.task::a.ntasks]

    if a.dry_run:
        tot = sum(j[5] for j in jobs)
        loads = [sum(j[5] for j in jobs[t::a.ntasks]) / 3600 for t in range(a.ntasks)]
        print(f"{len(jobs)} chunks, {sum(j[4] for j in jobs):,} trees, "
              f"{tot/3600:.1f} core-h total")
        print(f"per-task predicted hours: min {min(loads):.2f}  median "
              f"{np.median(loads):.2f}  max {max(loads):.2f}")
        return

    print(f"task {a.task}/{a.ntasks}: {len(mine)} chunks, "
          f"{sum(j[5] for j in mine)/3600:.2f} h predicted", flush=True)
    allf = []
    for n, rho, th, r0, k, s in mine:
        out, faults, sec = run_cell(n, rho, th, r0, k, a.seed, outdir)
        allf += faults
        print(f"  n={n:<6} rho={rho:<6g} th={th:<4g} reps {r0}-{r0+k-1}  "
              f"{sec:8.1f} s (pred {s:8.1f})  faults {len(faults)}  -> {out.name}", flush=True)
    print(f"task {a.task} done, {len(allf)} structure faults total")
    if allf:
        print("FAULTS:", allf[:20])


if __name__ == "__main__":
    sys.exit(main())
