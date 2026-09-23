#!/usr/bin/env python3
r"""
11_fig_pull_present.py -- FIGURE S3: the pull of the present, and what moves it.

An ILLUSTRATION, not a test (approved 2026-09-22 in that form; the closed-form and
matched-pair tests were proposed and declined).  Nothing here is concluded; every
number printed is descriptive.

What is plotted.  L(t) = lineages in the reconstructed tree at time t, t/T from 0
(the founding cell) to 1 (harvest); the curve is the mean of ln L(t) over the R
trees of one (n, rho, theta) cell, with a band over the middle 50% of trees.
Computed from stored branch lengths on a 400-point grid -- the library's own
50-point LTT is 0.02 T apart, too coarse where the effect lives.

  (a) the phenomenon      rho = 1, theta in {0, 0.5, 0.7}
  (b) capture pushes back theta = 0.5, rho in {1, 0.25, 0.02, 0.002}
  (c) both together       one small panel per rho in {1, 0.25, 0.02, 0.002},
                          theta in {0, 0.3, 0.5, 0.7} as lines

All at n = 63.  rho = 1 is not in the library, so those trees are simulated here
(R = 20 per theta, seconds) with the production simulator.

⚠ ROOT TIME IS NOT STORED.  The library keeps branch lengths with the root's set
to 0, so node times come back relative to the root.  Every tip sits at T, so the
root's absolute time is T minus any root-to-tip path length; recovered that way
and asserted equal across tips.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("bdtree", _HERE / "01_bdtree.py")
bd = importlib.util.module_from_spec(_spec)
sys.modules["bdtree"] = bd
_spec.loader.exec_module(bd)

FIG = _HERE.parent / "figures"
RES = _HERE.parent / "results"
LIB = RES / "tree_library"

N, R = 63, 20
GRID = np.linspace(0.0, 1.0, 401)
TIP_WINDOW = 0.1          # descriptive tip slope: d ln L / dt over the last 10% of T

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
# theta: one sequential blue ramp, light -> dark (magnitude, not identity)
TH_COL = {0.0: "#9ec5f4", 0.3: "#5598e7", 0.5: "#256abf", 0.7: "#0d366b"}
# rho: sequential orange->... no -- rho is ALSO a magnitude; use a second ramp only in (b)
RHO_COL = {1.0: "#0d366b", 0.25: "#256abf", 0.02: "#5598e7", 0.002: "#9ec5f4"}
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "pdf.fonttype": 42,
})


def abs_node_times(parent, branch, n, T=1.0):
    """Absolute times (0 = founding cell) of every node of one stored tree."""
    t = np.zeros(parent.size)
    for x in range(parent.size - 2, -1, -1):        # parents have larger ids
        t[x] = t[int(parent[x])] + branch[x]
    depth = t[:n]
    assert np.allclose(depth, depth[0], atol=1e-4), "tips not ultrametric"
    return t + (T - depth.mean())


def ln_ltt(internal_times):
    """ln L on GRID: one lineage (the stem) plus one per branch point at or before t."""
    return np.log(1 + np.searchsorted(np.sort(internal_times), GRID, side="right"))


def cell_from_library(rho, theta):
    curves = []
    for f in sorted(LIB.glob(f"trees_n{N}_rho{rho:g}_th{theta:g}_r*.npz")):
        with np.load(f) as d:
            for p, bl in zip(d["parent"], d["branch"]):
                t = abs_node_times(p, bl.astype(np.float64), N)
                curves.append(ln_ltt(t[N:]))
    if not curves:
        raise FileNotFoundError(f"no library cell n={N} rho={rho:g} theta={theta:g}")
    return np.array(curves)


def cell_rho1(theta, seed0=20260922):
    b, delta, T = bd.rate_params(N, 1.0, theta)
    curves = []
    for i in range(R):
        tr, _ = bd.simulate_tree(b, delta, T, 1.0, N, seed=seed0 + 1_000_003 * i)
        assert tr is not None and tr.n_tips == N
        curves.append(ln_ltt(tr.time[N:]))
    return np.array(curves)


def get(rho, theta, cache={}):
    k = (rho, theta)
    if k not in cache:
        cache[k] = cell_rho1(theta) if rho == 1.0 else cell_from_library(rho, theta)
    return cache[k]


def tip_slope(c):
    """Mean d ln L/dt over the last TIP_WINDOW of T: visible branching per lineage there."""
    m = c.mean(axis=0)
    i0 = np.searchsorted(GRID, 1 - TIP_WINDOW)
    return float((m[-1] - m[i0]) / (GRID[-1] - GRID[i0]))


def draw(ax, c, colour, label, lw=2.0, band=True, ls="-", z=3):
    m = c.mean(axis=0)
    if band:
        lo, hi = np.percentile(c, [25, 75], axis=0)
        ax.fill_between(GRID, np.exp(lo), np.exp(hi), color=colour, alpha=0.13, lw=0, zorder=1)
    ax.plot(GRID, np.exp(m), color=colour, lw=lw, ls=ls, label=label, zorder=z)


def style(ax, title=None):
    ax.set_yscale("log")
    ax.set_ylim(0.9, 80); ax.set_xlim(0, 1.0)
    ax.set_yticks([1, 2, 5, 10, 20, 63]); ax.set_yticklabels(["1", "2", "5", "10", "20", "63"])
    ax.minorticks_off()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=10.5)


def save(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / stem}.{{pdf,png}}")


def main():
    base = get(1.0, 0.0)
    out = {"n": N, "R": R, "tip_window": TIP_WINDOW, "tip_slope": {}}

    # (a) the phenomenon
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    for th in (0.0, 0.5, 0.7):
        draw(ax, get(1.0, th), TH_COL[th], f"turnover $\\theta$ = {th:g}")
    style(ax, "(a)  complete capture ($\\rho$ = 1): extinction bends the tree up")
    ax.set_xlabel("time since the founding cell,  $t/T$")
    ax.set_ylabel("lineages in the reconstructed tree,  $L(t)$")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    save(fig, "figS3a_pull")

    # (b) capture pushes back
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    draw(ax, base, MUTED, "baseline: $\\theta$ = 0, $\\rho$ = 1", lw=1.3, band=False, ls="--", z=2)
    for rho in (1.0, 0.25, 0.02, 0.002):
        draw(ax, get(rho, 0.5), RHO_COL[rho], f"capture $\\rho$ = {rho:g}")
    style(ax, "(b)  turnover $\\theta$ = 0.5: incomplete capture flattens it")
    ax.set_xlabel("time since the founding cell,  $t/T$")
    ax.set_ylabel("lineages in the reconstructed tree,  $L(t)$")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    save(fig, "figS3b_capture")

    # (c) both together
    rhos = (1.0, 0.25, 0.02, 0.002)
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.9), sharey=True)
    for ax, rho in zip(axes, rhos):
        draw(ax, base, MUTED, "baseline", lw=1.1, band=False, ls="--", z=2)
        for th in (0.0, 0.3, 0.5, 0.7):
            c = get(rho, th)
            draw(ax, c, TH_COL[th], f"$\\theta$ = {th:g}", band=False)
            out["tip_slope"][f"rho={rho:g},theta={th:g}"] = tip_slope(c)
        style(ax, f"capture $\\rho$ = {rho:g}")
        ax.set_xlabel("$t/T$")
    axes[0].set_ylabel("lineages,  $L(t)$")
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.tight_layout(w_pad=1.2)
    save(fig, "figS3c_grid")

    for th in (0.0, 0.5, 0.7):
        out["tip_slope"].setdefault(f"rho=1,theta={th:g}", tip_slope(get(1.0, th)))
    for rho in rhos:
        b, d, T = bd.rate_params(N, rho, 0.5)
        out.setdefault("rates_theta0.5", {})[f"rho={rho:g}"] = dict(b=b, delta=d, r=b - d)
    (RES / "figS3_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
