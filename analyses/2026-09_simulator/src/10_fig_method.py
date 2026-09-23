#!/usr/bin/env python3
r"""
10_fig_method.py -- FIGURE S1: how the tree layer of the simulator works.

Three standalone panels, each its own file (PDF for the LaTeX write-up, PNG for
slides), plus a JSON of every number the write-up quotes:

  figS1a_forward     one tree through the three stages: the full birth-death
                     process, capture of the survivors, the reconstructed tree
  figS1b_condition   conditioning on clone size: many forward attempts, the
                     Binomial capture draw at T, accept iff K = n.  Right half:
                     the law of K, Monte Carlo against the closed form
  figS1c_library     what the tree library holds: trees stocked per (n, rho, theta)

Write-up: ../writeup/simulator_figures.tex, section S1.

⚠ Panel (a) re-simulates rather than reading the library: the library stores
only RECONSTRUCTED trees, so the full genealogy with its deaths has to be
regenerated.  Capture is the simulator's own rule -- K ~ Binomial(N(T), rho),
accept iff K = n, captured set a uniform K-subset -- and the tree shown is the
first accepted attempt that is also small enough to draw.  That drawability
filter is an illustration choice and is stated in the caption.

⚠ Corrects the first draft of this script (2026-09-22, never committed): its
panel 3 kept full-tree heights, so single-child chains floated unconnected and
the "spliced" nodes were never actually spliced; and it drew exactly n of the
survivors instead of a Binomial(N, rho) capture, so 9 of 22 (41%) were shown
captured under a stated rho of 0.25.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("bdtree", _HERE / "01_bdtree.py")
bd = importlib.util.module_from_spec(_spec)
sys.modules["bdtree"] = bd
_spec.loader.exec_module(bd)

FIG = _HERE.parent / "figures"
RES = _HERE.parent / "results"
LIB = RES / "tree_library"

# reference palette (dataviz skill): slot 1 blue, slot 8 red, inks and a ghost grey
BLUE, RED = "#2a78d6", "#e34948"
INK, INK2, MUTED, GHOST = "#0b0b0b", "#52514e", "#898781", "#d6d5d0"
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# the one parameter set used in panels (a) and (b)
N_TIPS, RHO, THETA = 8, 0.25, 0.5


def save(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / stem}.{{pdf,png}}")


# ---------------------------------------------------------------------------
# closed forms for the capture count K
# ---------------------------------------------------------------------------

def k_law(b, delta, T, rho, kmax):
    r"""P(K = k), k = 0..kmax, for K | N(T) ~ Binomial(N(T), rho).

    With alpha, beta from bd_alpha_beta (P(N=0) = alpha,
    P(N=m) = (1-alpha)(1-beta) beta^{m-1}), thinning keeps the geometric form:
        c = 1 - beta (1 - rho),  beta' = rho beta / c,
        P(K = k) = (1 - alpha) (rho / c) (1 - beta') beta'^{k-1},  k >= 1.
    """
    a, bet = bd.bd_alpha_beta(b, delta, T)
    c = 1.0 - bet * (1.0 - rho)
    bp = rho * bet / c
    k = np.arange(kmax + 1)
    p = np.where(k >= 1, (1 - a) * (rho / c) * (1 - bp) * bp ** np.maximum(k - 1, 0), 0.0)
    p[0] = 1.0 - (1 - a) * rho / c
    return p, dict(alpha=a, beta=bet, c=c, beta_prime=bp)


# ---------------------------------------------------------------------------
# panel (a): one tree, three stages
# ---------------------------------------------------------------------------

def full_tree(b, delta, T, rng):
    """Forward Gillespie keeping death times -- the naive oracle's loop, plus tend."""
    parent, tb, tend, alive = [-1], [0.0], [None], [0]
    t, total, pbirth = 0.0, b + delta, b / (b + delta)
    while alive:
        n = len(alive)
        t += float(rng.exponential(1.0 / (n * total)))
        if t > T:
            break
        i = int(rng.integers(n)); x = alive[i]
        tend[x] = t
        if rng.random() < pbirth:
            for _ in range(2):
                parent.append(x); tb.append(t); tend.append(None)
            alive[i] = len(parent) - 2; alive.append(len(parent) - 1)
        else:
            alive[i] = alive[-1]; alive.pop()
    for x in alive:
        tend[x] = T
    return np.array(parent), np.array(tb), np.array(tend, float), alive


def full_layout(parent, nodes):
    """y per node for a rectangular drawing: leaves in DFS order, internals at the child mean."""
    kids = {}
    for x in nodes:
        if parent[x] != -1:
            kids.setdefault(int(parent[x]), []).append(x)
    y, nxt = {}, [0.0]

    def walk(x):
        c = kids.get(x, [])
        if not c:
            y[x] = nxt[0]; nxt[0] += 1.0
        else:
            for k in c:
                walk(k)
            y[x] = float(np.mean([y[k] for k in c]))
    walk(0)
    return y, kids


def visible_tree(parent, tend, anc, sampled):
    """The reconstructed tree, drawn from the full genealogy.

    Visible nodes = sampled tips + ancestors with two ancestor children.  For each
    visible node v: its time (T for a tip, its division time otherwise), its
    visible parent (-1 for the root, whose branch is the stem from t = 0), and the
    division times of the unary ancestors spliced out of the branch above it.
    """
    akids = {}
    for x in anc:
        if parent[x] != -1:
            akids.setdefault(int(parent[x]), []).append(x)
    sset = set(sampled)
    vis = {x for x in anc if x in sset or len(akids.get(x, [])) == 2}
    vparent, splices = {}, {}
    for v in vis:
        s, p = [], int(parent[v])
        while p != -1 and p not in vis:
            s.append(float(tend[p])); p = int(parent[p])
        vparent[v], splices[v] = p, s
    vkids = {}
    for v in vis:
        if vparent[v] != -1:
            vkids.setdefault(vparent[v], []).append(v)
    return vis, vparent, vkids, splices


def panel_a():
    b, delta, T = bd.rate_params(N_TIPS, RHO, THETA)
    rng = np.random.default_rng(20260922)
    attempts = 0
    while True:
        attempts += 1
        par, tb, tend, alive = full_tree(b, delta, T, rng)
        if not alive:
            continue
        K = int(rng.binomial(len(alive), RHO))              # the simulator's capture rule
        if K != N_TIPS:
            continue
        samp = [int(s) for s in rng.choice(alive, size=K, replace=False)]
        if len(par) <= 120:                                  # drawability filter only
            break
    nodes = list(range(len(par)))
    has_kid = {int(p) for p in par if p != -1}
    dead = {x for x in nodes if x not in has_kid and tend[x] < T - 1e-12}

    anc = set()
    for s in samp:
        x = s
        while x != -1 and x not in anc:
            anc.add(x); x = int(par[x])
    vis, vpar, vkids, splices = visible_tree(par, tend, anc, samp)
    vroot = [v for v in vis if vpar[v] == -1][0]
    vtime = {v: (T if v in samp else float(tend[v])) for v in vis}
    n_spliced = sum(len(s) for s in splices.values())

    # cross-check against the production reconstruction in 01_bdtree.py
    rt = bd.reconstruct(par, tb, np.array(alive), np.array(samp), T)
    assert rt is not None and rt.n_tips == K
    mine = np.sort([vtime[v] for v in vis if v not in samp])
    assert np.allclose(mine, np.sort(rt.time[K:])), "visible_tree disagrees with bd.reconstruct"
    assert rt.n_ancestors == len(anc)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
    Y, KIDS = full_layout(par, nodes)
    sset = set(samp)

    def draw_full(ax, colour, lw):
        for x in nodes:
            ax.plot([tb[x], tend[x]], [Y[x], Y[x]], color=colour(x), lw=lw(x),
                    solid_capstyle="butt", zorder=2 + (x in anc))
            # each half of the connector wears its CHILD's colour, so a captured
            # lineage's ancestry does not run blue into a sister that died
            for k in KIDS.get(x, []):
                ax.plot([tend[x]] * 2, [Y[x], Y[k]], color=colour(k), lw=lw(k),
                        solid_capstyle="butt", zorder=2 + (k in anc))

    draw_full(axes[0], lambda x: RED if x in dead else INK2, lambda x: 1.3)
    for x in dead:
        axes[0].plot(tend[x], Y[x], "x", color=RED, ms=4.5, mew=1.3, zorder=4)
    draw_full(axes[1], lambda x: BLUE if x in anc else (RED if x in dead else GHOST),
              lambda x: 2.4 if x in anc else 1.1)
    for s in samp:
        axes[1].plot(T, Y[s], "o", color=BLUE, ms=5.5, zorder=5)

    # stage 3 -- the reconstructed tree, laid out on ITSELF
    yv, nxt = {}, [0.0]

    def walk(v):
        c = vkids.get(v, [])
        if not c:
            yv[v] = nxt[0]; nxt[0] += 1.0
        else:
            for k in sorted(c, key=lambda z: Y[z]):
                walk(k)
            yv[v] = float(np.mean([yv[k] for k in c]))
    walk(vroot)
    ys = (max(Y.values()) - min(Y.values())) / max(nxt[0] - 1, 1)   # match panel height
    ax = axes[2]
    for v in vis:
        t0 = 0.0 if vpar[v] == -1 else vtime[vpar[v]]
        ax.plot([t0, vtime[v]], [yv[v] * ys] * 2, color=BLUE, lw=2.4, solid_capstyle="butt")
        c = vkids.get(v, [])
        if len(c) == 2:
            ax.plot([vtime[v]] * 2, [yv[c[0]] * ys, yv[c[1]] * ys], color=BLUE, lw=2.4)
        for ts in splices[v]:
            ax.plot(ts, yv[v] * ys, "o", ms=6.5, mfc="white", mec=INK, mew=1.3, zorder=5)
        if v in sset:
            ax.plot(T, yv[v] * ys, "o", color=BLUE, ms=5.5, zorder=5)

    titles = ["1  every division and death",
              f"2  capture each survivor with prob. $\\rho$ = {RHO:g}",
              "3  the reconstructed tree"]
    for a, ti in zip(axes, titles):
        a.set_title(ti, loc="left", fontsize=11)
        a.set_xlabel("time since the founding cell,  $t/T$")
        a.set_xlim(-0.02 * T, 1.05 * T); a.set_yticks([])
        a.set_ylim(min(Y.values()) - 1.2, max(Y.values()) + 0.8)
        for sp in ("top", "right", "left"):
            a.spines[sp].set_visible(False)
        a.axvline(T, color=MUTED, ls=":", lw=1, zorder=1)
    axes[0].legend(handles=[Line2D([], [], color=INK2, lw=1.8, label="lineage"),
                            Line2D([], [], color=RED, lw=1.8, marker="x", ms=5,
                                   label="died before $T$")],
                   frameon=False, fontsize=8.5, loc="upper left")
    axes[1].legend(handles=[Line2D([], [], color=BLUE, lw=2.4, marker="o", ms=5,
                                   label=f"captured ({K}) and its ancestors"),
                            Line2D([], [], color=GHOST, lw=1.6,
                                   label=f"survived, not captured ({len(alive) - K})"),
                            Line2D([], [], color=RED, lw=1.6, label="died")],
                   frameon=False, fontsize=8.5, loc="upper left")
    axes[2].legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=INK, mew=1.3,
                                   ms=6.5, label=f"division spliced out ({n_spliced})")],
                   frameon=False, fontsize=8.5, loc="upper left")
    fig.tight_layout(w_pad=2.5)
    save(fig, "figS1a_forward")

    return dict(b=b, delta=delta, T=T, n=N_TIPS, rho=RHO, theta=THETA,
                attempts_to_first_drawable_accept=attempts,
                full_tree_nodes=len(par), divisions=int((len(par) - 1) // 2),
                deaths=len(dead), survivors=len(alive), captured=K,
                ancestors=len(anc), spliced=n_spliced,
                recon_nodes=2 * K - 1,
                recon_branching_times=[round(float(x), 4) for x in np.sort(rt.time[K:])])


# ---------------------------------------------------------------------------
# panel (b): conditioning on clone size
# ---------------------------------------------------------------------------

def panel_b(n_show=60, n_mc=40000):
    b, delta, T = bd.rate_params(N_TIPS, RHO, THETA)
    kmax = 60
    p_exact, pars = k_law(b, delta, T, RHO, kmax)

    # left: n_show forward attempts, each with its own capture draw
    trajs = []
    for a in range(n_show):
        c_ss, s_ss, _ = np.random.SeedSequence([7, a]).spawn(3)
        tr = bd.simulate_counts(b, delta, T, np.random.default_rng(c_ss), record=True)
        k = int(np.random.default_rng(s_ss).binomial(tr.n_T, RHO)) if tr.n_T else 0
        trajs.append((tr, k))

    # right: the law of K by Monte Carlo, same seeding scheme, counts only
    ks = np.empty(n_mc, int)
    for a in range(n_mc):
        c_ss, s_ss, _ = np.random.SeedSequence([8, a]).spawn(3)
        tr = bd.simulate_counts(b, delta, T, np.random.default_rng(c_ss))
        ks[a] = int(np.random.default_rng(s_ss).binomial(tr.n_T, RHO)) if tr.n_T else 0
    emp = np.bincount(np.minimum(ks, kmax), minlength=kmax + 1) / n_mc
    se = np.sqrt(p_exact * (1 - p_exact) / n_mc)
    z = np.where(se > 0, (emp - p_exact) / np.where(se > 0, se, 1), 0)[: kmax]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 4.4),
                                 gridspec_kw=dict(width_ratios=[1.35, 1]))
    n_acc = n_ext = 0
    for tr, k in trajs:
        t = np.concatenate([[0.0], tr.times, [T]])
        N = np.concatenate([[1], 1 + np.cumsum(tr.steps.astype(int)), [tr.n_T]])
        if tr.n_T == 0:
            # log axis: draw up to the last cell, mark the extinction on the floor
            n_ext += 1
            ax.step(t[:-2], N[:-2], where="post", color=RED, lw=0.8, alpha=0.8, zorder=2)
            ax.plot([t[-3], t[-2]], [N[-3], N[-3]], color=RED, lw=0.8, alpha=0.8, zorder=2)
            ax.plot(t[-2], 0.75, "x", color=RED, ms=4, mew=1.1, zorder=3, clip_on=False)
        elif k == N_TIPS:
            n_acc += 1
            ax.step(t, N, where="post", color=BLUE, lw=2.0, zorder=4)
        else:
            ax.step(t, N, where="post", color=GHOST, lw=0.9, zorder=1)
    ax.axvline(T, color=MUTED, ls=":", lw=1)
    ax.axhline(N_TIPS / RHO, color=MUTED, ls="--", lw=0.9)
    ax.text(0.02 * T, N_TIPS / RHO * 1.08, f"$E[N(T)] = n/\\rho$ = {N_TIPS / RHO:.0f}",
            fontsize=8.5, color=INK2, va="bottom")
    ax.text(1.01 * T, 0.75, "extinct", fontsize=8, color=RED, va="center")
    ax.set_yscale("log")
    ax.set_ylim(0.75, None)
    ax.set_xlabel("time since the founding cell,  $t/T$")
    ax.set_ylabel("cells alive,  $N(t)$  (log scale)")
    ax.set_title(f"{n_show} forward attempts; keep one only if exactly $n$ = {N_TIPS} "
                 "cells are captured", loc="left", fontsize=10.5)
    ax.legend(handles=[Line2D([], [], color=BLUE, lw=2, label=f"accepted: $K = n$ ({n_acc})"),
                       Line2D([], [], color=GHOST, lw=1.5,
                              label=f"rejected: $K \\neq n$ ({n_show - n_acc - n_ext})"),
                       Line2D([], [], color=RED, lw=1, marker="x", ms=4,
                              label=f"extinct before $T$ ({n_ext})")],
              frameon=False, fontsize=8.5, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    kk = np.arange(1, kmax)
    bx.bar(kk, emp[1:kmax], width=0.8, color=GHOST, label=f"simulated, {n_mc:,} attempts")
    bx.bar([N_TIPS], [emp[N_TIPS]], width=0.8, color=BLUE, label=f"accepted value $k = n$")
    bx.plot(kk, p_exact[1:kmax], "o", ms=3.2, color=INK, label="closed form, eq. (7)")
    bx.set_xlabel("cells captured at harvest,  $K$")
    bx.set_ylabel("probability,  $P(K = k)$")
    bx.set_title(f"law of the capture count $K$  ($P(K=0)$ = {p_exact[0]:.3f}, not drawn)",
                 loc="left", fontsize=10.5)
    bx.set_xlim(0, 45)
    bx.legend(frameon=False, fontsize=8.5)
    for sp in ("top", "right"):
        bx.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=3)
    save(fig, "figS1b_condition")

    # the recorded check (analysis_log 2026-09-20): 32.1 attempts at n=8, rho=0.5, theta=0.3
    b2, d2, T2 = bd.rate_params(8, 0.5, 0.3)
    p2, _ = k_law(b2, d2, T2, 0.5, 8)
    a_n = float(p_exact[N_TIPS])
    return dict(**{k: float(v) for k, v in pars.items()},
                P_K0=float(p_exact[0]), P_accept=a_n, attempts_exact=1 / a_n,
                attempts_approx_e_n=float(np.e * N_TIPS / ((1 - pars["alpha"]) * RHO / pars["c"])),
                mc_attempts=n_mc, mc_P_accept=float(emp[N_TIPS]),
                mc_max_abs_z_k1_to_59=float(np.max(np.abs(z[1:]))),
                mc_P_K0=float(emp[0]),
                shown=n_show, shown_accepted=n_acc, shown_extinct=n_ext,
                recorded_check_n8_rho05_th03=float(1 / p2[8]))


# ---------------------------------------------------------------------------
# panel (c): the library inventory
# ---------------------------------------------------------------------------

PAT = re.compile(r"trees_n(\d+)_rho([\d.e-]+)_th([\d.]+)_r(\d+)\.npz")


def panel_c(target=20):
    cells = {}
    for f in LIB.glob("trees_*.npz"):
        m = PAT.fullmatch(f.name)
        if not m:
            continue
        key = (int(m[1]), float(m[2]), float(m[3]))
        with np.load(f) as d:
            cells[key] = cells.get(key, 0) + int(d["seed"].shape[0])
    rhos = sorted({k[1] for k in cells}, reverse=True)
    thetas = sorted({k[2] for k in cells})
    cmap = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#104281"])

    fig, axes = plt.subplots(len(thetas), 1, figsize=(10.5, 1.05 * len(thetas) + 1.0),
                             sharex=True)
    for ax, th in zip(axes, thetas):
        pts = [(n, rhos.index(r), c) for (n, r, t), c in cells.items() if t == th]
        n_, y_, c_ = map(np.array, zip(*pts))
        full = c_ >= target
        ax.scatter(n_[full], y_[full], marker="s", s=16, c=BLUE, linewidths=0)
        ax.scatter(n_[~full], y_[~full], marker="s", s=22, facecolors="white",
                   edgecolors=RED, linewidths=1.1)
        ax.set_yticks(range(len(rhos)))
        ax.set_yticklabels([f"{r:g}" for r in rhos], fontsize=7.5)
        ax.set_ylim(len(rhos) - 0.4, -0.6)
        ax.set_xscale("log")
        ax.text(1.005, 0.5, f"$\\theta$ = {th:g}", transform=ax.transAxes, fontsize=9,
                color=INK2, va="center")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[len(axes) // 2].set_ylabel("capture fraction  $\\rho$")
    axes[-1].set_xlabel("clone size  $n$  (cells, log scale)")
    axes[0].legend(handles=[Line2D([], [], marker="s", ls="", color=BLUE, ms=5,
                                   label=f"{target} trees"),
                            Line2D([], [], marker="s", ls="", mfc="white", mec=RED, ms=5,
                                   label=f"< {target} trees")],
                   frameon=False, fontsize=8, ncol=2, loc="lower left",
                   bbox_to_anchor=(0, 1.0))
    fig.tight_layout(h_pad=0.4)
    save(fig, "figS1c_library")

    under = sorted((k, c) for k, c in cells.items() if c < target)
    return dict(cells=len(cells), trees=int(sum(cells.values())),
                rhos=rhos, thetas=thetas,
                sizes=sorted({k[0] for k in cells}),
                under_filled=[dict(n=k[0], rho=k[1], theta=k[2], trees=c) for k, c in under])


def main():
    out = dict(a=panel_a(), b=panel_b(), c=panel_c())
    RES.mkdir(exist_ok=True)
    p = RES / "figS1_numbers.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("sizes",)}
                      for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
