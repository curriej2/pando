#!/usr/bin/env python3
r"""
16_fig_editing.py -- Fig S5: the editing layer, illustrated.  Three standalone panels.

Justin chose a, b, c on 2026-09-28; d (the homoplasy cost of a signal channel) is DEFERRED to a
fuller homoplasy test.  These are ILLUSTRATIONS, not tests -- the validation record is
15_validate_editing.py / README "Session 18".  Where a closed form exists it is drawn as a line
over the simulated points (notes §S4.9.5), so each panel shows the layer working in plain terms.

  S5a  one 8-cell library tree, 3 tapes per cell drawn as rows of 6 slots; every edit shaded by
       how many cells share it.  Cells inherit their ancestors' edits as a common prefix.
  S5b  shared slots between two cells against the time their lineages split, one comparison per
       (internal node, tape) -- one tip drawn from each side -- so no comparison is counted twice
       (README "Session 18": random cell pairs mostly re-compare the same edits).  Line: the
       FULL expectation -- the inherited depth E[min(Poisson(Lam_T t), N)] plus chance matches
       past the split, sum_m q^m P(both lineages write >= m more)^2.  Points sit at the mean split
       time within each bin, not the bin centre (splits bunch late in a bin).  ⚠ Both corrected
       after the first render, where drawing only the inherited depth at bin centres put the
       points above the line in 19 of 20 bins (by ~0.03 slots early, ~0.01 late).
       No signal channel in a and b (p = 0): pure lineage recording, q = q_L.
  S5c  a signal pulse (share 0.02 -> 0.30 -> 0.02, on from t = 0.4 to 0.7) at a CONSTANT editing
       rate, read from tapes that filled all 6 slots: top, when each slot was written (K_{j,6}
       in calendar time) under the shaded on-window; bottom, the signal share per slot, simulated
       points (per tip, as the closed form is stated -- NOT deduplicated) against the prediction.
       ⚠ Constant rate, unlike the validation's TV set-up with its zero-rate window, so the slot
       percentages differ from the 3.2 / ... / 15.6% quoted there.

  S5d  (added 2026-09-28 at Justin's request) the editing rate programmed over time, with the
       validation's pause-and-burst schedule (levels 1.2 / 0 / 1.6 / 0.8 on knots 0, .25, .35, .6,
       1, normalised; x Lam_T): top, the realised rate -- edits in a time bin / tape-time spent
       OPEN (not yet full) in it, summed over every lineage alive -- against the programmed step,
       plus the same edits / ALL tape-time, which sags as tapes fill; bottom, filled slots per tape
       along lineages against E[min(Poisson(Lam(t)), N)], with the constant-rate curve faint.
       Rate bins are 0.025 wide so every knot is a bin edge; the realised rate is a ratio of sums
       over trees, se with trees as clusters.  ⚠ Uses the simulator's TRUE edit times: it shows
       the rate can be IMPOSED, not that it can be RECOVERED from sequenced tapes (§S4.7).
       (The deferred homoplasy curve moves to the homoplasy figure rather than keep the letter d.)

Standard errors are across TREES (200 at n = 210, the validation's spread over rho and theta).
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.special import gammainc, gammaln
from scipy.stats import beta as beta_dist
from scipy.stats import poisson

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("editing", _HERE / "14_editing.py")
ed = importlib.util.module_from_spec(_spec)
sys.modules["editing"] = ed
_spec.loader.exec_module(ed)

FIG = _HERE.parent / "figures"
RES = _HERE.parent / "results"
SEED = 24092026
N_TIPS = 210
TREE_CELLS = [(rho, th) for rho in (0.0005, 0.002, 0.02, 0.1, 0.25) for th in (0.3, 0.7)]
CARTOON = dict(n=8, rho=0.25, theta=0.5, rep=0, k=3)
PULSE = dict(comp_knots=(0.0, 0.4, 0.7, 1.0), shares=((0.02,), (0.30,), (0.02,)))
RATE = dict(rate_knots=(0.0, 0.25, 0.35, 0.6, 1.0), rate_levels=(1.2, 0.0, 1.6, 0.8))
RATE_BIN = 0.025
N_BINS_B = 20

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID_C, BAND = "#e6e5e1", "#ecebe7"
LINE = "#0d366b"
# ordered quantities -> one sequential blue ramp, light -> dark (never lighter than step 250)
RAMP6 = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
SHARE_LEVELS = [(1, 1, "1 cell"), (2, 2, "2"), (3, 4, "3–4"), (5, 7, "5–7"), (8, 8, "all 8")]
SHARE_COL = ["#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK, "pdf.fonttype": 42,
})


def style(ax, grid=True):
    ax.minorticks_off()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if grid:
        ax.grid(axis="y", color=GRID_C, lw=0.6, zorder=0)


def save(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / stem}.{{pdf,png}}", flush=True)


def children_of(parent):
    kids = [[] for _ in range(parent.size)]
    for x in range(parent.size - 1):
        kids[parent[x]].append(x)
    return kids


def tips_below(parent, n):
    """Boolean (nn, n): tip c lies below node v."""
    nn = parent.size
    below = np.zeros((nn, n), bool)
    below[np.arange(n), np.arange(n)] = True
    for x in range(nn - 1):                       # children before parents (smaller ids)
        below[parent[x]] |= below[x]
    return below


def kernel_density_t(j, d, N, Lam_T, t):
    """Density, in calendar time t (constant rate), of the j-th edit's time on a tape that ends at
    depth d at t = 1 -- K_{j,d} of §S4.9.5 with x = Lam_T t."""
    x, mu = Lam_T * t, Lam_T
    if d < N:
        dens = beta_dist.pdf(x / mu, j, d + 1 - j)
    else:
        dens = np.exp((j - 1) * np.log(np.maximum(x, 1e-300)) - x - gammaln(j))
        if N - j > 0:
            dens = dens * gammainc(N - j, np.maximum(mu - x, 0.0))
    return dens / np.trapezoid(dens, t)


def expected_shared(tA, Lam_T, N, q):
    """E[shared prefix of two cells whose lineages split at tA] (constant rate, constant xi):
    the depth D_A ~ min(Poisson(Lam_T tA), N) is shared exactly; past it each lineage writes
    min(Poisson(Lam_T (1 - tA)), N - D_A) more, independently, and slot D_A + m matches by chance
    only if every earlier one did, so P(extra >= m) = q^m P(both write >= m)^2."""
    tA = np.atleast_1d(np.asarray(tA, float))
    out = np.zeros(tA.size)
    for i, t in enumerate(tA):
        pA = poisson.pmf(np.arange(N), Lam_T * t)
        pA = np.append(pA, 1 - pA.sum())
        mu2 = Lam_T * (1 - t)
        tot = 0.0
        for d in range(N + 1):
            extra = 0.0
            for m in range(1, N - d + 1):
                p_ge = 1 - poisson.cdf(m - 1, mu2)          # P(min(Pois, N-d) >= m), m <= N-d
                extra += q ** m * p_ge ** 2
            tot += pA[d] * (d + extra)
        out[i] = tot
    return out


def panel_a(numbers):
    c = CARTOON
    rec = ed.Recorder(k=c["k"], shares=((0.0,),))
    parent, branch = ed.load_library_trees(c["n"], c["rho"], c["theta"])[c["rep"]]
    out = ed.decorate(parent, branch, c["n"], rec, np.random.default_rng(SEED))
    pd, tsym, ttime, torig = ed.assemble(out)
    n, nn, N, k = c["n"], parent.size, rec.N, rec.k
    t, t0 = out["t"], out["t0"]
    below = tips_below(parent, n)
    share = below.sum(1)                                           # cells under each branch
    kids = children_of(parent)
    # tip order from a depth-first walk, so the tree draws without crossings
    order, stack = [], [nn - 1]
    while stack:
        x = stack.pop()
        if x < n:
            order.append(x)
        stack.extend(sorted(kids[x], key=lambda y: -t[y]))
    y = np.zeros(nn)
    y[order] = np.arange(n)[::-1]
    for x in range(n, nn):                                         # parents after children
        y[x] = np.mean([y[ch] for ch in kids[x]])

    def level(m):
        return next(i for i, (lo, hi, _) in enumerate(SHARE_LEVELS) if lo <= m <= hi)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for x in range(nn):
        ax.plot([t0[x], t[x]], [y[x], y[x]], color=MUTED, lw=1.4, zorder=2, solid_capstyle="butt")
        if kids[x]:
            ys = [y[ch] for ch in kids[x]]
            ax.plot([t[x], t[x]], [min(ys), max(ys)], color=MUTED, lw=1.4, zorder=2)
    for x in range(nn):                                            # the edits, on their branch
        m = out["bcount"][x] > 0
        for z in np.nonzero(m)[0]:
            for jj in range(out["bcount"][x, z]):
                ax.scatter(out["btime"][x, z, jj], y[x], s=26, zorder=4,
                           color=SHARE_COL[level(share[x])], edgecolors="white", linewidths=0.8)
    # the tapes, right of the tips
    x0, w, h, gap = 1.06, 0.07, 0.24, 0.28
    for tip in range(n):
        for z in range(k):
            yy = y[tip] + (1 - z) * gap - h / 2
            for j in range(N):
                o = torig[tip, z, j]
                fc = "white" if o < 0 else SHARE_COL[level(share[o])]
                ax.add_patch(Rectangle((x0 + j * w, yy), w * 0.9, h, facecolor=fc,
                                       edgecolor=GRID_C if o < 0 else "white", lw=0.8, zorder=3))
    for j in range(N):
        ax.text(x0 + j * w + w * 0.45, n - 0.45, str(j + 1), ha="center", va="bottom",
                color=MUTED, fontsize=8)
    ax.text(x0 + N * w / 2, n - 0.05, "slot", ha="center", va="bottom", color=INK2, fontsize=8.5)
    ax.set_xlim(-0.02, x0 + N * w + 0.02)
    ax.set_ylim(-0.7, n + 0.25)
    ax.set_yticks([])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("time (fraction of the experiment)")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=8, color=col, label=lab)
               for col, (_, _, lab) in zip(SHARE_COL, SHARE_LEVELS)]
    handles.append(plt.Line2D([], [], marker="s", ls="", ms=8, markerfacecolor="white",
                              markeredgecolor=MUTED, label="empty slot"))
    ax.legend(handles=handles, title="edit shared by", frameon=False, fontsize=8,
              title_fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 1.02))
    ax.set_title("a   Cells inherit their ancestors' edits as a shared prefix",
                 loc="left", fontsize=10.5)
    fig.text(0.0, -0.17, f"one library tree ({n} cells, capture {c['rho']}, turnover "
             f"{c['theta']}), {k} tapes per cell shown; dots = edits on the branch that wrote "
             f"them, boxes = each cell's tapes at harvest\nexpected edits per tape over the "
             f"experiment $\\Lambda_T$ = {rec.Lam_T}; tapes empty at the clone founder",
             color=MUTED, fontsize=8, ha="left", va="top", transform=ax.transAxes)
    save(fig, "figS5a_inheritance")
    numbers["a"] = {"tree": c, "seed": SEED, "edits_on_tree": int(out["bcount"].sum()),
                    "tip_depths": pd[:n].tolist()}


def panel_b(trees, numbers):
    rec = ed.Recorder(k=30, shares=((0.0,),))
    N, k = rec.N, rec.k
    q = float(rec.inner(0.5, 0.5))
    edges = np.linspace(0, 1, N_BINS_B + 1)
    per_tree = np.full((len(trees), N_BINS_B), np.nan)
    per_tree_t = np.full((len(trees), N_BINS_B), np.nan)
    excess_sum, excess_n = 0.0, 0
    ss = np.random.SeedSequence([SEED, 5]).spawn(len(trees))
    for ti, (parent, branch) in enumerate(trees):
        g_dec, g_pick = [np.random.default_rng(s) for s in ss[ti].spawn(2)]
        out = ed.decorate(parent, branch, N_TIPS, rec, g_dec)
        pd, tsym, _, _ = ed.assemble(out)
        below = tips_below(parent, N_TIPS)
        kids = children_of(parent)
        nodes = np.arange(N_TIPS, parent.size)
        a_tip = np.array([g_pick.choice(np.nonzero(below[kids[v][0]])[0]) for v in nodes])
        b_tip = np.array([g_pick.choice(np.nonzero(below[kids[v][1]])[0]) for v in nodes])
        s1, s2 = tsym[a_tip], tsym[b_tip]                          # (m, k, N)
        eq = (s1 == s2) & (s1 >= 0)
        cp = np.cumprod(eq, axis=-1).sum(-1)                       # shared prefix, (m, k)
        excess_sum += float((cp - pd[nodes]).sum()); excess_n += cp.size
        tb = np.clip(np.searchsorted(edges, out["t"][nodes], side="right") - 1, 0, N_BINS_B - 1)
        for b in range(N_BINS_B):
            m = tb == b
            if m.any():
                per_tree[ti, b] = cp[m].mean()
                per_tree_t[ti, b] = out["t"][nodes][tb == b].mean()
    ok = np.sum(np.isfinite(per_tree), 0) >= 20
    mid = 0.5 * (edges[:-1] + edges[1:])
    xpos = np.nanmean(per_tree_t, 0)
    mean = np.nanmean(per_tree, 0)
    se = np.nanstd(per_tree, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(per_tree), 0))
    tg = np.linspace(0, 1, 401)
    cf = expected_shared(tg, rec.Lam_T, N, q)
    cf_inh = expected_shared(tg, rec.Lam_T, N, 0.0)

    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.axhline(N, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
    ax.text(0.01, N + 0.08, f"tape full ({N} slots)", color=MUTED, fontsize=8.5)
    ax.plot(tg, cf, color=LINE, lw=2, zorder=3, label="prediction")
    ax.errorbar(xpos[ok], mean[ok], yerr=se[ok], fmt="o", ms=5, color=RAMP6[2], ecolor=RAMP6[2],
                elinewidth=1, capsize=0, mec="white", mew=0.8, zorder=4,
                label="simulated (mean ± se across trees)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, N + 0.5)
    ax.set_xlabel("when the two lineages split (fraction of the experiment)")
    ax.set_ylabel("slots the two cells share, per tape")
    ax.legend(frameon=False, loc="lower right", fontsize=8.5)
    style(ax)
    ax.set_title("b   The later two lineages split, the more slots their cells share",
                 loc="left", fontsize=10.5)
    fig.text(0.0, -0.17, f"200 library trees of {N_TIPS} cells, 30 tapes, lineage symbols only; "
             f"one comparison per (split, tape), one cell from each side\nline = slots filled "
             f"before the split, plus chance matches after it (at most "
             f"{np.max(cf - cf_inh):.3f} slots); $\\Lambda_T$ = {rec.Lam_T} edits per tape",
             color=MUTED, fontsize=8, ha="left", va="top", transform=ax.transAxes)
    save(fig, "figS5b_shared_prefix")
    pred_x = expected_shared(xpos, rec.Lam_T, N, q)
    numbers["b"] = {"bin_mid": mid.round(3).tolist(), "mean_split_time": np.round(xpos, 4).tolist(),
                    "mean": np.round(mean, 4).tolist(),
                    "se": np.round(se, 4).tolist(), "bin_tested": ok.tolist(),
                    "prediction_at_mean_split_time": np.round(pred_x, 4).tolist(),
                    "z": np.round((mean - pred_x) / se, 2).tolist(), "q": q,
                    "max_chance_contribution_slots": float(np.max(cf - cf_inh)),
                    "chance_excess_slots": excess_sum / excess_n, "comparisons": excess_n}


def panel_c(trees, numbers):
    rec = ed.Recorder(k=30, **PULSE)
    N, A = rec.N, rec.A
    per_tree = np.full((len(trees), N), np.nan)
    full_frac = []
    ss = np.random.SeedSequence([SEED, 6]).spawn(len(trees))
    for ti, (parent, branch) in enumerate(trees):
        out = ed.decorate(parent, branch, N_TIPS, rec, np.random.default_rng(ss[ti]))
        pd, tsym, _, _ = ed.assemble(out)
        td, ts = pd[:N_TIPS], tsym[:N_TIPS]
        full = td == N
        full_frac.append(full.mean())
        sig = (ts >= 0) & (rec.channel_of[np.maximum(ts, 0)] < A)
        if full.any():
            per_tree[ti] = sig[full].mean(0)
    mean = np.nanmean(per_tree, 0)
    se = np.nanstd(per_tree, 0, ddof=1) / np.sqrt(np.sum(np.isfinite(per_tree), 0))
    tg = np.linspace(1e-6, 1 - 1e-6, 4001)
    ptil = rec.share_at(tg).sum(-1)
    dens = [kernel_density_t(j, N, N, rec.Lam_T, tg) for j in range(1, N + 1)]
    pred = np.array([np.trapezoid(dd * ptil, tg) for dd in dens])
    tmean = np.array([np.trapezoid(dd * tg, tg) for dd in dens])
    p_off, p_on = PULSE["shares"][0][0], PULSE["shares"][1][0]
    on0, on1 = PULSE["comp_knots"][1], PULSE["comp_knots"][2]

    fig, (top, ax) = plt.subplots(2, 1, figsize=(5.8, 5.8),
                                  gridspec_kw={"height_ratios": [1.15, 2.0], "hspace": 0.55})
    top.axvspan(on0, on1, color=BAND, zorder=0, lw=0)
    top.text((on0 + on1) / 2, 1.02, f"signal on: {p_on:.0%} of edits", transform=top.get_xaxis_transform(),
             ha="center", va="bottom", color=INK2, fontsize=8.5)
    top.text(0.2, 1.02, f"off: {p_off:.0%}", transform=top.get_xaxis_transform(), ha="center",
             va="bottom", color=MUTED, fontsize=8.5)
    top.text(0.85, 1.02, f"off: {p_off:.0%}", transform=top.get_xaxis_transform(), ha="center",
             va="bottom", color=MUTED, fontsize=8.5)
    for j, dd in enumerate(dens):
        top.plot(tg, dd, color=RAMP6[j], lw=1.6, zorder=3)
        pk = int(np.argmax(dd))
        top.text(tg[pk], dd[pk] + 0.12, str(j + 1), color=INK2, fontsize=8, ha="center")
    top.set_xlim(0, 1)
    top.set_ylim(0, max(dd.max() for dd in dens) * 1.3)
    top.set_yticks([])
    top.set_xlabel("when the slot was written (fraction of the experiment)", fontsize=9)
    top.set_ylabel("how often", fontsize=9)
    style(top, grid=False)
    top.spines["left"].set_visible(False)
    top.set_title("c   A signal pulse is written into slot order", loc="left", fontsize=10.5,
                  pad=18)

    xs = np.arange(1, N + 1)
    ax.axhline(p_on, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.axhline(p_off, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(N + 0.45, p_on, "if written\nentirely on", color=MUTED, fontsize=8, va="center")
    ax.text(N + 0.45, p_off, "entirely off", color=MUTED, fontsize=8, va="center")
    ax.scatter(xs, pred, marker="_", s=420, color=LINE, lw=2, zorder=3, label="prediction")
    for j in range(N):
        ax.errorbar(xs[j], mean[j], yerr=se[j], fmt="o", ms=7, color=RAMP6[j], ecolor=RAMP6[j],
                    mec="white", mew=0.8, zorder=4)
    ax.scatter([], [], marker="o", s=40, color=RAMP6[3], label="simulated (mean ± se across trees)")
    ax.set_xticks(xs)
    ax.set_xlim(0.4, N + 0.4)
    ax.set_ylim(0, p_on * 1.15)
    ax.set_xlabel("slot on the tape (1 = written first)")
    ax.set_ylabel("share of slots holding the signal symbol")
    ax.legend(frameon=False, loc="upper left", fontsize=8.5)
    style(ax)
    fig.text(0.0, -0.2, f"tapes that filled all {N} slots ({np.mean(full_frac):.0%} of tapes); "
             f"200 library trees of {N_TIPS} cells, 30 tapes; constant editing rate,\n"
             f"$\\Lambda_T$ = {rec.Lam_T} edits per tape; top = when each slot was written on "
             f"those tapes ($K_{{j,6}}$ of notes §S4.9.5)",
             color=MUTED, fontsize=8, ha="left", va="top", transform=ax.transAxes)
    save(fig, "figS5c_signal_in_slot_order")
    numbers["c"] = {"pulse": PULSE, "slot": xs.tolist(), "simulated": np.round(mean, 4).tolist(),
                    "se": np.round(se, 5).tolist(), "prediction": np.round(pred, 4).tolist(),
                    "mean_write_time": np.round(tmean, 4).tolist(),
                    "fraction_of_tapes_full": float(np.mean(full_frac)),
                    "z": np.round((mean - pred) / se, 2).tolist()}


def panel_d(trees, numbers):
    rec = ed.Recorder(k=30, shares=((0.0,),), **RATE)
    N, k = rec.N, rec.k
    edges = np.arange(0.0, 1.0 + RATE_BIN / 2, RATE_BIN)
    B = edges.size - 1
    lo, hi = edges[:-1], edges[1:]
    tg = np.linspace(0.0, 1.0, 41)
    E_t = np.zeros((len(trees), B)); O_t = np.zeros((len(trees), B)); X_t = np.zeros((len(trees), B))
    depth_t = np.zeros((len(trees), tg.size))
    ss = np.random.SeedSequence([SEED, 7]).spawn(len(trees))
    for ti, (parent, branch) in enumerate(trees):
        out = ed.decorate(parent, branch, N_TIPS, rec, np.random.default_rng(ss[ti]))
        t0, t1 = out["t0"], out["t"]
        d0 = out["bslot0"].astype(int)                                 # depth on arrival (nn, k)
        ne = out["bcount"].astype(int)
        bt = out["btime"].astype(float)
        last = np.take_along_axis(bt, np.maximum(ne - 1, 0)[..., None], axis=2)[..., 0]
        start = np.broadcast_to(t0[:, None], d0.shape)
        close = np.where(d0 >= N, start,                               # arrived full: never open
                         np.where(d0 + ne >= N, last, t1[:, None]))    # filled here: closes then
        a, b = start.ravel(), close.ravel()
        O_t[ti] = np.clip(np.minimum(b[:, None], hi) - np.maximum(a[:, None], lo), 0, None).sum(0)
        a2 = start.ravel(); b2 = np.broadcast_to(t1[:, None], d0.shape).ravel()
        X_t[ti] = np.clip(np.minimum(b2[:, None], hi) - np.maximum(a2[:, None], lo), 0, None).sum(0)
        E_t[ti] = np.histogram(bt[np.isfinite(bt)], bins=edges)[0]
        for g, tt in enumerate(tg):                                    # filled slots at time tt
            span = (t0 <= tt) & ((t1 > tt) | ((tt >= 1.0) & (t1 >= 1.0)))
            dep = d0[span] + np.sum(bt[span] <= tt, axis=-1)
            depth_t[ti, g] = dep.mean()

    def ratio(E, D):
        R = E.sum(0) / D.sum(0)
        se = np.sqrt(np.sum((E - R * D) ** 2, 0)) / D.sum(0)          # trees as clusters
        return R, se

    r_open, se_open = ratio(E_t, O_t)
    r_all, _ = ratio(E_t, X_t)
    mid = 0.5 * (lo + hi)
    programmed = rec.Lam_T * rec.omega[np.clip(np.searchsorted(rec.kappa, mid, side="right"), 1,
                                                rec.omega.size) - 1]
    dmean = depth_t.mean(0)
    dse = depth_t.std(0, ddof=1) / np.sqrt(len(trees))
    tf = np.linspace(0, 1, 801)

    def mean_depth(mu):
        return (np.minimum(np.arange(60), N)[None, :] *
                poisson.pmf(np.arange(60)[None, :], np.asarray(mu)[:, None])).sum(1)

    pred = mean_depth(rec.Lam_T * rec.W_of(tf))
    pred_const = mean_depth(rec.Lam_T * tf)
    pred_g = mean_depth(rec.Lam_T * rec.W_of(tg))

    fig, (top, ax) = plt.subplots(2, 1, figsize=(5.8, 6.6), sharex=True,
                                  gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.18})
    for axx in (top, ax):
        axx.axvspan(0.25, 0.35, color=BAND, lw=0, zorder=0)
    top.text(0.30, 1.01, "pause", transform=top.get_xaxis_transform(), ha="center", va="bottom",
             color=INK2, fontsize=8.5)
    top.text(0.475, 1.01, "burst", transform=top.get_xaxis_transform(), ha="center", va="bottom",
             color=INK2, fontsize=8.5)
    kx = np.repeat(rec.kappa, 2)[1:-1]
    ky = np.repeat(rec.Lam_T * rec.omega, 2)
    top.plot(kx, ky, color=LINE, lw=2, zorder=3, label="programmed rate")
    top.errorbar(mid, r_open, yerr=se_open, fmt="o", ms=4.5, color=RAMP6[2], ecolor=RAMP6[2],
                 elinewidth=1, mec="white", mew=0.7, zorder=4, label="realised, per open tape")
    top.scatter(mid, r_all, s=16, facecolors="white", edgecolors=MUTED, lw=1, zorder=4,
                label="realised, per tape (open or full)")
    top.set_ylabel("editing rate\n(edits per tape per experiment)", fontsize=9)
    top.set_ylim(0, max(ky.max(), r_open.max()) * 1.5)
    top.legend(frameon=False, loc="upper right", fontsize=8, bbox_to_anchor=(1.0, 1.0))
    style(top)
    top.set_title("d   The editing rate can be programmed over time", loc="left", fontsize=10.5,
                  pad=16)

    ax.axhline(N, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(0.01, N + 0.08, f"tape full ({N} slots)", color=MUTED, fontsize=8.5)
    ax.plot(tf, pred_const, color=MUTED, lw=1.2, ls=(0, (2, 2)), zorder=2,
            label="constant rate, same total")
    ax.plot(tf, pred, color=LINE, lw=2, zorder=3, label="prediction")
    ax.errorbar(tg, dmean, yerr=dse, fmt="o", ms=4.5, color=RAMP6[2], ecolor=RAMP6[2],
                elinewidth=1, mec="white", mew=0.7, zorder=4, label="simulated (mean ± se)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, N + 0.5)
    ax.set_xlabel("time (fraction of the experiment)")
    ax.set_ylabel("filled slots per tape\n(average along lineages)", fontsize=9)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    style(ax)
    fig.text(0.0, -0.25, f"200 library trees of {N_TIPS} cells, 30 tapes, lineage symbols only; "
             f"$\\Lambda_T$ = {rec.Lam_T} edits per tape; rate bins {RATE_BIN} wide\n"
             f"⚠ uses the simulator's true edit times: shows the rate can be imposed, not that it "
             f"can be recovered from sequenced tapes (notes §S4.7)",
             color=MUTED, fontsize=8, ha="left", va="top", transform=ax.transAxes)
    save(fig, "figS5d_programmed_rate")
    live = se_open > 0
    numbers["d"] = {"schedule": {"knots": rec.kappa.tolist(),
                                 "rate_edits_per_tape": (rec.Lam_T * rec.omega).round(4).tolist()},
                    "bin_mid": mid.round(4).tolist(), "programmed": programmed.round(4).tolist(),
                    "realised_open": r_open.round(4).tolist(), "se_open": se_open.round(4).tolist(),
                    "realised_all": r_all.round(4).tolist(),
                    "z_open": np.where(live, (r_open - programmed) / np.where(live, se_open, 1),
                                       0.0).round(2).tolist(),
                    "pause_edits": int(E_t[:, (mid > 0.25) & (mid < 0.35)].sum()),
                    "depth_time": tg.round(3).tolist(), "depth_sim": dmean.round(4).tolist(),
                    "depth_se": dse.round(4).tolist(), "depth_pred": pred_g.round(4).tolist(),
                    "z_depth": np.where(dse > 0, (dmean - pred_g) / np.where(dse > 0, dse, 1),
                                        0.0).round(2).tolist()}


def main():
    trees = []
    for rho, th in TREE_CELLS:
        trees += ed.load_library_trees(N_TIPS, rho, th)[:20]
    numbers = {"seed": SEED, "trees": len(trees), "n": N_TIPS}
    panel_a(numbers)
    panel_b(trees, numbers)
    panel_c(trees, numbers)
    panel_d(trees, numbers)
    (RES / "figS5_numbers.json").write_text(json.dumps(numbers, indent=1, default=float))
    print(json.dumps({k: numbers[k] for k in ("b", "c")}, default=float)[:1500])


if __name__ == "__main__":
    main()
