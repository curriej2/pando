#!/usr/bin/env python3
r"""
13_fig_switch_evidence.py -- Fig S4: how far can a lineage-specific switch be detected?

Reads results/switch_evidence.json (written by 12) and draws two standalone panels.
Nothing is computed here beyond medians and ranges over the 30 (n, rho) settings.

  S4a  the AGREED headline: fraction of all switch locations callable, against switch
       duration, one line per fold change.
  S4b  the same result read by WHEN the switch starts: callable fraction by start-time
       decile (bottom), under a strip showing where switch locations fall (top).
       ⚠ Added after the run (README "Session 17"): switch locations are uniform along
       tree length, and 43-92% of that length sits in the last 30% of the experiment,
       so S4a is mostly a statement about late switches.

Fixed in both: k = 30 tapes/cell, p1 = 0.3, edits placed only to their branch
(the realistic mode), threshold 3 + ln H nats.  Line = median over the 30 settings,
band = full range (the settings are a design grid, not a distribution).
"""
from __future__ import annotations

import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
FIG = _HERE.parent / "figures"
RES = _HERE.parent / "results"

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
BAR = "#c9c7c1"
# fold and duration are magnitudes: one sequential blue ramp, light -> dark,
# starting at step 250 so the lightest line keeps 2:1 against white
FOLD_COL = {2.0: "#86b6ef", 5.0: "#3987e5", 15.0: "#1c5cab", 24.0: "#0d366b"}
DUR_COL = {0.5: "#86b6ef", 1.0: "#2a78d6", 2.0: "#0d366b"}
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "pdf.fonttype": 42,
})

MODE_BRANCH, P_HEAD, K_NOTE = 1, 0.3, 30


def arr(x):
    return np.array(json.loads(json.dumps(x).replace("null", "NaN")), dtype=float)


def style(ax):
    ax.minorticks_off()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#e6e5e1", lw=0.6, zorder=0)


def save(fig, stem):
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG / stem}.{{pdf,png}}")


def main():
    J = json.loads((RES / "switch_evidence.json").read_text())
    S, cells = J["settings"], J["cells"]
    lam = np.array(S["lam_e_edits_per_tape"]); folds = S["fold_F"]
    p = S["p1"].index(P_HEAD)
    assert S["k_head"] == K_NOTE
    H = np.array([arr(c["frac_callable_k30_tree_mean"]) for c in cells])[:, MODE_BRANCH, :, :, p]
    TB = np.array([arr(c["frac_callable_k30_by_time"]) for c in cells])[:, :, MODE_BRANCH, :, :, p]
    share = np.array([arr(c["time_bin_share"]) for c in cells])
    edges = np.array(cells[0]["time_bins_left_edge"] + [1.0])
    mid = 0.5 * (edges[:-1] + edges[1:])
    note = (f"k = {K_NOTE} tapes per cell, signal share when on p₁ = {P_HEAD}, "
            "edits placed only to their branch\n"
            "line = median over 30 (clone size, capture fraction) settings; band = full range")
    numbers = {"settings": {"k": K_NOTE, "p1": P_HEAD, "mode": "branch_only"}}

    # ---- S4a: the agreed headline -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    x = np.log2(lam)
    numbers["a"] = {"lam_e": lam.tolist(), "series": {}}
    for f, F in enumerate(folds):
        y = H[:, :, f]
        med, lo, hi = np.median(y, 0), y.min(0), y.max(0)
        col = FOLD_COL[F]
        ax.fill_between(x, lo, hi, color=col, alpha=0.15, lw=0, zorder=1)
        ax.plot(x, med, color=col, lw=2, marker="o", ms=5, zorder=3, label=f"{F:g}×")
        ax.annotate(f"{F:g}×", (x[-1], med[-1]), xytext=(7, 0), textcoords="offset points",
                    va="center", color=INK2, fontsize=9)
        numbers["a"]["series"][f"{F:g}x"] = {"median": med.round(4).tolist(),
                                             "min": lo.round(4).tolist(),
                                             "max": hi.round(4).tolist()}
    ax.axhline(0.5, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
    ax.text(x[0], 0.515, "pre-set positive bar (50% of locations)", color=MUTED, fontsize=8.5)
    ax.set_xticks(x, [f"{v:g}\n({100 * v / S['Lam_T_edits_per_tape']:.1f}%)" for v in lam])
    ax.set_xlim(x[0] - 0.2, x[-1] + 0.45)
    ax.set_ylim(0, 1)
    ax.set_xlabel("switch duration, edits per tape (share of the experiment)")
    ax.set_ylabel("fraction of all switch locations callable")
    ax.legend(title="fold change", frameon=False, loc="upper left", fontsize=8.5,
              title_fontsize=8.5, bbox_to_anchor=(0.0, 0.93))
    style(ax)
    ax.set_title("a   Long, strong switches can be seen; short or weak ones cannot",
                 loc="left", fontsize=10.5)
    fig.text(0.0, -0.24, note, color=MUTED, fontsize=8, ha="left", va="top",
             transform=ax.transAxes)
    save(fig, "figS4a_callable_by_duration")

    # ---- S4b: read by start time ------------------------------------------------------
    F24 = folds.index(24.0)
    fig, (top, ax) = plt.subplots(2, 1, figsize=(5.8, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [1, 2.6], "hspace": 0.12})
    smed, slo, shi = np.median(share, 0), share.min(0), share.max(0)
    top.bar(mid, smed, width=0.085, color=BAR, zorder=2)
    top.errorbar(mid, smed, yerr=[smed - slo, shi - smed], fmt="none", ecolor=MUTED,
                 elinewidth=0.9, capsize=2, zorder=3)
    top.set_ylabel("share of all\nswitch locations", fontsize=9)
    top.set_ylim(0, max(0.55, shi.max() * 1.05))
    style(top)
    top.set_title("b   The early history is visible; the late history is not",
                  loc="left", fontsize=10.5)
    numbers["b"] = {"start_bin_mid": mid.round(3).tolist(), "fold": 24,
                    "location_share": {"median": smed.round(4).tolist(),
                                       "min": slo.round(4).tolist(),
                                       "max": shi.round(4).tolist()},
                    "series": {}}
    for L in (0.5, 1.0, 2.0):
        d = int(np.where(lam == L)[0][0])
        y = TB[:, :, d, F24]
        med, lo, hi = np.nanmedian(y, 0), np.nanmin(y, 0), np.nanmax(y, 0)
        col = DUR_COL[L]
        ax.fill_between(mid, lo, hi, color=col, alpha=0.15, lw=0, zorder=1)
        ax.plot(mid, med, color=col, lw=2, marker="o", ms=5, zorder=3,
                label=f"{L:g} edits/tape ({100 * L / S['Lam_T_edits_per_tape']:.0f}%)")
        numbers["b"]["series"][f"lam_e_{L:g}"] = {"median": med.round(4).tolist(),
                                                 "min": lo.round(4).tolist(),
                                                 "max": hi.round(4).tolist()}
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, 1)
    ax.set_xlabel("when the switch starts (fraction of the experiment)")
    ax.set_ylabel("fraction callable (fold change 24×)")
    ax.legend(title="switch duration", frameon=False, loc="lower left", fontsize=8.5,
              title_fontsize=8.5)
    style(ax)
    fig.text(0.0, -0.18, note, color=MUTED, fontsize=8, ha="left", va="top",
             transform=ax.transAxes)
    save(fig, "figS4b_callable_by_time")

    (RES / "figS4_numbers.json").write_text(json.dumps(numbers, indent=1))
    print(f"wrote {RES / 'figS4_numbers.json'}")


if __name__ == "__main__":
    main()
