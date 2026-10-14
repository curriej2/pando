"""Validation figure 1: our engine's likelihood vs SciPhy's own Java likelihood (check V6).

Replots results/validate_sciphy.json (written by 06_validate_sciphy.py, 2026-10-12); computes nothing new.
One panel: total log-likelihood of each simulated dataset (30 tapes summed), SciPhy on x, ours on y.
Only the V6a coding (blank tapes entered as observed blanks, SciPhy's own coding) is plotted.
Our engine is evaluated at the special case SciPhy computes: no technical dropout (pi = 0), no silencing
(Gam_T = 0), blanks read as blank (no depth-0 censoring), editing rate constant in time (per-tape rates
r_z allowed, as in SciPhy), one constant composition xi, blank founder. README "Session 25".
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# palette: reference categorical slot 1 (validated, light mode), as in src/07
BLUE = "#2a78d6"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})

d = json.load(open(ROOT / "results" / "validate_sciphy.json"))
S = d["summary"]["V6a_observed_blank"]
obs = sorted((r for r in d["rows"] if r["coding"] == "V6a_observed_blank"), key=lambda r: r["config"])
assert len(obs) == S["configs"] == 50
n_tapes = S["tape_comparisons"] // S["configs"]
worst = max(r["worst_rel"] for r in obs)

xs = np.array([-r["sum_ll_sciphy"] for r in obs])
ys = np.array([-r["sum_ll_engine"] for r in obs])
lo, hi = xs.min() / 1.4, xs.max() * 1.4

fig, ax = plt.subplots(figsize=(6.2, 6.4))
ax.plot([lo, hi], [lo, hi], color=INK2, lw=1, ls="--", zorder=1, label="identity line (y = x)")
ax.scatter(xs, ys, s=40, color=BLUE, edgecolors=SURF, lw=1, zorder=3,
           label=f"one simulated dataset ({n_tapes} tapes)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel(f"SciPhy's code: −log-likelihood, {n_tapes} tapes summed (nats)")
ax.set_ylabel("our engine: −log-likelihood (nats)")
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.legend(loc="upper left", fontsize=9, frameon=False, handletextpad=0.3)
ax.text(0.97, 0.04, f"largest per-tape difference\n{worst:.1e} (relative), over {S['tape_comparisons']:,} tapes",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=INK)
ax.set_title("Our likelihood reproduces SciPhy's on its own model", loc="left")
ax.text(0, -0.13, f"{len(obs)} simulated datasets: trees of 25 / 64 / 210 cells, 1–8 edits/tape, 0–79% of tapes\n"
        "full, 0–47% blank. Our engine set to SciPhy's model: no dropout, no silencing, blanks\n"
        "read as blank, editing rate constant in time, one symbol composition.",
        transform=ax.transAxes, fontsize=8.5, color=INK2, va="top")

fig.tight_layout()
f = OUT / "fig_val1_sciphy.png"
fig.savefig(f, dpi=200, bbox_inches="tight")
print(f"wrote {f}; worst per-tape rel {worst:.2e} over {S['tape_comparisons']} tapes")
