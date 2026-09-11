#!/usr/bin/env python3
r"""
Fig 4e -- structured dropout WITHOUT events.  Three standalone panels.

THE ARGUMENT.  Every earlier route to "dropout is heritable" had to define an
event: pick a clade, score it, threshold it, collapse overlaps, control
multiplicity.  These two statistics need none of that -- no clades, no
thresholds, no attribution -- which is exactly what makes them usable as a
SIMULATOR CALIBRATION TARGET: push simulated data through the identical code and
match the curve.

⚠⚠ THE CONFOUND BOTH DEFEAT, and it must be visible in the caption.  Lineage
relatedness is read from the edit data, and dropout decides which edits are
readable, so two cells that both lack tapes 1-50 look related BECAUSE their
dropout matches.  ⇒ DISJOINT TAPE SPLIT: relatedness from half A, dropout from
half B, over 10 random splits.  Panel a's null is the same statistic with the
cell labels permuted on the half-B residuals.

PANEL a  variogram: mean Pearson-residual co-deviation vs lineage relatedness.
PANEL b  the conditional read-off: observed dropout rate against what fraction of
         a cell's relatives lack the tape, WITHIN strata of the model's own
         predicted probability -- so cell and tape quality are held fixed.
PANEL c  held-out likelihood gain, nats per cell, by k and arm.

⚠ Panel a is plotted as obs - null, not obs.  The pair-weighted mean of the curve
is pinned at the all-pairs reference (~0) because sum_c r_cz = 0 within a clone,
so the negative low bins are the ARITHMETIC COMPLEMENT of the positive high bins,
not a second finding.  The content is the SLOPE.  The zero line is therefore
drawn as the reference it is, and the panel is annotated to say so.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Initial", "Subclone", "Mouse2", "Mouse1", "Mouse3"]
NICE = {"Initial": "Pre-TX", "Mouse1": "Mouse 1", "Mouse2": "Mouse 2",
        "Mouse3": "Mouse 3", "Subclone": "Subclone"}
COL = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
       "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
CLONE = {"Mouse3": 110, "Mouse1": 36, "Mouse2": 76, "Initial": 7, "Subclone": 2}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
NEUT, NEUT2 = "#cfcec6", "#8f8e86"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

# ---------------------------------------------------------------- panel a ----
fig, ax = plt.subplots(figsize=(6.6, 4.4))
for arm in ARMS:
    f = RES / f"variogram_{arm}_pooled.npz"
    if not f.exists():
        print(f"  [skip] {f.name}"); continue
    z = np.load(f)
    x = np.nanmean(z["relmeans"], 0)
    d = z["obs"] - z["null"]
    y = np.nanmean(d, 0)
    e = np.nanstd(d, 0, ddof=1) / np.sqrt(d.shape[0])
    ax.plot(x, np.nanmean(z["null"], 0), color=NEUT2, lw=1.0, zorder=2)  # drawn, not claimed
    ax.errorbar(x, y, yerr=e, color=COL[arm], lw=1.8, marker="o", ms=3.4,
                capsize=0, elinewidth=1.0, zorder=3, label=NICE[arm])
ax.plot([], [], color=NEUT2, lw=1.0, label="null (cell labels permuted)")
ax.set_xlabel("lineage relatedness  =  mean shared prefix depth per half-A tape")
ax.set_ylabel("dropout co-deviation, obs $-$ null\n"
              r"$\frac{1}{|B|}\sum_{z\in B} r^{*}_{cz}\,r^{*}_{c'z}$")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.text(0.985, 0.03,
        "relatedness from tape half A · dropout from half B · within clone\n"
        r"$r^{*}=(X-\tilde p)/\sqrt{\tilde p(1-\tilde p)}$,  "
        r"$\tilde p=\sigma(\alpha_c+\beta_z+\gamma_{Cz})$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6, color=INK2)
fig.savefig(FIG / "fig4e_a_variogram.png"); plt.close(fig)
print("wrote figures/fig4e_a_variogram.png")

# ---------------------------------------------------------------- panel b ----
ARM_B, TAG_B = "Mouse2", "c76"
j = json.loads((RES / f"profilepred_{ARM_B}_{TAG_B}.json").read_text())
tb = j["table"]; rows = tb["rows"]; FL = tb["f_labels"]
xs = np.arange(len(FL))
fig, ax = plt.subplots(figsize=(6.6, 4.4))
shades = ["#cfd9ea", "#8fb0dc", "#4f83c8", "#1f4f8f"]
for i, r in enumerate(rows):
    ax.axhline(100 * r["p_mid"], color=shades[i], lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.plot(xs, [100 * v for v in r["rates"]], color=shades[i], lw=1.9,
            marker="o", ms=4.5, zorder=3,
            label=f"model says {100*r['p_mid']:.0f}%")
    for xi, (v, c) in enumerate(zip(r["rates"], r["counts"])):
        if not c:
            continue
        if xi == len(FL) - 1:
            # ⚠ all four series converge at "all" (92-99%), so a count placed at the
            # data's own y always collides.  Put them on a FIXED ladder in the right
            # margin, colour-coded to the series -- the colour carries the mapping.
            ax.text(xi + 0.12, 58 + i * 9, f"n={c:,}", ha="left", va="center",
                    fontsize=6.2, color=shades[i])
        elif xi == 0 or c < 1000:
            ax.text(xi, 100 * v + 3.2, f"n={c:,}", ha="center", fontsize=6.2, color=shades[i])
ax.set_xticks(xs); ax.set_xticklabels(FL); ax.set_xlim(-0.45, len(FL) + 0.32)
ax.set_xlabel("fraction of the cell's $k=20$ nearest relatives that lack this tape\n"
              "(relatives found on half A, self excluded; tape read on half B)")
ax.set_ylabel("observed dropout rate (%)")
ax.set_ylim(-4, 108)
ax.legend(frameon=False, fontsize=8.5, ncol=4, loc="lower center",
          bbox_to_anchor=(0.5, 1.0), columnspacing=1.4, handlelength=1.6,
          title="held fixed — the model's own predicted dropout for the cell and tape",
          title_fontsize=8.5)
ax.text(0.985, 0.03, f"{NICE[ARM_B]}, clone {j['clone']} · {j['n_cells']:,} cells\n"
        "dashed line = the model's own prediction for that stratum",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6, color=INK2)
fig.savefig(FIG / "fig4e_b_prediction.png"); plt.close(fig)
print("wrote figures/fig4e_b_prediction.png")

# ---------------------------------------------------------------- panel c ----
fig, ax = plt.subplots(figsize=(6.0, 4.2))
for arm in ARMS:
    f = RES / f"profilepred_{arm}_c{CLONE[arm]}.json"
    if not f.exists():
        print(f"  [skip] {f.name}"); continue
    d = json.loads(f.read_text())
    ks = d["k_list"]
    y, e = [], []
    for k in ks:
        o = np.array([r["nats_per_cell"] for r in d["runs"] if r["k"] == k and r["kind"] == "obs"])
        u = np.array([r["nats_per_cell"] for r in d["runs"] if r["k"] == k and r["kind"] == "null"])
        y.append((o - u).mean()); e.append((o - u).std(ddof=1) / np.sqrt(o.size))
    ax.errorbar(ks, y, yerr=e, color=COL[arm], lw=1.8, marker="o", ms=4.2,
                capsize=0, elinewidth=1.0,
                label=f"{NICE[arm]} c{CLONE[arm]} ({d['n_cells']:,} cells)")
ax.axhline(0, color=NEUT2, lw=1.0)
ax.set_xscale("log"); ax.set_xticks([5, 20, 50]); ax.set_xticklabels(["5", "20", "50"])
ax.minorticks_off()
ax.set_xlabel("$k$, number of nearest relatives used")
ax.set_ylabel("held-out likelihood gain, nats per cell\n(observed $-$ null)")
ax.legend(frameon=False, fontsize=8.0, ncol=3, loc="upper center",
          bbox_to_anchor=(0.5, -0.20), columnspacing=1.2, handlelength=1.6)
ax.text(0.985, 0.03, "one clone per arm, the largest · 5 tape splits\n"
        r"$\operatorname{logit}\Pr(X_{cz}{=}1)=\eta_{cz}+w\,u_{cz}$, one free scalar $w$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6, color=INK2)
fig.savefig(FIG / "fig4e_c_nats.png"); plt.close(fig)
print("wrote figures/fig4e_c_nats.png")
