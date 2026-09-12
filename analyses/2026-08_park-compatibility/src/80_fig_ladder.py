#!/usr/bin/env python3
r"""
FIGS A, B, C -- the talk set for "does lineage information improve dropout prediction".
Three standalone PNGs.  Source: ladder_{arm}_mincl100_k20.json (77) and
relvsrnd_{arm}_mincl100_k20.json (79).

EVERY NUMBER PLOTTED, AND WHERE IT COMES FROM
---------------------------------------------
Each ladder replicate stores l_0..l_4, the MEAN HELD-OUT LOG-LIKELIHOOD PER TEST
ENTRY under each nested model, plus l_4^null for the matched null.  All five are
computed on the SAME held-out entries, so differences between them are like-for-like.

  l_j = (1/|V|) * sum over held-out (c,z) of [ X log p_j + (1-X) log(1-p_j) ]

FIG A -- cost against benefit
  y, gain of rung j   = K * (l_j - l_{j-1})        [K = 166 tapes]
      Intuition: l is per ENTRY, so multiplying by the number of tapes puts it on a
      per-CELL scale, the unit every earlier result in this project uses.  For j=4
      we use K * (l_4 - l_4^null) instead, because the raw rung is biased: the M3
      score equation forces the clone's residuals to sum to zero, so a neighbour set
      that excludes the cell carries -1/(n_C - 1) of the cell's own residual.  The
      null carries the identical term, so the difference is the clean quantity.
  x, parameters added at rung j.  Derived from the model definitions, not counted:
      M1  mu_C per clone, over one global mu          ->  nC - 1
      M2  alpha_c per cell, one constraint per clone  ->  n - nC
      M3  beta_z per tape per clone, absorbing mu_C   ->  nC * (|B| - 1)
      M4  one scalar w                                ->  1
      |B| = K//2 = 83, because the model is fitted on half the tapes (the other half
      defines relatedness and must stay disjoint from it).
      Intuition: the x-axis is "how many knobs did you turn", the y-axis is "how much
      better did the prediction get".  A point at the far left and high up is buying
      a lot for almost nothing, which is the whole argument.

FIG B -- the honest zoom
  The saturated log-likelihood is 0 (with p = X exactly, every term is log 1), so
  -l_0 is the TOTAL deviance available to explain, and every share below is a
  fraction of it.
      technical share  = 1 - l_3/l_0          (all of clone + cell + tape-in-clone)
      lineage share    = (l_4 - l_4^null) / (-l_0)
      unexplained      = 1 - technical - lineage
      Q                = lineage / (1 - technical) = (l_4 - l_4^null) / (-l_3)
      Intuition: the first three rungs dominate, and a plain stacked bar would make
      the lineage slice look trivial.  Q re-bases the question onto what is LEFT
      after every technical correction, which is the question actually being asked.
  ⚠ The M1/M2/M3 rungs are NOT shown separately here: M2 is negative on two arms and
  a negative segment in a stacked bar is unreadable.  Fig A carries the breakdown.

FIG C -- relatives against strangers
  From 79.  For each HELD-OUT entry (c,z): f = the fraction of c's k neighbours that
  are missing tape z, computed over the neighbours' TRAINING entries only.  Two
  neighbour sets: the k nearest relatives (found on tape half A, self excluded) and
  k random clone-mates (self excluded) -- the latter IS the ladder's null, so the
  flat series in this panel is the null drawn from data rather than asserted.
  Entries are binned by f WITHIN quartiles of the model's own predicted probability
  p~, so cell quality and tape quality are held fixed and cannot explain the contrast.
      plotted rate = (number missing in that bin) / (number of entries in that bin)
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Subclone", "Initial", "Mouse2", "Mouse3", "Mouse1"]
NICE = {"Initial": "Pre-TX", "Mouse1": "Mouse 1", "Mouse2": "Mouse 2",
        "Mouse3": "Mouse 3", "Subclone": "Subclone"}
COL = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
       "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
NEUT, NEUT2 = "#cfcec6", "#8f8e86"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GRID,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.dpi": 160, "savefig.bbox": "tight"})

D = {}
for a in ARMS:
    d = json.loads((RES / f"ladder_{a}_mincl100_k20.json").read_text())
    L = np.array([r["ll"] for r in d["runs"]]); L4n = np.array([r["ll_null4"] for r in d["runs"]])
    K, n, nC = d["K"], d["n_cells"], d["n_clones"]; B = K // 2
    D[a] = dict(K=K, n=n, nC=nC,
                pars=[nC - 1, n - nC, nC * (B - 1), 1],
                gain=[K * (L[:, 1] - L[:, 0]).mean(), K * (L[:, 2] - L[:, 1]).mean(),
                      K * (L[:, 3] - L[:, 2]).mean(), K * (L[:, 4] - L4n).mean()],
                se=[0, 0, 0, K * (L[:, 4] - L4n).std(ddof=1) / np.sqrt(len(L4n))],
                tech=(1 - L[:, 3] / L[:, 0]).mean(),
                lin=((L[:, 4] - L4n) / (-L[:, 0])).mean(),
                Q=((L[:, 4] - L4n) / (-L[:, 3])).mean())

# ----------------------------------------------------------------- FIG A ----
# ⚠ First version was a log scatter of gain against parameters with the four rungs
# of an arm joined by a line.  Two faults: the line implied a TRAJECTORY through a
# space the rungs do not move through, and M3's +48 nats squashed M4's +1.5..+8.3
# against the axis.  Grouped bars put the rungs side by side, show the NEGATIVE M2
# honestly, and carry the parameter cost as a printed number under each group --
# which is where it belongs, since "1" and "38,604" do not share a usable axis.
RUNG = ["M1\nwhich clone", "M2\ncell capture", "M3\ntape, within clone", "M4\nlineage"]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
W = 0.15
ax.axhline(0, color=INK2, lw=1.0, zorder=3)
for i, a in enumerate(ARMS):
    d = D[a]
    xs = np.arange(4) + (i - 2) * W
    ax.bar(xs, d["gain"], width=W, color=COL[a], edgecolor=SURFACE, linewidth=0.6,
           zorder=2, label=NICE[a])
    ax.errorbar(xs[3], d["gain"][3], yerr=d["se"][3], color=INK, lw=1.0,
                capsize=2.0, zorder=4)
for j in range(4):
    lo = min(D[a]["pars"][j] for a in ARMS); hi = max(D[a]["pars"][j] for a in ARMS)
    lab = "1 parameter" if j == 3 else f"{lo:,}–{hi:,} parameters"
    ax.text(j, -6.5, lab, ha="center", fontsize=8.0,
            color=INK if j == 3 else INK2, weight="bold" if j == 3 else "normal")
ax.set_xticks(np.arange(4)); ax.set_xticklabels(RUNG, fontsize=9)
ax.tick_params(axis="x", length=0, pad=18)
ax.set_ylabel("held-out gain, nats per cell")
ax.set_ylim(-8.5, 52)
ax.legend(frameon=False, fontsize=8.4, ncol=5, loc="upper center",
          bbox_to_anchor=(0.5, 1.10), columnspacing=1.2, handlelength=1.2)
ax.text(0.015, 0.95, "M2 goes NEGATIVE where the cell filter is strict\n"
        "(Subclone and Pre-TX keep cells with \u2265100 tapes, the mice \u226520)\n\n"
        "M4 is scored against its own null:\n"
        "k random clone-mates, self excluded",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.8, color=INK2,
        linespacing=1.45)
fig.savefig(FIG / "fig6a_cost_benefit.png"); plt.close(fig)
print("wrote figures/fig6a_cost_benefit.png")

# ----------------------------------------------------------------- FIG B ----
fig, ax = plt.subplots(figsize=(7.0, 3.6))
yy = np.arange(len(ARMS))[::-1]
for i, a in enumerate(ARMS):
    d = D[a]; y = yy[i]
    t, l = 100 * d["tech"], 100 * d["lin"]
    ax.barh(y, t, color=NEUT, height=0.62, edgecolor=SURFACE, linewidth=0.8)
    ax.barh(y, l, left=t, color=COL[a], height=0.62, edgecolor=SURFACE, linewidth=0.8)
    ax.barh(y, 100 - t - l, left=t + l, color="#efeee8", height=0.62,
            edgecolor=SURFACE, linewidth=0.8)
    ax.text(t + l + 1.6, y, f"{100*d['Q']:.1f}%", va="center", fontsize=9.2,
            color=COL[a], weight="bold")
ax.set_yticks(yy); ax.set_yticklabels([NICE[a] for a in ARMS])
ax.set_xlim(0, 100); ax.set_xlabel("share of total deviance (%)")
ax.set_ylim(-0.75, len(ARMS) - 0.25)
h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (NEUT, COL["Subclone"], "#efeee8")]
ax.legend(h, ["explained by clone + cell + tape", "explained by lineage",
              "still unexplained"], frameon=False, fontsize=8.2, ncol=3,
          loc="upper center", bbox_to_anchor=(0.5, 1.16), handlelength=1.3)
ax.text(0.5, -0.30, "the percentage is the lineage share of what the technical model "
        "left unexplained", transform=ax.transAxes, ha="center", va="top",
        fontsize=8.0, color=INK2)
fig.savefig(FIG / "fig6b_deviance.png"); plt.close(fig)
print("wrote figures/fig6b_deviance.png")

# ----------------------------------------------------------------- FIG C ----
ARM_C = "Subclone"
f = RES / f"relvsrnd_{ARM_C}_mincl100_k20.json"
if not f.exists():
    print(f"  [skip fig C] {f.name} not written yet")
else:
    d = json.loads(f.read_text()); FL = d["f_labels"]; xs = np.arange(len(FL))
    cnt = {tg: np.array(d["counts"][tg]) for tg in ("rel", "rnd")}
    hit = {tg: np.array(d["hits"][tg]) for tg in ("rel", "rnd")}
    # ⚠ FACET ALL FOUR STRATA, do not pick one.  The first version auto-picked the
    # best-supported stratum, which is the one where the model ALREADY predicts 77%
    # dropout -- and there the two curves nearly coincide, because a tape the model
    # expects to be missing is missing for everyone and lineage has nothing to add.
    # The effect lives where the model expects the tape to be PRESENT.  That
    # modulation is a result, so the figure shows it rather than cropping to it.
    fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.2), sharex=True, sharey=True)
    for i, ax in enumerate(axs.ravel()):
        pm = 100 * d["p_mid"][i]
        ax.axhline(pm, color=NEUT2, lw=1.0, ls=(0, (4, 3)), zorder=1)
        for tg, lab, col in (("rnd", f"{d['k']} random clone-mates (the null)", NEUT2),
                             ("rel", f"{d['k']} nearest relatives", COL[ARM_C])):
            r = np.where(cnt[tg][i] > 0, hit[tg][i] / np.maximum(cnt[tg][i], 1), np.nan)
            sz = 16 + 90 * (cnt[tg][i] / max(cnt[tg][i].max(), 1)) ** 0.5
            ax.plot(xs, 100 * r, color=col, lw=1.8, zorder=3)
            ax.scatter(xs, 100 * r, s=sz, color=col, zorder=4,
                       edgecolor=SURFACE, linewidth=0.7, label=lab if i == 0 else None)
            xi = len(FL) - 1                      # only "all" -- the thin, load-bearing bin
            if np.isfinite(r[xi]):
                ax.text(xi, 100 * r[xi] + (7 if tg == "rel" else -11),
                        f"n={int(cnt[tg][i][xi]):,}", ha="right", fontsize=6.4, color=col)
        ax.set_title(f"the model predicts {pm:.0f}%", fontsize=9.2, color=INK, pad=6)
        ax.set_ylim(-14, 116)
        ax.set_xticks(xs)
        ax.set_xticklabels(FL, fontsize=8, rotation=30, ha="right")
    for ax in axs[:, 0]:
        ax.set_ylabel("observed dropout rate (%)")
    fig.supxlabel("fraction of those neighbours missing this tape", fontsize=10, y=0.02)
    axs[0, 0].legend(frameon=False, fontsize=8.4, loc="upper left", handlelength=1.2)
    fig.suptitle("Subclone · held-out entries · dashed line = the model's own prediction, held fixed",
                 fontsize=8.6, color=INK2, y=0.975)
    fig.tight_layout(rect=[0, 0.035, 1, 0.955])
    fig.savefig(FIG / "fig6c_relatives.png"); plt.close(fig)
    print("wrote figures/fig6c_relatives.png")
