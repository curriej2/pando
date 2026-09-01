#!/usr/bin/env python3
r"""
FIGURE 3, PANEL c -- the shelf, and the admission threshold that hides it.

Panel a showed Mouse1's R_c is bimodal: a mode near 120 plus a flat SHELF running
from the admission cut at 20 up to ~100.  This panel puts all five arms on that
axis and marks the two thresholds the paper actually used -- >=20 recovered tapes
for Mouse1-3, >=100 for Initial and Subclone.

CLAIM.  The shelf is not absent from Pre-TX and Subclone; it was DELETED from them
by QC and KEPT in the mice.  Every apparent arm difference in rho_cell is that
choice: cut all five arms at >=100 and rho_cell collapses to 0.012-0.024
everywhere, with the mice the most homogeneous of the set.

WHAT THIS DOES NOT SAY.  The mouse analysis population IS admitted at >=20, so
rho_cell = 0.13-0.18 is the real number any inference on those data must handle.
The >=100 comparison establishes what KIND of thing it is -- a data-quality
population with no biological content, hence marginalisable -- as against the tape
axis of panel b, which is informative and must be modelled explicitly.

DIAGNOSTIC RUN FIRST (printed, and it could have changed the claim).  If the shelf
concentrated in particular harvest sites it would be a dissection/prep batch
effect, not per-cell capture.  Tested with the same variance machinery used for
rho throughout: for the indicator Z_c = 1[20 <= R_c < 100], rho_sample is the
fraction of Var(Z) explained by which sample the cell came from,

    rho_sample = [Var_s(mean Z per sample) - noise] / [s(1-s)],  noise = s(1-s)/n_s

pooled and per arm.  Small rho_sample => the shelf is a per-cell property.

Output: figures/fig3c_shelf.png, results/fig3c_shelf.json
"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
ARMS = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]
SERIES = {"Mouse1": "#2a78d6", "Mouse2": "#1baf7a", "Mouse3": "#eda100",
          "Initial": "#e87ba4", "Subclone": "#4a3aa7"}
LABEL = {"Initial": "Pre-TX", "Subclone": "Subclone", "Mouse1": "Mouse 1",
         "Mouse2": "Mouse 2", "Mouse3": "Mouse 3"}
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
LO, HI, BINW = 20, 100, 4

R, SM, SN = {}, {}, {}
for a in ARMS:
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y = z["recovered"]
    R[a] = Y.sum(1); SM[a] = z["sample"]; SN[a] = [str(x) for x in z["sample_names"]]

def rho_group(Z, g):
    """Fraction of Var(Z) explained by the grouping g, noise-corrected."""
    s = float(Z.mean()); tot = s * (1 - s)
    ids = np.unique(g)
    means = np.array([Z[g == i].mean() for i in ids], dtype=float)
    ns = np.array([(g == i).sum() for i in ids], dtype=float)
    w = ns / ns.sum()
    between = float(((means - s) ** 2 * w).sum())
    noise = float((tot / ns * w).sum())          # sampling var of a group mean
    return (between - noise) / tot, len(ids), s

out = {"shelf_range": [LO, HI], "arms": {}}
print("=== diagnostic: is the shelf a batch (harvest-site) effect? ===")
print("    rho_sample = fraction of Var(1[cell in shelf]) explained by sample id\n")
for a in ARMS:
    Z = ((R[a] >= LO) & (R[a] < HI)).astype(float)
    rs, k, s = rho_group(Z, SM[a])
    per = {SN[a][i]: {"n": int((SM[a] == i).sum()),
                      "mean_R": float(R[a][SM[a] == i].mean()),
                      "shelf_frac": float(Z[SM[a] == i].mean())}
           for i in range(len(SN[a])) if (SM[a] == i).sum() >= 100}
    out["arms"][a] = {"n": int(R[a].size), "shelf_frac": s,
                      "frac_below_thr": float((R[a] < THR[a]).mean()),
                      "rho_sample": float(rs), "n_samples": k, "by_sample": per}
    print(f"{a:9s} shelf {100*s:5.1f}%  rho_sample {rs:+.4f} over {k} samples")
    for nm, d in per.items():
        print(f"            {nm:12s} n={d['n']:6,d}  mean R_c {d['mean_R']:6.1f}  "
              f"shelf {100*d['shelf_frac']:5.1f}%")

Zp = np.concatenate([((R[a] >= LO) & (R[a] < HI)).astype(float) for a in ARMS[:3]])
gp = np.concatenate([SM[a] + 10 * i for i, a in enumerate(ARMS[:3])])
rs_pool, k_pool, s_pool = rho_group(Zp, gp)
out["pooled_mice"] = {"rho_sample": float(rs_pool), "n_samples": k_pool,
                      "shelf_frac": float(s_pool)}
print(f"\npooled mice: shelf {100*s_pool:.1f}%, rho_sample {rs_pool:+.4f} "
      f"over {k_pool} harvest samples")

marg = json.loads((RES / "dropout_margins.json").read_text())
fol = json.loads((RES / "dropout_followup.json").read_text())
rho_raw = {a: marg[a]["all"]["cell"]["rho"] for a in ARMS}
rho_cut = {a: fol["shelf"][a]["cut100"]["rho_cell"] for a in ARMS}
out["rho_cell_raw"], out["rho_cell_cut100"] = rho_raw, rho_cut
print("\nrho_cell   raw -> at a common >=100 cut")
for a in ARMS:
    print(f"  {LABEL[a]:9s} {rho_raw[a]:.4f} -> {rho_cut[a]:.4f}")
(RES / "fig3c_shelf.json").write_text(json.dumps(out, indent=1))

# ------------------------------------------------------------------- the panel
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SHADE = "#eceae1"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": INK2,
    "ytick.color": INK2, "font.size": 12, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "axes.titlesize": 14})
fig, ax = plt.subplots(figsize=(9.0, 5.6))

edges = np.arange(0, 168 + BINW, BINW)
ctr = edges[:-1] + BINW / 2
ax.axvspan(LO, HI, color=SHADE, lw=0, zorder=1)
top = 0.0
for a in ARMS:
    h = np.histogram(R[a], bins=edges)[0] / R[a].size
    top = max(top, h.max())
    ax.plot(ctr, h, color=SERIES[a], lw=2.2, zorder=4,
            label=f"{LABEL[a]}      {100*out['arms'][a]['shelf_frac']:.0f}%")
ax.set_ylim(0, top * 1.30)

for x, txt in ((LO, "mice admitted\nat $\\geq$20 tapes"),
               (HI, "Pre-TX and Subclone\nadmitted at $\\geq$100")):
    ax.axvline(x, color=INK2, lw=1.3, ls=(0, (4, 3)), zorder=3)
    ax.text(x + 2.5, top * 1.26, txt, color=INK2, fontsize=10.5, va="top", ha="left")

leg = ax.legend(frameon=True, facecolor="white", edgecolor=AXIS, framealpha=1.0,
                fontsize=11, loc="upper left",
                bbox_to_anchor=(0.015, 0.72), handlelength=1.7, labelspacing=0.42,
                title="share of cells in the shaded band",
                title_fontsize=10.5, alignment="left", borderpad=0.6)
leg.get_frame().set_linewidth(0.8)
leg.get_title().set_color(INK2)
for t in leg.get_texts():
    t.set_color(INK)

ax.set_xlim(0, 167)
ax.set_xlabel("tapes recovered per cell   $R_c$   (of 166)")
ax.set_ylabel(f"fraction of cells   ({BINW}-tape bins)")
ax.set_title("c   The shelf is a QC choice, not a difference between arms",
             loc="left", pad=30)
ax.text(0, 1.015, "cut all five arms at $\\geq$100 tapes and $\\rho_{cell}$ falls "
        "0.13–0.18  $\\rightarrow$  0.012–0.024", transform=ax.transAxes,
        ha="left", va="bottom", fontsize=11.5, color=INK2)

fig.tight_layout()
fig.savefig(FIG / "fig3c_shelf.png", dpi=200, facecolor=SURF)
print("\nwrote figures/fig3c_shelf.png")
