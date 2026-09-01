#!/usr/bin/env python3
r"""
Follow-ups to 24, driven by what the R_c distribution actually looked like.

1. THE SHELF.  R_c is not merely overdispersed, it is bimodal in the mice: a mode
   near 120 plus a flat shelf from ~20 to ~100.  Initial/Subclone were admitted at
   >=100 tapes, which removes exactly that shelf -- so the per-arm difference in
   rho_cell may be the QC threshold rather than biology.  Re-measure every arm at
   a COMMON >=100 cut and see what survives.

2. IS THE SHELF THE SAME CELLS THAT LOST THEIR CLONAL BARCODE?  44.9% of Mouse1
   cells carry no ClonalBC.  If low R_c and barcode loss are one phenomenon
   (per-cell capture quality) then the analysis population is selected on capture,
   and P(no barcode | R_c) should fall steeply with R_c.

3. IS THE TAPE AXIS REPRODUCIBLE?  rho_tape ~ 0.20-0.26 in all five arms.  If the
   per-tape recovery rate is an intrinsic property of the tape it must agree
   across five independently prepared libraries; if it is stochastic it will not.
   Cross-arm correlation of the 166 per-tape rates, plus the within-arm split by
   replicate/harvest sample.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ARMS = ["Initial", "Subclone", "Mouse1", "Mouse2", "Mouse3"]
Y, CL, SM, SN, TP = {}, {}, {}, {}, {}
for a in ARMS:
    z = np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
    Y[a], CL[a], SM[a] = z["recovered"], z["clone"], z["sample"]
    SN[a] = [str(s) for s in z["sample_names"]]
    TP[a] = [str(t) for t in z["tapes"]]

assert all(TP[a] == TP[ARMS[0]] for a in ARMS), "tape order differs between arms"
out = {}

def rho(R, m, p):
    v, null = R.var(ddof=1), m * p * (1 - p)
    return (v / null - 1) / (m - 1)

print("=== 1. the shelf, and a common QC cut ===")
out["shelf"] = {}
for a in ARMS:
    R = Y[a].sum(1); k = Y[a].shape[1]
    row = {"n": int(R.size), "frac_below_100": float((R < 100).mean()),
           "frac_below_60": float((R < 60).mean())}
    for cut in (0, 20, 100):
        s = R >= cut
        if s.sum() < 50: continue
        p = float(Y[a][s].mean())
        row[f"cut{cut}"] = {
            "n": int(s.sum()), "p": p, "sd": float(R[s].std(ddof=1)),
            "rho_cell": float(rho(R[s].astype(float), k, p)),
            "rho_tape": float(rho(Y[a][s].sum(0).astype(float), int(s.sum()), p)),
            "sigma_cell": float(R[s].std(ddof=1) / k),
            "sigma_tape": float((Y[a][s].mean(0)).std(ddof=1))}
    out["shelf"][a] = row
    print(f"{a:9s} <100 tapes: {100*row['frac_below_100']:5.1f}%  <60: "
          f"{100*row['frac_below_60']:5.1f}%  |  rho_cell  raw {row['cut0']['rho_cell']:.4f}"
          f"  at>=100 {row.get('cut100',{}).get('rho_cell',float('nan')):.4f}"
          f"  |  sigma_tape {row['cut0']['sigma_tape']:.3f}")

print("\n=== 2. clonal-barcode loss vs R_c ===")
out["barcode"] = {}
for a in ARMS:
    R = Y[a].sum(1); nobc = CL[a] < 0
    edges = [0, 20, 40, 60, 80, 100, 120, 140, 167]
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (R >= lo) & (R < hi)
        prof.append({"lo": lo, "hi": hi, "n": int(s.sum()),
                     "p_nobc": float(nobc[s].mean()) if s.sum() else None})
    out["barcode"][a] = {"overall_nobc": float(nobc.mean()),
                         "mean_R_with_bc": float(R[~nobc].mean()),
                         "mean_R_without_bc": float(R[nobc].mean()), "profile": prof}
    print(f"{a:9s} no-BC {100*nobc.mean():5.1f}%  mean R_c  with BC {R[~nobc].mean():6.1f}"
          f"  without {R[nobc].mean():6.1f}   P(no BC | R_c): "
          + " ".join(f"{b['lo']}-{b['hi']}:{100*b['p_nobc']:.0f}%" for b in prof if b["n"] > 30))

print("\n=== 3. is the per-tape rate reproducible? ===")
rates = {a: Y[a].mean(0) for a in ARMS}
M = np.array([rates[a] for a in ARMS])
C = np.corrcoef(M)
out["tape_rate_corr"] = {"arms": ARMS, "pearson": C.tolist(),
                         "rates": {a: rates[a].tolist() for a in ARMS}}
print("        " + "".join(f"{a:>10s}" for a in ARMS))
for i, a in enumerate(ARMS):
    print(f"{a:9s}" + "".join(f"{C[i,j]:10.3f}" for j in range(len(ARMS))))
print("per-tape rate range (min..max across the 166 tapes):")
for a in ARMS:
    r = rates[a]
    print(f"  {a:9s} {r.min():.3f} .. {r.max():.3f}   sd {r.std(ddof=1):.3f}"
          f"   decile 1 {np.quantile(r,0.1):.3f}  decile 9 {np.quantile(r,0.9):.3f}")

print("\nwithin-arm, by sample (independent libraries where replicates exist):")
out["by_sample"] = {}
for a in ARMS:
    d = {}
    for si, s in enumerate(SN[a]):
        m = SM[a] == si
        if m.sum() < 100: continue
        d[s] = {"n": int(m.sum()), "p": float(Y[a][m].mean()),
                "mean_R": float(Y[a][m].sum(1).mean()),
                "rate": Y[a][m].mean(0).tolist()}
    keys = list(d)
    cc = np.corrcoef(np.array([d[s]["rate"] for s in keys])) if len(keys) > 1 else None
    out["by_sample"][a] = {"samples": {s: {k: v for k, v in d[s].items() if k != "rate"}
                                       for s in keys},
                           "tape_rate_corr": cc.tolist() if cc is not None else None}
    print(f"  {a:9s} " + "  ".join(f"{s}: n={d[s]['n']:,} meanR={d[s]['mean_R']:.0f}" for s in keys))
    if cc is not None:
        print(f"            per-tape rate corr between samples: "
              f"{np.min(cc[np.triu_indices(len(keys),1)]):.3f}"
              f"..{np.max(cc[np.triu_indices(len(keys),1)]):.3f}")

(RES / "dropout_followup.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/dropout_followup.json")
