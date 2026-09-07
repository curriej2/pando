#!/usr/bin/env python3
r"""
================================================================================
 61 -- THE VALIDATION LADDER for the recovered z -> s(z) map
================================================================================

The map itself is never the evidence -- a spurious map has no reason to be
consistent, injective, or reproducible.  Its STRUCTURE is.  All 166 TargetBCs are
identical across the five arms (verified 2026-09-03: same engineered line, same
integrations), so tape z must name the SAME symbol in Pre-TX, Mouse 1, 2, 3 and
Subclone, recovered from disjoint sets of cells.  That is what this script tests.

 0. ⚠⚠ WITHIN-ARM CLONE SPLIT -- and this, not cross-arm, is the clean test.
    Measured 2026-09-07: the arms SHARE CLONES.  285 of Mouse 1's 296 clones
    (96.3%) also appear in Initial, because the mice were transplanted from the
    pre-TX pool; Mouse1 x Mouse2 share 42.7%, Mouse2 x Mouse3 31.3%.  Silencing
    is heritable and predates transplantation, so a clone measured in two arms
    carries THE SAME losses -- cross-arm agreement is the same events re-measured,
    not a replication.  Two disjoint clone splits of ONE arm share no cell, no
    clone and no lineage, so agreement between them tests the only thing that
    matters: whether z -> s(z) is a property of the LINE.  Only Subclone is
    near-disjoint from the mice (1, 4 and 1 shared clones), so Subclone x mouse is
    the one honest cross-arm pair.

 1. CROSS-ARM CONCORDANCE -- confounded by clone sharing, reported with it.  Agreement between two independently
    estimated maps, against a permutation null that permutes tape labels within
    one arm.  ⚠ The null must be the PERMUTATION, not 1/|alphabet|: both maps
    have pile-ups (Mouse 2 sends 15 tapes to AAGCGGA, an artefact of its single
    dominant clone making the design matrix near rank-1), and a naive 1/98 would
    credit that as signal.
 2. CROSS-METHOD -- 60 (deconvolution over clades and clones) against 58 (the
    cell-level scan).  The two fail differently: 60 by collinearity of the loss
    sets, 58 by collapse bias and cell-capture structure.  Agreement between them
    is therefore worth more than either alone.
 3. INJECTIVITY.  166 integrations drawing NNNN uniformly from 4^4 = 256 predicts
    256(1 - (1 - 1/256)^166) = 122 distinct symbols; restricted to the n tapes an
    arm can actually call, and to the A symbols conforming there, the prediction
    is A(1 - (1 - 1/A)^n).  A map far below that is piling up, i.e. unresolved.
 4. DOSE.  A symbol named by j integrations should carry about j x the base xi.
 5. corr(beta_z, xi_s(z)) > 0 -- tape recovery rate and symbol frequency are two
    readouts of one locus's expression level.

Usage: 61_concordance.py [--nperm 2000] [--zmax -3] [--gapmin 0.5]
Output: results/concordance.json
================================================================================
"""
import json, glob, itertools, collections, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_w = __import__("57_writes")

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 2000
ZMAX = float(sys.argv[sys.argv.index("--zmax") + 1]) if "--zmax" in sys.argv else -3.0
GAPMIN = float(sys.argv[sys.argv.index("--gapmin") + 1]) if "--gapmin" in sys.argv else 0.5
rng = np.random.default_rng(20260907)

# ---- load the two independent maps ----------------------------------------
dec, cell, meta = {}, {}, {}
for f in sorted(glob.glob(str(RES / "depmap_*_both.json"))):
    d = json.load(open(f)); a = d["arm"]
    dec[a] = {m["tape"]: (m["symbol"], m["z"], m["gap"]) for m in d["map"] if m["callable"]}
    meta[a] = {"n_units": d["n_units"], "rank": d["design"]["rank"],
               "n_symbols": d["n_symbols"], "vif_max": d["design"]["vif_max"]}
for f in sorted(glob.glob(str(RES / "tape_cis_Z_*.npz"))):
    a = Path(f).stem.replace("tape_cis_Z_", "")
    z = np.load(f, allow_pickle=False)
    Z, syms, tps, Nt = z["Z"], z["symbols"], z["tapes"], z["Nt"]
    # ⚠ the cell-level scan looks for ENRICHMENT of s(t) at tape t (a cell that
    # kept tape t also kept its pegRNA), the opposite sign to 60's depletion.
    bi = np.argmax(Z, axis=1); bz = Z[np.arange(len(tps)), bi]
    rest = Z.copy(); rest[np.arange(len(tps)), bi] = -np.inf
    gap = bz - np.max(rest, axis=1)
    cell[a] = {str(tps[t]): (str(syms[bi[t]]), float(bz[t]), float(gap[t]))
               for t in range(len(tps)) if Nt[t] >= 100}


def agree(m1, m2, tapes, perm=None):
    """Agreement count; with `perm`, m2's tape labels are permuted."""
    if perm is None:
        return sum(m1[t][0] == m2[t][0] for t in tapes)
    return sum(m1[t][0] == m2[p][0] for t, p in zip(tapes, perm))


def pair(m1, m2, label, filt=None):
    common = sorted(set(m1) & set(m2))
    if filt:
        common = [t for t in common if filt(m1[t]) and filt(m2[t])]
    n = len(common)
    if n < 5:
        return {"pair": label, "n": n, "note": "too few"}
    obs = agree(m1, m2, common)
    null = np.array([agree(m1, m2, common, rng.permutation(common)) for _ in range(NPERM)])
    p = (1 + int((null >= obs).sum())) / (NPERM + 1)
    return {"pair": label, "n": n, "agree": obs, "rate": obs / n,
            "null_mean": float(null.mean()), "null_p95": float(np.percentile(null, 95)),
            "null_max": int(null.max()), "enrichment": obs / max(null.mean(), 1e-9), "p": p}


out = {"nperm": NPERM, "filters": {"zmax": ZMAX, "gapmin": GAPMIN}, "meta": meta}

# ---- clone sharing, so every cross-arm row can be read against it ---------
import numpy as _np
share = {}
_used = {}
for a in list(dec):
    try:
        z = _np.load(RES / f"dropout_matrix_{a}.npz", allow_pickle=False)
        cn, cl = z["clone_names"], z["clone"]
        _used[a] = set(str(cn[c]) for c in _np.unique(cl[cl >= 0]))
    except FileNotFoundError:
        pass
for a, b in itertools.combinations(sorted(_used), 2):
    s = _used[a] & _used[b]
    share[f"{a} x {b}"] = len(s) / max(min(len(_used[a]), len(_used[b])), 1)
out["clone_sharing"] = share

print("=== 0. WITHIN-ARM CLONE SPLIT (no shared cell, clone or lineage) ===")
out["split"] = []
for f in sorted(glob.glob(str(RES / "depmap_*_both_s0of2.json"))):
    a = json.load(open(f))["arm"]
    g = RES / f"depmap_{a}_both_s1of2.json"
    if not g.exists():
        continue
    m0 = {m["tape"]: (m["symbol"], m["z"], m["gap"])
          for m in json.load(open(f))["map"] if m["callable"]}
    m1 = {m["tape"]: (m["symbol"], m["z"], m["gap"])
          for m in json.load(open(g))["map"] if m["callable"]}
    for lab, filt in (("all", None),
                      ("strict", lambda v: v[1] <= ZMAX and v[2] >= GAPMIN)):
        r = pair(m0, m1, f"{a} half0 x half1 [{lab}]", filt=filt)
        out["split"].append(r)
        if "note" in r:
            print(f"{r['pair']:30s} {r['n']:5d}   (too few)"); continue
        print(f"{r['pair']:30s} {r['n']:5d} {r['agree']:6d} {r['rate']:7.3f} "
              f"{r['null_mean']:8.2f} {r['null_max']:8d} {r['enrichment']:6.2f} {r['p']:8.4f}")

print("\n=== 1. CROSS-ARM CONCORDANCE, deconvolution map (60) "
      "-- READ WITH THE CLONE-SHARING COLUMN ===")
print(f"{'pair':22s} {'n':>5} {'agree':>6} {'rate':>7} {'null mu':>8} {'nullmax':>8} {'x':>6} {'p':>8}")
out["cross_arm_dec"] = []
for a, b in itertools.combinations(sorted(dec), 2):
    r = pair(dec[a], dec[b], f"{a} x {b}")
    out["cross_arm_dec"].append(r)
    if "note" in r:
        print(f"{r['pair']:22s} {r['n']:5d}   (too few)"); continue
    print(f"{r['pair']:22s} {r['n']:5d} {r['agree']:6d} {r['rate']:7.3f} "
          f"{r['null_mean']:8.2f} {r['null_max']:8d} {r['enrichment']:6.2f} {r['p']:8.4f}"
          f"   shared clones {100*share.get(r['pair'], float('nan')):5.1f}%")

print("\n=== 1b. restricted to confident calls (z <= %.1f, gap >= %.1f) ===" % (ZMAX, GAPMIN))
out["cross_arm_dec_strict"] = []
for a, b in itertools.combinations(sorted(dec), 2):
    r = pair(dec[a], dec[b], f"{a} x {b}", filt=lambda v: v[1] <= ZMAX and v[2] >= GAPMIN)
    out["cross_arm_dec_strict"].append(r)
    if "note" in r:
        print(f"{r['pair']:22s} {r['n']:5d}   (too few)"); continue
    print(f"{r['pair']:22s} {r['n']:5d} {r['agree']:6d} {r['rate']:7.3f} "
          f"{r['null_mean']:8.2f} {r['null_max']:8d} {r['enrichment']:6.2f} {r['p']:8.4f}")

print("\n=== 2. CROSS-ARM CONCORDANCE, cell-level scan (58) ===")
out["cross_arm_cell"] = []
for a, b in itertools.combinations(sorted(cell), 2):
    r = pair(cell[a], cell[b], f"{a} x {b}")
    out["cross_arm_cell"].append(r)
    if "note" in r:
        continue
    print(f"{r['pair']:22s} {r['n']:5d} {r['agree']:6d} {r['rate']:7.3f} "
          f"{r['null_mean']:8.2f} {r['null_max']:8d} {r['enrichment']:6.2f} {r['p']:8.4f}")

print("\n=== 3. CROSS-METHOD, within arm: 60 vs 58 ===")
out["cross_method"] = []
for a in sorted(set(dec) & set(cell)):
    r = pair(dec[a], cell[a], f"{a}: dec x cell")
    out["cross_method"].append(r)
    if "note" in r:
        print(f"{r['pair']:22s} {r['n']:5d}   (too few)"); continue
    print(f"{r['pair']:22s} {r['n']:5d} {r['agree']:6d} {r['rate']:7.3f} "
          f"{r['null_mean']:8.2f} {r['null_max']:8d} {r['enrichment']:6.2f} {r['p']:8.4f}")

print("\n=== 4. INJECTIVITY: distinct symbols called vs the collision prediction ===")
out["injectivity"] = []
for a in sorted(dec):
    n = len(dec[a]); A = meta[a]["n_symbols"]
    d = len(set(v[0] for v in dec[a].values()))
    exp = A * (1 - (1 - 1.0 / A) ** n)
    out["injectivity"].append({"arm": a, "n_called": n, "distinct": d,
                               "expected_if_injective": exp, "ratio": d / exp})
    print(f"  {a:9s} {n:3d} tapes -> {d:3d} distinct   (injective draws predict "
          f"{exp:5.1f})   ratio {d/exp:.2f}")

(RES / "concordance.json").write_text(json.dumps(out, indent=1))
print("\n-> results/concordance.json")
