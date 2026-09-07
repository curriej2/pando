#!/usr/bin/env python3
r"""
================================================================================
 Pool the per-combo nulls and calibrate the MARGIN
================================================================================

Input: results/percombo/percombo_{arm}_p*.npz from 53_percombo.py.

THE STATISTIC, per (clade, tape) combo:

    m_obs = Lambda_obs - max_{b >= HOLD} Lambda_b            [nats]

"how far above the best of its own ~997 permutations the real labelling falls".
⚠ Reported in NATS, not as a rank: with a per-combo null tail below 1e-6 the
rank saturates at 1/(B+1) for every real event and cannot order them, whereas
the margin is uncensored and is exactly what makes the question answerable.

WHY NO CLADE-SIZE STRATA HERE.  40_event_catalogue.py needs six size strata to
stop a 5-cell clade being judged against a null dominated by 500-cell clades.
Each combo here is judged against its OWN permutations, which already condition
on its clade size, its own p_tilde values and its clone's composition -- at full
resolution. The "a 9-cell Pre-TX clade caps at 12.4 nats" problem also
disappears, because the cap applies to that clade's null too: topping out at 6
nats and scoring 11 is decisive whatever any global threshold says.

MULTIPLICITY, calibrated exactly.  53 holds out the first HOLD permutations, so

    m_null = Lambda_{b < HOLD} - max_{b >= HOLD} Lambda_b     [HOLD x NC draws]

is the margin of a candidate over the max of the SAME permutation set -- exactly
exchangeable with m_obs under H0. Hence

    FDR(t) = ( #{m_null >= t} / HOLD ) / #{m_obs >= t}

monotonised BH-style by a running minimum UP the grid (an event at margin t may
be reported at any threshold <= t).

⚠ The margin, not a z-score. For a small clone the null is genuinely coarse -- a
4-cell clade in a 6-cell clone has only C(6,4)=15 distinct compositions, so
sd can be exactly 0 and z explodes. sd is reported as context only.

Usage: 54_percombo_merge.py <arm> [--maxd 4] [--budget 2] [--eps 0.01]
Outputs: results/percombo_{arm}.json          (aggregate, committable)
         results/percombo_{arm}.tsv.gz        (significant combos; gitignored)
================================================================================
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4

arm = sys.argv[1]
MAXD = int(sys.argv[sys.argv.index("--maxd") + 1]) if "--maxd" in sys.argv else 4
BUDGET = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 2.0
tag = "" if MAXD == 4 else f"_d{MAXD}"

parts = sorted((RES / "percombo").glob(f"percombo_{arm}{tag}_p*.npz"))
assert parts, f"no parts for {arm}{tag}"
obs = valid = size = None
n_ge = mx = s1 = s2 = None
n_acc, HOLD, pseudo = 0, None, {}
for f in parts:
    z = np.load(f, allow_pickle=False)
    if obs is None:
        obs, valid, size = z["obs"], z["valid"], z["size"]
        n_ge = z["n_ge"].astype(np.int32).copy()
        mx = z["mx"].copy(); s1 = z["s1"].copy(); s2 = z["s2"].copy()
    else:
        assert np.array_equal(z["obs"], obs), f"{f.name}: observed Lambda differs"
        assert np.array_equal(z["valid"], valid), f"{f.name}: layout differs"
        n_ge += z["n_ge"]; s1 += z["s1"]; s2 += z["s2"]
        np.maximum(mx, z["mx"], out=mx)
    m = z["meta"]; n_acc += int(m[4]); HOLD = int(m[7])
    for b in z["hold_ids"]:
        pseudo[int(b)] = z[f"pseudo{int(b)}"]
    del z
B = n_acc
assert len(pseudo) == HOLD, f"expected {HOLD} held-out permutations, got {sorted(pseudo)}"
print(f"{arm}{tag}: {len(parts)} parts, {B:,} accumulated permutations, "
      f"{HOLD} held out, {obs.size:,} combo slots ({valid.sum():,} scorable)", flush=True)

V = valid & np.isfinite(mx) & (mx > -1e29)
m_obs = np.where(V, obs - mx, -np.inf)
m_null = [np.where(V, pseudo[b] - mx, -np.inf) for b in sorted(pseudo)]
NCV = int(V.sum())

grid = np.unique(np.round(np.concatenate(
    [np.arange(0.0, 5.0, 0.05), np.geomspace(5.0, 2000.0, 300)]), 4))
o_n = np.array([(m_obs >= t).sum() for t in grid], float)
u_n = np.array([np.mean([(v >= t).sum() for v in m_null]) for t in grid])
raw = np.where(o_n > 0, u_n / np.maximum(o_n, 1.0), 1.0)
q = np.minimum(np.minimum.accumulate(raw), 1.0)


def pick(crit, val):
    ok = np.flatnonzero((u_n <= val) if crit == "abs" else (q <= val))
    return int(ok[0]) if ok.size else -1


print(f"\n  margin FDR curve (nats: observed / null / FDR)")
for t in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0):
    j = int(np.searchsorted(grid, t, side="right") - 1)
    print(f"   >={grid[j]:>7.2f}  {o_n[j]:>12,.0f}  {u_n[j]:>12,.1f}  {100*q[j]:>8.4f}%")
rows = {}
for nm, crit, val in (("fdr05", "fdr", 0.05), ("fdr01", "fdr", 0.01),
                      ("budget", "abs", BUDGET)):
    j = pick(crit, val)
    rows[nm] = dict(j=j, thr=(float(grid[j]) if j >= 0 else None),
                    n=(float(o_n[j]) if j >= 0 else 0.0),
                    exp_false=(float(u_n[j]) if j >= 0 else 0.0),
                    fdr=(float(q[j]) if j >= 0 else None))
    r = rows[nm]
    print(f"  {nm:>7}: margin >= {r['thr']}  ->  {r['n']:,.0f} combos, "
          f"{r['exp_false']:,.1f} expected false, FDR "
          f"{100*r['fdr']:.4f}%" if j >= 0 else f"  {nm:>7}: unreachable")

# ---- rebuild the layout so a slot index maps back to (clone, depth, anchor, tape)
z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s_ = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s_], clone[s_]
codes = np.load(RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4
                       else f"prefix_codes_{arm}.npz"), allow_pickle=False)["codes"]
n, K = Y.shape
cl_names, g_clone = np.unique(clone, return_inverse=True)
miss = ~Y
J = rows["budget"]["j"]
THRM = grid[J] if J >= 0 else np.inf
hits, off = [], 0
for d in range(1, MAXD + 1):
    for a_ in range(K):
        cd = codes[:, a_, d - 1]
        ok = cd >= 0
        if ok.sum() < MIN_CLADE:
            continue
        idx = np.flatnonzero(ok)
        key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
        _, sub = np.unique(key, return_inverse=True)
        G = int(sub.max() + 1)
        ns = np.bincount(sub, minlength=G)
        keep = np.flatnonzero(ns >= MIN_CLADE)
        if keep.size == 0:
            continue
        blk = m_obs[off:off + keep.size * K].reshape(keep.size, K)
        gg, zz = np.nonzero(blk >= THRM)
        if gg.size:
            o2 = np.argsort(sub, kind="stable")
            cb = np.concatenate([[0], np.cumsum(ns)])
            owner = np.zeros(G, dtype=np.int64); owner[sub] = g_clone[idx]
            for gi, z_ in zip(gg, zz):
                gk = int(keep[gi]); cells = idx[o2[cb[gk]:cb[gk + 1]]]
                sl = off + gi * K + z_
                hits.append(dict(clone=int(owner[gk]), depth=d, anchor=a_, tape=int(z_),
                                 size=int(ns[gk]), lam=float(obs[sl]),
                                 nullmax=float(mx[sl]), margin=float(m_obs[sl]),
                                 n_ge=int(n_ge[sl]),
                                 mean=float(s1[sl] / max(B, 1)),
                                 sd=float(np.sqrt(max(s2[sl] / max(B, 1)
                                                      - (s1[sl] / max(B, 1)) ** 2, 0.0))),
                                 n_missing=int(miss[cells, z_].sum()), cells=cells))
        off += keep.size * K
assert off == obs.size, (off, obs.size)

hits.sort(key=lambda h: -h["margin"])
kept, seen = [], {}
for h in hits:
    k = (h["clone"], h["tape"]); ss = set(h["cells"].tolist())
    if any(ss & prev for prev in seen.get(k, [])):
        continue
    seen.setdefault(k, []).append(ss); kept.append(h)
cxt = sum(h["n_missing"] for h in kept)
tot = int(miss.sum())
print(f"\n  at the budget margin ({THRM:.2f} nats): {len(hits):,} combos -> "
      f"{len(kept):,} events after the overlap collapse")
print(f"  cell x tape entries inside them: {cxt:,} of {tot:,} = {100*cxt/max(tot,1):.2f}%")
if kept:
    mg = np.array([h["margin"] for h in kept]); sz = np.array([h["size"] for h in kept])
    print(f"  margin: median {np.median(mg):.1f}, max {mg.max():.0f} nats | "
          f"clade size median {np.median(sz):.0f}, min {sz.min()}")
    print(f"  combos whose margin beats EVERY permutation (n_ge = 0): "
          f"{sum(1 for h in kept if h['n_ge'] == 0):,} of {len(kept):,}")

out = {"arm": arm, "max_depth": MAXD, "B": B, "hold_out": HOLD,
       "n_combo_slots": int(obs.size), "n_scorable": NCV,
       "criteria": rows, "budget": BUDGET,
       "n_events": len(kept), "n_combos_above": len(hits),
       "cellxtape_in_events": cxt, "total_missing_entries": tot,
       "frac_missing_in_events": cxt / max(tot, 1),
       "margin_curve": {"nats": grid.tolist(), "observed": o_n.tolist(),
                        "null": u_n.tolist(), "fdr": q.tolist()}}
(RES / f"percombo_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"percombo_{arm}{tag}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\tdepth\tanchor\ttape\tclade_cells\tn_missing\tlambda_nats\t"
             "null_max\tmargin_nats\tn_ge\tnull_mean\tnull_sd\n")
    for h in kept:
        fh.write(f"{h['clone']}\t{cl_names[h['clone']]}\t{h['depth']}\t{h['anchor']}\t"
                 f"{h['tape']}\t{h['size']}\t{h['n_missing']}\t{h['lam']:.3f}\t"
                 f"{h['nullmax']:.3f}\t{h['margin']:.3f}\t{h['n_ge']}\t"
                 f"{h['mean']:.3f}\t{h['sd']:.3f}\n")
print(f"\nwrote results/percombo_{arm}{tag}.json and .tsv.gz")
