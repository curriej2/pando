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

⚠⚠ CORRECTION (2026-09-07, measured).  The first version of this script applied
ONE GLOBAL margin threshold and claimed the per-combo null "dissolves" the
clade-size problem.  HALF RIGHT, and the half that was wrong matters: the NULL
conditions on clade size, but a GLOBAL THRESHOLD ON THE MARGIN does not.  A
5-cell clade's Lambda range is compressed, so it cannot reach a 36-nat margin
however complete the loss.  Measured on Mouse1: clades of 4-9 cells are 55% of
scorable combos and 0.0% of called events, while 200+ cell clades are 1.3% of
combos and 36.7% of events.

So the margin is calibrated PER CLADE-SIZE STRATUM here.  That is the synthesis,
not a retreat: the per-combo null conditions on the individual clade (its size,
its own p_tilde values, its clone's composition) which no six-bin scheme can do,
and per-stratum calibration then makes the THRESHOLD comparable across sizes,
which a single global margin cut cannot.  Each piece does what it is good at.
The "a 9-cell Pre-TX clade caps at 12.4 nats" problem is genuinely gone, because
the cap now applies to that clade's null AND to its stratum's threshold.

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

SIZE_EDGES = np.array([4, 6, 10, 20, 50, 200, 10 ** 9])
NB = len(SIZE_EDGES) - 1
SIZE_LABEL = ["4-5", "6-9", "10-19", "20-49", "50-199", "200+"]
strat = np.clip(np.searchsorted(SIZE_EDGES, size, side="right") - 1, 0, NB - 1)

grid = np.unique(np.round(np.concatenate(
    [np.arange(0.0, 5.0, 0.05), np.geomspace(5.0, 2000.0, 300)]), 4))
o_all = np.array([(m_obs >= t).sum() for t in grid], float)
u_all = np.array([np.mean([(v >= t).sum() for v in m_null]) for t in grid])
q_all = np.minimum(np.minimum.accumulate(
    np.where(o_all > 0, u_all / np.maximum(o_all, 1.0), 1.0)), 1.0)
print(f"\n  margin FDR curve, POOLED (nats: observed / null / FDR)")
for t in (0.0, 0.5, 2.0, 5.0, 25.0):
    j = int(np.searchsorted(grid, t, side="right") - 1)
    print(f"   >={grid[j]:>7.2f}  {o_all[j]:>13,.0f}  {u_all[j]:>13,.1f}  {100*q_all[j]:>8.4f}%")

# ---- per-stratum calibration of the margin
THRS = np.full(NB, np.inf)
QS, ON, UN, js = [], [], [], np.full(NB, -1, dtype=int)
print(f"\n  per-stratum margin threshold at expected false <= {BUDGET:g} per stratum:")
print(f"   {'clade':>8} {'combos':>12} {'thr(nats)':>10} {'combos>=':>10} "
      f"{'exp.false':>10} {'FDR%':>9}")
for i in range(NB):
    m = V & (strat == i)
    o_i = np.array([(m_obs[m] >= t).sum() for t in grid], float)
    u_i = np.array([np.mean([(v[m] >= t).sum() for v in m_null]) for t in grid])
    q_i = np.minimum(np.minimum.accumulate(
        np.where(o_i > 0, u_i / np.maximum(o_i, 1.0), 1.0)), 1.0)
    ON.append(o_i); UN.append(u_i); QS.append(q_i)
    ok = np.flatnonzero(u_i <= BUDGET)
    if ok.size and o_i[int(ok[0])] > 0:
        js[i] = int(ok[0]); THRS[i] = grid[js[i]]
        print(f"   {SIZE_LABEL[i]:>8} {int(m.sum()):>12,} {THRS[i]:>10.2f} "
              f"{o_i[js[i]]:>10,.0f} {u_i[js[i]]:>10.2f} {100*q_i[js[i]]:>8.4f}%")
    else:
        print(f"   {SIZE_LABEL[i]:>8} {int(m.sum()):>12,} {'NONE':>10} "
              f"{'--':>10} {'--':>10} {'--':>9}")
QS, ON, UN = np.array(QS), np.array(ON), np.array(UN)
rows = {"budget": dict(
    per_stratum={SIZE_LABEL[i]: (float(THRS[i]) if js[i] >= 0 else None)
                 for i in range(NB)},
    n=float(sum(ON[i][js[i]] for i in range(NB) if js[i] >= 0)),
    exp_false=float(sum(UN[i][js[i]] for i in range(NB) if js[i] >= 0)))}
print(f"  => {rows['budget']['n']:,.0f} combos above their stratum's margin threshold, "
      f"{rows['budget']['exp_false']:.1f} expected false in total")

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
THRM_BY_SLOT = THRS[strat]        # each combo judged at its own stratum's threshold
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
        sl0 = slice(off, off + keep.size * K)
        blk = m_obs[sl0].reshape(keep.size, K)
        thrblk = THRM_BY_SLOT[sl0].reshape(keep.size, K)
        gg, zz = np.nonzero(blk >= thrblk)
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
                                 n_ge=int(n_ge[sl]), strat=SIZE_LABEL[strat[sl]],
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
print(f"\n  at the per-stratum margins: {len(hits):,} combos -> "
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
       "size_labels": SIZE_LABEL, "size_edges": SIZE_EDGES.tolist(),
       "strat_margin_thresholds": {SIZE_LABEL[i]: (float(THRS[i]) if js[i] >= 0 else None)
                                   for i in range(NB)},
       "strat_combos": {SIZE_LABEL[i]: int((V & (strat == i)).sum()) for i in range(NB)},
       "pooled_margin_curve": {"nats": grid.tolist(), "observed": o_all.tolist(),
                               "null": u_all.tolist(), "fdr": q_all.tolist()},
       "n_events": len(kept), "n_combos_above": len(hits),
       "cellxtape_in_events": cxt, "total_missing_entries": tot,
       "frac_missing_in_events": cxt / max(tot, 1),
       "strat_margin_curves": {SIZE_LABEL[i]: {"observed": ON[i].tolist(),
                                               "null": UN[i].tolist(),
                                               "fdr": QS[i].tolist()} for i in range(NB)}}
(RES / f"percombo_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"percombo_{arm}{tag}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\tdepth\tanchor\ttape\tclade_cells\tn_missing\tlambda_nats\t"
             "null_max\tmargin_nats\tn_ge\tsize_stratum\tnull_mean\tnull_sd\n")
    for h in kept:
        fh.write(f"{h['clone']}\t{cl_names[h['clone']]}\t{h['depth']}\t{h['anchor']}\t"
                 f"{h['tape']}\t{h['size']}\t{h['n_missing']}\t{h['lam']:.3f}\t"
                 f"{h['nullmax']:.3f}\t{h['margin']:.3f}\t{h['n_ge']}\t"
                 f"{h['strat']}\t{h['mean']:.3f}\t{h['sd']:.3f}\n")
print(f"\nwrote results/percombo_{arm}{tag}.json and .tsv.gz")
