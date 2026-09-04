#!/usr/bin/env python3
r"""
================================================================================
 An event catalogue: candidate heritable dropout events, one row each
================================================================================

Turns "rho is elevated" into a list of concrete (clade, tape) losses, each
individually inspectable, plus a decomposition saying how much of all missingness
they account for.

MECHANISM THIS IS TESTING (verified from the paper, 2026-09-03).  Park's recorder
is one piggyBac cassette, `PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC`, and the
tape is read from cDNA -- 10x 3' v4 with "Feature cDNA Primer 3 ... so that the
Read 2N primer was included and TAPE cDNA co-amplified".  So a tape is recovered
only if ITS INTEGRATION IS TRANSCRIBED.  Epigenetic silencing of an integration
therefore (i) removes the tape from the readout and (ii) removes the co-integrated
pegRNA's symbol from the cell's writing pool, both heritably.  Cells were sorted
mRFP+, but mRFP comes from all 166 integrations, so silencing one leaves the cell
in the data -- no selection confound.

--------------------------------------------------------------------------------
 THE STATISTIC: a log-likelihood ratio for a Dollo loss on the clade's stem
--------------------------------------------------------------------------------
For clade S and tape z, two stories for what we see in each cell:

  H1 (loss)   the integration was silenced on S's stem, so P(missing) = 1 - eps
              for every cell in S.  eps is leakage (ambient RNA, mis-assignment,
              partial silencing).
  H0 (none)   each cell independently, P(missing) = p_hat_cz = logit^-1(alpha_c +
              beta_z) -- which already knows how good the cell and tape are.

    Lambda(S,z) = SUM_{c in S, missing} log[(1-eps) / p_hat_cz]
                + SUM_{c in S, present} log[eps / (1 - p_hat_cz)]

Two properties, and they are exactly what we need:

 1. CAPTURE-INDEPENDENCE IS BUILT IN.  A missing cell earns evidence in proportion
    to how surprising its missingness is, and "surprising" means precisely "this
    cell recovered plenty of other tapes and this tape is normally reliable".
    p_hat=0.05 missing -> +2.98 nats; p_hat=0.90 missing -> +0.09.  The score is
    dominated by cells transcript sampling cannot explain.
 2. IT DEMANDS ALL-OR-NONE.  One present cell costs up to -4.55 nats, cancelling
    ~1.5 strongly-surprising missing cells.  A clade missing the tape in 70% of
    cells scores poorly however large it is.  That is what separates a discrete
    loss from an elevated rate.  Note the penalty is LARGER for a good cell
    (-4.55 at p_hat=0.05) than a bad one (-2.30 at 0.90): a well-captured cell
    would have shown the tape if it were there, so seeing it is near-decisive.

Lambda is in nats, and is the log-likelihood a per-tape absorbing state would
gain from this clade -- so the catalogue total answers "what does modelling this
buy?" in the units the model uses.  It is NOT self-calibrating (eps is fixed, the
clade was chosen by search), so the threshold comes from permutation.

--------------------------------------------------------------------------------
 CLADES
--------------------------------------------------------------------------------
For anchor tape a and depth d, cells sharing the depth-d prefix of tape a AND a
clone form a candidate clade; cells that did not reach depth d, or lost tape a,
are excluded.  Lambda is then evaluated for every tape z != a (cross-tape, so the
clade is never defined by the thing being tested).  Clades are deduplicated by
cell set, since many anchors induce the same partition.  Among overlapping
passing clades for one (clone, tape) the highest-Lambda one is kept, so one loss
event is one row.

⚠ Prefix-defined clades are only the clades some anchor tape happens to resolve,
so the catalogue is a LOWER BOUND on the number of events.

--------------------------------------------------------------------------------
 NULL
--------------------------------------------------------------------------------
Whole cell rows shuffled WITHIN CLONE: preserves each cell's entire missingness
profile, each tape's marginal, and the clone's own rate for that tape, destroying
only the alignment between cells and subclades.  A clone-WIDE loss is invariant
under it and so is automatically excluded from this catalogue -- correctly, since
such a loss may predate the clone.  Clone-wide cases are counted separately.
Running the identical pipeline on permuted data gives the expected number of
candidates at each threshold, hence an empirical FDR.

--------------------------------------------------------------------------------
 B = 1,000 PERMUTATIONS: parts, a fixed grid, p and q  (added 2026-09-03)
--------------------------------------------------------------------------------
The whole scan is redone per permutation, so B = 1,000 is split across array
tasks with --permpart i/N.  Two rules make the parts addable:

  * A FIXED threshold grid, independent of the data.  Each task writes only the
    vector #{candidates >= t} for t in that grid -- a few hundred integers per
    permutation -- never a candidate list.
  * Permutation b is seeded from (SEED, b), NOT from a stream advanced in order.
    So permutation b is the same object however B is split, parts can be merged
    in any order, and a duplicated part is detectable instead of silently
    doubling the null.

--permpart writes results/permparts/permcounts_40_{arm}_p{i}.npz and exits.
46_perm_merge.py pools those into results/permnull_40_{arm}.npz; the observed
run is then redone ONCE with --nullfile <that>, which gives

  * a global permutation p-value, p = (1 + #{b: C_b >= C_obs}) / (B + 1), where
    C is the total candidate count above a threshold, and
  * a per-event q-value, attached as a column of events_{arm}.tsv.gz.

Usage: 40_event_catalogue.py <arm> [--screen] [--eps 0.01] [--nperm 10]
       [--anchors N] [--lam 10] [--seed S] [--permpart i/N] [--nullfile F]
Outputs: results/events_{arm}.tsv.gz   (gitignored: carries clone barcodes)
         results/event_catalogue_{arm}.json  (aggregate, committable)
================================================================================
"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
DEPTHS, MIN_CLADE = [1, 2, 3, 4], 4

arm = sys.argv[1]
SCREEN = "--screen" in sys.argv
EPS = float(sys.argv[sys.argv.index("--eps") + 1]) if "--eps" in sys.argv else 0.01
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 10
NANCH = int(sys.argv[sys.argv.index("--anchors") + 1]) if "--anchors" in sys.argv else 166
LAM0 = float(sys.argv[sys.argv.index("--lam") + 1]) if "--lam" in sys.argv else 10.0
SEED = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260903
PART = sys.argv[sys.argv.index("--permpart") + 1] if "--permpart" in sys.argv else None
NULLF = sys.argv[sys.argv.index("--nullfile") + 1] if "--nullfile" in sys.argv else None
tag = "_screen" if SCREEN else ""

# The fixed grid.  It must NOT depend on the data: count vectors from different
# array tasks are only addable if every task scored the same thresholds.
GRID = np.unique(np.round(np.geomspace(LAM0, 1.0e7, 600), 4))


def counts_ge(v):
    """#{v >= t} for every t in GRID."""
    sv = np.sort(v)
    return (sv.size - np.searchsorted(sv, GRID, side="left")).astype(np.int64)

z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
codes = np.load(RES / f"prefix_codes_{arm}.npz", allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0]
if SCREEN:
    good = ~np.load(RES / f"collision_flags_{arm}.npz", allow_pickle=False)["flag"]
    Y, clone, codes = Y[good], clone[good], codes[good]
n, K = Y.shape
cl_names, g_clone = np.unique(clone, return_inverse=True)
Gc = g_clone.max() + 1
miss = ~Y
Rc = Y.sum(1)
print(f"{arm}{' +screen' if SCREEN else ''}: {n:,} cells, {K} tapes, {Gc:,} clones, "
      f"eps={EPS}, Lambda threshold {LAM0} nats", flush=True)


def fit_twoway(M, iters=40):
    cm = np.clip(M.mean(0), 1e-3, 1 - 1e-3)
    beta = np.log(cm / (1 - cm)); alpha = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :])))
            if ax == 0:
                g = (M - p).sum(1); h = (p * (1 - p)).sum(1)
                alpha = np.clip(alpha + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
                alpha -= alpha.mean()
            else:
                g = (M - p).sum(0); h = (p * (1 - p)).sum(0)
                beta = np.clip(beta + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(alpha[:, None] + beta[None, :]))), 1e-6, 1 - 1e-6)


P = fit_twoway(miss.astype(float))            # P(missing), cell x tape margins

# ⚠⚠ THE THIRD MARGIN, and the first version of this script was wrong without it.
# Scored against the arm-wide p_hat alone, a tape lost CLONE-WIDE makes every
# subset of that clone -- real or permuted -- look like a spectacular event, so
# such candidates flooded both counts and the FDR never fell below ~64% at any
# threshold.  Clone-wide losses are real but are not datable within the experiment
# and are exactly what the within-clone permutation reproduces.
#
# So add a per-(clone, tape) offset gamma, fitted to match the clone's own
# observed missing count for that tape, holding alpha and beta fixed:
#
#     SUM_{c in C} logit^-1(alpha_c + beta_z + gamma_{C,z}) = #missing in C
#
# the natural third margin alongside the cell and tape margins.  Consequences:
#   * a clone-wide loss gives p_tilde -> 1, hence Lambda ~ 0 for every subclade
#     -- correctly excluded rather than dominating;
#   * a loss confined to part of a clone leaves p_tilde intermediate, so the
#     subclade still scores;
#   * for a tiny clone gamma absorbs nearly everything, so no subclade can be
#     called -- conservative and correct, since with 4 cells clone-wide and
#     subclade-restricted are indistinguishable.
eta = np.log(P / (1 - P))                      # alpha_c + beta_z
Pt = np.empty_like(P)
order0 = np.argsort(g_clone, kind="stable")
bnd0 = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])
for a_, b_ in zip(bnd0[:-1], bnd0[1:]):
    ix = order0[a_:b_]
    e = eta[ix]                                # (m, K)
    tgt = miss[ix].sum(0).astype(float)        # observed missing per tape
    gam = np.zeros(e.shape[1])
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        g = p.sum(0) - tgt
        h = (p * (1 - p)).sum(0)
        gam -= np.clip(g / np.maximum(h, 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    Pt[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)
P = Pt
print(f"  fitted per-(clone,tape) offsets; clone margins matched to "
      f"max abs error {np.abs(np.array([miss[order0[a_:b_]].sum(0) - P[order0[a_:b_]].sum(0) for a_, b_ in zip(bnd0[:-1], bnd0[1:])])).max():.3g}",
      flush=True)
# per-entry contribution to Lambda: one array, so Lambda = a group sum
W = np.where(miss, np.log(1 - EPS) - np.log(P), np.log(EPS) - np.log1p(-P))
anchors = np.unique(np.linspace(0, K - 1, min(NANCH, K)).astype(int))
tape_ix = np.arange(K)


def scan(Wm, missm, collect):
    """Sweep anchors x depths; return (hits, all candidate Lambda values).

    The Lambda values of EVERY candidate are kept, for observed and permuted runs
    alike, so the FDR can be read as a function of threshold in one pass instead
    of being guessed: FDR(t) = mean #null candidates >= t / #observed >= t."""
    hits, lams = [], []
    for d in DEPTHS:
        for a_ in anchors:
            cd = codes[:, a_, d - 1]
            ok = cd >= 0
            if ok.sum() < MIN_CLADE:
                continue
            idx = np.flatnonzero(ok)
            key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
            _, sub = np.unique(key, return_inverse=True)
            G = sub.max() + 1
            ns = np.bincount(sub, minlength=G)
            L = np.empty((G, K))
            for k in range(K):
                L[:, k] = np.bincount(sub, weights=Wm[idx, k], minlength=G)
            L[:, a_] = -np.inf                      # never score the anchor itself
            L[ns < MIN_CLADE] = -np.inf
            gg, zz = np.nonzero(L >= LAM0)
            lams.append(L[gg, zz])
            if not collect or gg.size == 0:
                continue
            Mk = np.empty((G, K))
            Pk = np.empty((G, K))
            for k in range(K):
                Mk[:, k] = np.bincount(sub, weights=missm[idx, k].astype(float), minlength=G)
                Pk[:, k] = np.bincount(sub, weights=P[idx, k], minlength=G)
            owner = np.zeros(G, dtype=np.int64); owner[sub] = g_clone[idx]
            for g_, z_ in zip(gg, zz):
                cells = idx[sub == g_]
                hits.append(dict(clone=int(owner[g_]), anchor=int(a_), depth=int(d),
                                 tape=int(z_), size=int(ns[g_]), lam=float(L[g_, z_]),
                                 n_missing=int(Mk[g_, z_]), expected=float(Pk[g_, z_]),
                                 cells=cells))
    return hits, (np.concatenate(lams) if lams else np.zeros(0))


order = np.argsort(g_clone, kind="stable")
bnds = np.concatenate([[0], np.cumsum(np.bincount(g_clone))])


def perm_inv(b):
    """Row permutation for permutation index b, blocked WITHIN CLONE.

    Whole cell ROWS move; labels, prefix codes and P stay at their positions.
    Seeded from (SEED, b) rather than from a stream advanced b times, so
    permutation b is the same object however B is split across array tasks."""
    rb = np.random.default_rng([SEED, b])
    pm = order.copy()
    for a_, b_ in zip(bnds[:-1], bnds[1:]):
        pm[a_:b_] = rb.permutation(pm[a_:b_])
    inv = np.empty_like(pm)
    inv[order] = pm
    return inv


# ---- --permpart i/N: run a slice of the B permutations and write counts only
if PART is not None:
    ip, npart = (int(x) for x in PART.split("/"))
    assert 1 <= ip <= npart, PART
    lo, hi = ((ip - 1) * NPERM) // npart, (ip * NPERM) // npart
    idx = list(range(lo, hi))
    print(f"  part {ip}/{npart}: permutations {lo}..{hi-1} of B={NPERM}", flush=True)
    rows, t0 = [], time.time()
    for j, b in enumerate(idx):
        inv = perm_inv(b)
        _, nl = scan(W[inv], miss[inv], False)
        rows.append(counts_ge(nl))
        print(f"    perm {b}: {nl.size:,} candidates "
              f"[{j+1}/{len(idx)}, {time.time()-t0:.0f}s]", flush=True)
    outd = RES / "permparts"
    outd.mkdir(exist_ok=True)
    fn = outd / f"permcounts_40_{arm}{tag}_p{ip:03d}.npz"
    np.savez_compressed(fn, grid=GRID, counts=np.array(rows, dtype=np.int64),
                        perm_index=np.array(idx, dtype=np.int64),
                        meta=np.array([NPERM, npart, ip, SEED, LAM0], dtype=float))
    print(f"wrote {fn.relative_to(ROOT)}  ({len(idx)} permutations, "
          f"{time.time()-t0:.0f}s)")
    sys.exit(0)

hits, obs_lam = scan(W, miss, True)
n_cand = obs_lam.size
print(f"  candidates at Lambda >= {LAM0}: {n_cand:,}", flush=True)
grid = GRID
obs_n = counts_ge(obs_lam).astype(float)

# ---- the null, either pooled from parts or run here
if NULLF is not None:
    zn = np.load(NULLF, allow_pickle=False)
    assert np.array_equal(zn["grid"], GRID), "null file was built on a different grid"
    null_counts = zn["counts"].astype(float)                 # (B, G)
    print(f"  null: {null_counts.shape[0]:,} pooled permutations from "
          f"{Path(NULLF).name}", flush=True)
else:
    rows = []
    for b in range(NPERM):
        inv = perm_inv(b)
        _, nl = scan(W[inv], miss[inv], False)
        rows.append(counts_ge(nl))
        print(f"    perm {b+1}/{NPERM}: {nl.size:,} candidates", flush=True)
    null_counts = np.array(rows, dtype=float)
B = int(null_counts.shape[0])
nul_n = null_counts.mean(0)

# BH-style monotonisation.  An event at Lambda may be reported at ANY threshold
# t <= Lambda, so its q is the smallest FDR over those rejection regions: a
# running minimum UP the grid.  ⚠ Not down from the top -- at high thresholds
# where nothing is observed the ratio takes the convention 1.0, and a downward
# running minimum would drag that 1.0 across the whole curve.
fdr_raw = np.where(obs_n > 0, nul_n / np.maximum(obs_n, 1.0), 1.0)
fdr_curve = np.minimum(np.minimum.accumulate(fdr_raw), 1.0)
ok = np.flatnonzero(fdr_curve <= 0.05)
j0 = int(ok[0]) if ok.size else len(grid) - 1
LAM = float(grid[j0])
fdr = float(fdr_curve[j0])
# global permutation p-value: how often does a permuted scan produce as many
# candidates above the threshold as the real one?
p_at = lambda j: (1 + int((null_counts[:, j] >= obs_n[j]).sum())) / (B + 1)
p_thresh, p_floor = p_at(j0), p_at(0)
show = np.searchsorted(grid, [10, 15, 20, 30, 50, 75, 100, 200, 500, 1000])
show = np.unique(np.clip(show, 0, len(grid) - 1))
print(f"\n  FDR curve: " + "  ".join(
    f"{grid[i]:.0f}n:{100*fdr_curve[i]:.0f}%" for i in show))
print(f"  => threshold {LAM:.1f} nats at FDR {100*fdr:.1f}% "
      f"({int((obs_lam>=LAM).sum()):,} candidates, {nul_n[j0]:,.1f} expected null)")
print(f"  global permutation p (B={B:,}): {p_thresh:.4g} at the threshold, "
      f"{p_floor:.4g} at the scan floor {LAM0:.0f} nats "
      f"(obs {obs_n[0]:,.0f} vs null mean {nul_n[0]:,.1f})", flush=True)

# ---- deduplicate / enforce one row per loss event, at the chosen threshold
hits = [h for h in hits if h["lam"] >= LAM]
hits.sort(key=lambda h: -h["lam"])
kept, by_key = [], {}
for h in hits:
    k = (h["clone"], h["tape"])
    s = set(h["cells"].tolist())
    if any(s & prev for prev in by_key.get(k, [])):
        continue
    by_key.setdefault(k, []).append(s)
    kept.append(h)
print(f"  after overlap/nesting collapse: {len(kept):,} events", flush=True)

# ---- summaries
tot_missing = int(miss.sum())
explained = sum(h["n_missing"] - h["expected"] for h in kept)
ev_cells = np.unique(np.concatenate([h["cells"] for h in kept])) if kept else np.array([], int)
qj = np.searchsorted(grid, np.array([h["lam"] for h in kept], float), side="right") - 1
for h, j_ in zip(kept, np.clip(qj, 0, len(grid) - 1)):
    h["q"] = float(fdr_curve[j_])
out = {"arm": arm, "screen": SCREEN, "eps": EPS, "lambda_scan_floor": LAM0,
       "n_cells": int(n), "n_clones": int(Gc), "n_candidates": int(n_cand),
       "n_events": len(kept), "fdr": fdr, "lambda_chosen": LAM, "n_perm": B,
       "seed": SEED, "null_source": (Path(NULLF).name if NULLF else "inline"),
       "p_global_at_threshold": p_thresh, "p_global_at_scan_floor": p_floor,
       "p_floor_resolution": 1.0 / (B + 1),
       "null_mean_at_threshold": float(nul_n[j0]),
       "null_max_at_threshold": float(null_counts[:, j0].max()),
       "null_mean_at_scan_floor": float(nul_n[0]),
       "null_max_at_scan_floor": float(null_counts[:, 0].max()),
       "max_q_of_kept_events": float(max((h["q"] for h in kept), default=1.0)),
       "fdr_curve": {"lambda": grid.tolist(), "observed": obs_n.tolist(),
                     "null": nul_n.tolist(), "fdr": fdr_curve.tolist()},
       "distinct_tapes": int(len({h["tape"] for h in kept})),
       "distinct_clones": int(len({h["clone"] for h in kept})),
       "total_missing_entries": tot_missing,
       "excess_missing_explained": float(explained),
       "frac_missing_explained": float(explained / max(tot_missing, 1)),
       "cells_in_events": int(ev_cells.size),
       "lambda_total_nats": float(sum(h["lam"] for h in kept)),
       "lam_quantiles": {}, "size_quantiles": {}, "inside_rate_quantiles": {},
       "Rc_all_median": float(np.median(Rc)),
       "Rc_event_cells_median": float(np.median(Rc[ev_cells])) if ev_cells.size else None}
if kept:
    lam = np.array([h["lam"] for h in kept]); sz = np.array([h["size"] for h in kept])
    ins = np.array([h["n_missing"] / h["size"] for h in kept])
    exr = np.array([h["expected"] / h["size"] for h in kept])
    for nm, v in (("lam_quantiles", lam), ("size_quantiles", sz),
                  ("inside_rate_quantiles", ins)):
        out[nm] = {str(q): float(np.quantile(v, q)) for q in (0.1, 0.5, 0.9, 1.0)}
    out["expected_rate_median"] = float(np.median(exr))
    print(f"\n  events {len(kept):,} on {out['distinct_tapes']} tapes in "
          f"{out['distinct_clones']} clones")
    print(f"  clade size: median {np.median(sz):.0f}, max {sz.max()}")
    print(f"  inside missing rate: median {np.median(ins):.2f} vs expected "
          f"{np.median(exr):.2f}")
    print(f"  Lambda: median {np.median(lam):.1f}, max {lam.max():.1f} nats; "
          f"total {out['lambda_total_nats']:,.0f}")
    print(f"  excess missingness explained: {explained:,.0f} of {tot_missing:,} "
          f"entries = {100*out['frac_missing_explained']:.2f}%")
    print(f"  R_c median: all cells {out['Rc_all_median']:.0f}, "
          f"cells in events {out['Rc_event_cells_median']:.0f}", flush=True)

(RES / f"event_catalogue_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"events_{arm}{tag}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\tanchor_tape\tdepth\ttape\tclade_cells\tn_missing\t"
             "expected\tinside_rate\tlambda_nats\tq_value\tmedian_Rc\n")
    for h in kept:
        fh.write(f"{h['clone']}\t{cl_names[h['clone']]}\t{h['anchor']}\t{h['depth']}\t"
                 f"{h['tape']}\t{h['size']}\t{h['n_missing']}\t{h['expected']:.3f}\t"
                 f"{h['n_missing']/h['size']:.4f}\t{h['lam']:.3f}\t{h['q']:.3g}\t"
                 f"{np.median(Rc[h['cells']]):.0f}\n")
print(f"\nwrote results/event_catalogue_{arm}{tag}.json and events_{arm}{tag}.tsv.gz")
