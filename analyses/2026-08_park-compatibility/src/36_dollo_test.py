#!/usr/bin/env python3
r"""
================================================================================
 Dollo, or a heritable propensity?  -- the discriminator the gradient cannot make
================================================================================

35_lineage_depth.py showed concordance follows the topology: closer relatives
agree more about which tapes are missing.  TWO mechanisms predict that gradient
and they need different models:

  DOLLO                  the tape is lost ONCE in an ancestor, irreversibly
                         (CNV deletion, silencing of the integration site).
                         Every descendant lacks it => within a subclade below the
                         loss, the absence fraction a is 1, EXACTLY.
                         Model: one absorbing state per tape in the pruning.
                         Constant factor.

  HERITABLE PROPENSITY   a lineage merely has a lower recovery RATE at that locus;
                         cells lose it stochastically.  a sits at an intermediate
                         value and is similar across nested subclades.
                         Model: a per-lineage rate multiplier on dropout -- a
                         latent field over the tree.  Expensive; a real SBI case.

THE TEST.  For subclade S (cells sharing the depth-d prefix of an anchor tape,
within a clone) and feature tape z != anchor, let

    a = (cells of S missing z) / |S|

Dollo predicts a piles up at a = 1 EXACTLY, far beyond chance, and that the excess
SURVIVES as |S| grows -- a propensity with rate pi < 1 produces complete absence
with probability pi^|S|, which vanishes fast.  So the discriminating quantity is
the enrichment of complete absence AS A FUNCTION OF SUBCLADE SIZE.

THE NULL IS EXACT, no simulation needed for the extremes:

    P(all of S missing z)  = prod_{c in S} p_hat_cz          -> exp(sum log p_hat)
    P(none of S missing z) = prod_{c in S} (1 - p_hat_cz)

with p_hat = sigma(alpha_c + beta_z) from the same two-way fit used throughout, so
per-cell capture and per-tape quality are already accounted for.  Summing those
products over all (S, z) gives the expected COUNT under independence directly.
For the full distribution of a we draw Bernoulli realisations from the same
p_hat and histogram them the same way.

WELL-RECOVERED TAPES.  Complete absence is unremarkable for a tape recovered in 4%
of cells.  Everything is therefore also reported restricted to tapes with
beta_z >= BETA_OK, where a whole subclade missing it is surprising.

Usage: 36_dollo_test.py <arm> [--screen] [--anchors N] [--nsim K]
Output: results/dollo_{arm}{_screen}.json
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
DEPTHS = [1, 2, 3, 4]
BANDS = [(10, 20), (20, 50), (50, 200), (200, 10**9)]
BETA_OK, MIN_SUB = 0.50, 10
JUNK = 2

arm = sys.argv[1]
SCREEN = "--screen" in sys.argv
NANCH = int(sys.argv[sys.argv.index("--anchors") + 1]) if "--anchors" in sys.argv else 166
NSIM = int(sys.argv[sys.argv.index("--nsim") + 1]) if "--nsim" in sys.argv else 3
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 3
rng = np.random.default_rng(20260903)

z = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z["recovered"], z["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
codes = np.load(RES / f"prefix_codes_{arm}.npz", allow_pickle=False)["codes"]
assert codes.shape[0] == Y.shape[0]
if SCREEN:
    good = ~np.load(RES / f"collision_flags_{arm}.npz", allow_pickle=False)["flag"]
    Y, clone, codes = Y[good], clone[good], codes[good]
n, K = Y.shape
_, g_clone = np.unique(clone, return_inverse=True)
miss = (~Y).astype(np.float64)
beta_rec = Y.mean(0)
ok_tape = beta_rec >= BETA_OK
print(f"{arm}{' +screen' if SCREEN else ''}: {n:,} cells, {K} tapes, "
      f"{g_clone.max()+1:,} clones; {int(ok_tape.sum())} tapes with beta >= {BETA_OK}",
      flush=True)


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


P = fit_twoway(miss)                      # P = P(missing)
logP, log1P = np.log(P), np.log1p(-P)
NB = 20
edges = np.linspace(0, 1, NB + 1)


def gsum(x, lab, G):
    out = np.empty((G, x.shape[1]))
    for k in range(x.shape[1]):
        out[:, k] = np.bincount(lab, weights=x[:, k], minlength=G)
    return out


# "inherited" = the whole CLONE is already missing that tape, so every subclade
# inside it is trivially missing it and says nothing about a loss WITHIN the clone.
# "within" = the clone still carries the tape but this subclade has lost it all --
# the only category that evidences a loss event below the clone.
acc = {d: {b: dict(pairs=0, obs_all=0.0, exp_all=0.0, obs_none=0.0, exp_none=0.0,
                   pairs_ok=0, obs_all_ok=0.0, exp_all_ok=0.0,
                   obs_within=0.0, exp_within=0.0, obs_within_ok=0.0, exp_within_ok=0.0,
                   perm_within=0.0, perm_within_ok=0.0,
                   pairs_within=0, hist=np.zeros(NB), hist_null=np.zeros(NB),
                   hist_perm=np.zeros(NB))
           for b in range(len(BANDS))} for d in DEPTHS}
anchors = np.linspace(0, K - 1, min(NANCH, K)).astype(int)
sims = [(rng.random((n, K)) < P).astype(np.float64) for _ in range(NSIM)]

for d in DEPTHS:
    used = 0
    for zi in anchors:
        cd = codes[:, zi, d - 1]
        valid = cd >= 0
        if valid.sum() < MIN_SUB:
            continue
        key = g_clone[valid].astype(np.int64) * (cd.max() + 2) + cd[valid]
        _, sub = np.unique(key, return_inverse=True)
        G = sub.max() + 1
        ns = np.bincount(sub).astype(float)
        big = ns >= MIN_SUB
        if not big.any():
            continue
        used += 1
        idx = np.flatnonzero(valid)
        # clone-level completeness on exactly this cell set
        cl = g_clone[idx]
        _, cli = np.unique(cl, return_inverse=True)
        Gk = cli.max() + 1
        nc = np.bincount(cli).astype(float)
        clone_all = gsum(miss[idx], cli, Gk) >= nc[:, None]      # (Gk, K) bool
        owner = np.zeros(G, dtype=np.int64); owner[sub] = cli
        inherited = clone_all[owner]                             # (G, K) bool
        Sm = gsum(miss[idx], sub, G)
        Slp = gsum(logP[idx], sub, G)
        Sl1 = gsum(log1P[idx], sub, G)
        Ssim = [gsum(s[idx], sub, G) for s in sims]
        a = Sm / ns[:, None]
        a_sim = [S / ns[:, None] for S in Ssim]
        a_perm = []
        order = np.argsort(cli, kind="stable")
        bnds = np.concatenate([[0], np.cumsum(np.bincount(cli))])
        for _ in range(NPERM):
            pm = order.copy()
            for a_, b_ in zip(bnds[:-1], bnds[1:]):
                pm[a_:b_] = rng.permutation(pm[a_:b_])
            back = np.empty_like(pm); back[order] = pm
            a_perm.append(gsum(miss[idx][back], sub, G) / ns[:, None])
        pall, pnone = np.exp(Slp), np.exp(Sl1)
        for bi, (lo, hi) in enumerate(BANDS):
            rows = np.flatnonzero(big & (ns >= lo) & (ns < hi))
            if rows.size == 0:
                continue
            cols = np.arange(K) != zi
            A = acc[d][bi]
            aa = a[np.ix_(rows, cols)]
            A["pairs"] += aa.size
            A["obs_all"] += float((aa >= 1.0).sum())
            A["obs_none"] += float((aa <= 0.0).sum())
            A["exp_all"] += float(pall[np.ix_(rows, cols)].sum())
            A["exp_none"] += float(pnone[np.ix_(rows, cols)].sum())
            colsok = cols & ok_tape
            if colsok.any():
                ao = a[np.ix_(rows, colsok)]
                A["pairs_ok"] += ao.size
                A["obs_all_ok"] += float((ao >= 1.0).sum())
                A["exp_all_ok"] += float(pall[np.ix_(rows, colsok)].sum())
            # the decisive category: complete absence in a subclade whose clone
            # still carries the tape
            inh = inherited[np.ix_(rows, cols)]
            wa = (aa >= 1.0) & ~inh
            A["pairs_within"] += int((~inh).sum())
            A["obs_within"] += float(wa.sum())
            A["exp_within"] += float(pall[np.ix_(rows, cols)][~inh].sum())
            inhok = inherited[np.ix_(rows, colsok)] if colsok.any() else None
            if inhok is not None:
                aok = a[np.ix_(rows, colsok)]
                A["obs_within_ok"] += float(((aok >= 1.0) & ~inhok).sum())
                A["exp_within_ok"] += float(pall[np.ix_(rows, colsok)][~inhok].sum())
            # histogram over NON-inherited pairs only: is a bimodal at the technical
            # rate and 1, or unimodal at an intermediate value?
            for ap in a_perm:
                app = ap[np.ix_(rows, cols)]
                A["perm_within"] += float(((app >= 1.0) & ~inh).sum()) / NPERM
                A["hist_perm"] += np.histogram(app[~inh], bins=edges)[0] / NPERM
                if colsok.any():
                    apo = ap[np.ix_(rows, colsok)]
                    A["perm_within_ok"] += float(((apo >= 1.0) & ~inhok).sum()) / NPERM
            A["hist"] += np.histogram(aa[~inh], bins=edges)[0]
            for asim in a_sim:
                A["hist_null"] += np.histogram(
                    asim[np.ix_(rows, cols)][~inh], bins=edges)[0] / NSIM
    print(f"  depth {d}: {used} anchors", flush=True)

out = {"arm": arm, "screen": SCREEN, "n_cells": int(n), "beta_ok": BETA_OK,
       "n_tapes_ok": int(ok_tape.sum()), "bands": BANDS, "hist_edges": edges.tolist(),
       "depths": {}}
for d in DEPTHS:
    rec = {}
    for bi, (lo, hi) in enumerate(BANDS):
        A = acc[d][bi]
        if A["pairs"] == 0:
            continue
        e = lambda o, x: (o / x if x > 0 else float("inf"))
        rec[f"{lo}-{hi if hi < 10**8 else 'inf'}"] = {
            "pairs": int(A["pairs"]),
            "obs_all": A["obs_all"], "exp_all": A["exp_all"],
            "enrich_all": e(A["obs_all"], A["exp_all"]),
            "obs_none": A["obs_none"], "exp_none": A["exp_none"],
            "enrich_none": e(A["obs_none"], A["exp_none"]),
            "pairs_ok": int(A["pairs_ok"]), "obs_all_ok": A["obs_all_ok"],
            "exp_all_ok": A["exp_all_ok"],
            "enrich_all_ok": e(A["obs_all_ok"], A["exp_all_ok"]),
            "pairs_within": int(A["pairs_within"]),
            "obs_within": A["obs_within"], "exp_within": A["exp_within"],
            "perm_within": A["perm_within"], "perm_within_ok": A["perm_within_ok"],
            "obs_within_ok": A["obs_within_ok"], "exp_within_ok": A["exp_within_ok"],
            "hist_perm": A["hist_perm"].tolist(),
            "hist": A["hist"].tolist(), "hist_null": A["hist_null"].tolist()}
    out["depths"][str(d)] = rec

def fmt(o, x):
    if x >= 0.05:
        return f"{o/x:>9.1f}x"
    return f"  exp<0.05" if o == 0 else f" >{o/0.05:>7.0f}x"

print("\n  complete absence (a = 1) in a subclade -- observed vs the independence null")
print("  'within' = the subclade's own clone still carries the tape, so this is a")
print("  loss BELOW the clone; 'all' also counts absences inherited from the clone.\n")
print("  the comparator is PERM (whole cell rows shuffled within clone); the Bernoulli")
print("  null EXP assumes the additive-logit fit and is shown only for reference.\n")
print("  d  band       pairs_within |   within obs      perm  vs perm |      exp | "
      "beta>=%.1f: obs     perm" % BETA_OK)
for d in DEPTHS:
    for k, v in out["depths"][str(d)].items():
        pw = max(v["pairs_within"], 1)
        print(f"  {d}  {k:>9} {v['pairs_within']:>14,} | {v['obs_within']:>11,.0f} "
              f"{v['perm_within']:>9,.1f} {fmt(v['obs_within'], v['perm_within'])} | "
              f"{v['exp_within']:>8,.1f} | {v['obs_within_ok']:>9,.0f} "
              f"{v['perm_within_ok']:>8.1f}", flush=True)

tag = "_screen" if SCREEN else ""
(RES / f"dollo_{arm}{tag}.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote results/dollo_{arm}{tag}.json")
