#!/usr/bin/env python3
r"""
================================================================================
 Exact permutation moments + analytic p-values.  No calibration draws, no strata.
================================================================================

Runs on the STORED scan output of 62 -- no rescanning.

--------------------------------------------------------------------------------
 WHY THIS REPLACES THE MARGIN ROUTE
--------------------------------------------------------------------------------
63 standardised the score by the MODEL variance V = sum_S p~(1-p~) and then
compared it to the max over 1,000 permutations.  Two problems, both measured:

 1. V IS THE WRONG SCALE.  Under the PERMUTATION null with a FITTED p~ the
    standard deviation of z = (k-E)/sqrt(V) is 0.51-0.67, not 1.  (I verified
    Var = V by Monte Carlo under the MODEL null with p~ FIXED, which is the
    wrong null and could only confirm its own premise.)
 2. The margin needed a max, the max needs calibration draws, the draws need a
    budget, and the budget needs strata.  All of that is downstream of (1).

--------------------------------------------------------------------------------
 THE EXACT MOMENTS
--------------------------------------------------------------------------------
The gamma_{C,z} fit forces  sum_{c in C} (X_cz - p~_cz) = 0  per (clone, tape).
So with r_c = X_cz - p~_cz, the clone's residuals sum to zero BY CONSTRUCTION
(measured: max |sum r| = 2.1e-4).  Under the within-clone permutation a clade is
a simple random sample WITHOUT REPLACEMENT of size m from those residuals, so
with n = population size, rbar = population mean, sig2 = population variance:

    E_pi[T]   = m * rbar
    Var_pi[T] = m * (n-m)/(n-1) * sig2

⚠ The population is NOT the whole clone: a clade is drawn from the clone's cells
that reached the block's depth on its anchor tape.  So n, rbar and sig2 are
computed per (block, clone, tape).  Using the whole clone instead overstates the
SD by 3-6% (conservative) -- verified against the stored s1/s2.

The finite-population factor (n-m)/(n-1) is what V lacked: at m = n the clade is
the whole population, T is identically 0, and the variance correctly vanishes.

--------------------------------------------------------------------------------
 THE TAIL
--------------------------------------------------------------------------------
⚠⚠ EDGEWORTH WAS THE WRONG TOOL and is not used.  Its skewness correction
phi(z)*g1*(z^2-1)/6 is a CENTRAL approximation: at z=6, g1=1 the normal tail is
9.9e-10 and the correction is 3.6e-8 -- 37x the leading term, i.e. the expansion
has diverged exactly where decisions are made.

Instead, two p-values per combo, and the CONSERVATIVE (larger) one is used:

  p_norm   normal tail on the standardised score.  Uses within-clone capture
           heterogeneity (through E and sig2).  Validated against the stored
           1,000 permutations wherever they resolve: conservative by 1.04-2.78x.
  p_hyper  EXACT permutation p-value in the limit of constant alpha_c within the
           clone: the missing count among m cells drawn without replacement from
           a population of n containing kpop missing ones is Hypergeometric.
           This is Fisher's exact test on (clade membership) x (missingness).
           Discrete and skew-exact, so it does not degrade in the far tail.

--------------------------------------------------------------------------------
 MULTIPLICITY
--------------------------------------------------------------------------------
One global threshold p <= EFALSE / N_test (Bonferroni at level EFALSE), so the
expected number of false COMBOS is <= EFALSE.  Events are unions of combos and
the collapse only merges or drops, so expected false EVENTS <= EFALSE too --
the candidate-vs-event mismatch that forced the old budget simply does not arise.

⚑ The six clade-size strata are GONE and are not needed: Var_pi depends on m and
n explicitly, so p is already conditioned on clade size.  The strata were a patch
for an unstandardised statistic.

Usage: 67_exact_merge.py <arm> [--maxd 6] [--efalse 1.0] [--validate]
Outputs: results/exact_{arm}_d{D}.json / .tsv.gz
================================================================================
"""
import gzip, json, sys
from pathlib import Path
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4


def arg(f, d, c=str):
    return c(sys.argv[sys.argv.index(f) + 1]) if f in sys.argv else d


arm = sys.argv[1]
MAXD = arg("--maxd", 6, int)
EFALSE = arg("--efalse", 1.0, float)
VALIDATE = "--validate" in sys.argv
tag = f"{arm}_d{MAXD}"

# ---------------------------------------------------------------- stored scan
parts = sorted((RES / "percombo_soft").glob(f"percombo_soft_{tag}_p*.npz"))
assert parts, f"no scan parts for {tag}"
A = {}; B = 0
for f in parts:
    z = np.load(f, allow_pickle=False)
    if not A:
        for k in ("obs_k", "obs_e", "obs_h", "size", "valid"):
            A[k] = z[k]
        if VALIDATE:
            A["s1"] = z["s1"].copy(); A["s2"] = z["s2"].copy()
            A["n_ge"] = z["n_ge_z"].astype(np.int64).copy()
            A["obs_z"] = z["obs_z"]
    elif VALIDATE:
        A["s1"] += z["s1"]; A["s2"] += z["s2"]; A["n_ge"] += z["n_ge_z"]
    B += int(z["meta"][4]); del z
NC = A["size"].size
print(f"{tag}: {len(parts)} parts, {NC:,} combo slots, {A['valid'].sum():,} valid, "
      f"B={B} stored permutations{' (validation on)' if VALIDATE else ''}", flush=True)

# ---------------------------------------------------------------- refit p~
z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[sel], clone[sel]
n_, K = Y.shape
cl_names, g = np.unique(clone, return_inverse=True)
miss = ~Y
MISSF = miss.astype(np.float64)


def fit_twoway(M, iters=40):
    cm = np.clip(M.mean(0), 1e-3, 1 - 1e-3)
    be = np.log(cm / (1 - cm)); al = np.zeros(M.shape[0])
    for _ in range(iters):
        for ax in (0, 1):
            p = 1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :])))
            if ax == 0:
                gg = (M - p).sum(1); h = (p * (1 - p)).sum(1)
                al = np.clip(al + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
                al -= al.mean()
            else:
                gg = (M - p).sum(0); h = (p * (1 - p)).sum(0)
                be = np.clip(be + np.clip(gg / np.maximum(h, 1e-9), -2, 2), -12, 12)
    return np.clip(1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :]))), 1e-6, 1 - 1e-6)


P = fit_twoway(MISSF); eta = np.log(P / (1 - P))
o0 = np.argsort(g, kind="stable"); b0 = np.concatenate([[0], np.cumsum(np.bincount(g))])
for a_, b_ in zip(b0[:-1], b0[1:]):
    ix = o0[a_:b_]
    e = eta[ix]; tgt = MISSF[ix].sum(0); gam = np.zeros(K)
    for _ in range(60):
        p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
        gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
        gam = np.clip(gam, -15, 15)
    P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)
R = MISSF - P
chk = np.zeros((g.max() + 1, K))
for c in range(g.max() + 1):
    chk[c] = R[g == c].sum(0)
print(f"  gamma constraint: max |sum_c r_c| over (clone,tape) = {np.abs(chk).max():.2e}",
      flush=True)

# ---- per slot: the BLOCK-CONDITIONAL population moments and the clade's cells
npop = np.zeros(NC, np.int32); kpop = np.zeros(NC, np.int32)
rbar = np.zeros(NC, np.float64); sig2 = np.zeros(NC, np.float64)
codes = np.load(RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4
                       else f"prefix_codes_{arm}.npz"), allow_pickle=False)["codes"]
blocks, off = [], 0
for d in range(1, MAXD + 1):
    for a_ in range(K):
        cd = codes[:, a_, d - 1]; ok = cd >= 0
        if ok.sum() < MIN_CLADE:
            continue
        idx = np.flatnonzero(ok)
        key = g[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
        _, sub = np.unique(key, return_inverse=True)
        G = int(sub.max() + 1); ns = np.bincount(sub, minlength=G)
        keep = np.flatnonzero(ns >= MIN_CLADE)
        if keep.size == 0:
            continue
        # population = this block's cells, grouped by CLONE (not by clade)
        cl_ids, cl_sub = np.unique(g[idx], return_inverse=True)
        GC = cl_ids.size
        npc = np.bincount(cl_sub, minlength=GC).astype(np.float64)
        S1 = np.empty((GC, K)); S2 = np.empty((GC, K)); KP = np.empty((GC, K))
        Ri, Mi = R[idx], MISSF[idx]
        for t in range(K):
            S1[:, t] = np.bincount(cl_sub, weights=Ri[:, t], minlength=GC)
            S2[:, t] = np.bincount(cl_sub, weights=Ri[:, t] ** 2, minlength=GC)
            KP[:, t] = np.bincount(cl_sub, weights=Mi[:, t], minlength=GC)
        mb = S1 / npc[:, None]
        vb = np.maximum(S2 / npc[:, None] - mb ** 2, 0.0)
        owner_cl = np.zeros(G, np.int64); owner_cl[sub] = cl_sub      # clade -> clone slot
        oc = owner_cl[keep]
        sl = slice(off, off + keep.size * K)
        npop[sl] = np.repeat(npc[oc], K).astype(np.int32)
        kpop[sl] = KP[oc].ravel().astype(np.int32)
        rbar[sl] = mb[oc].ravel()
        sig2[sl] = vb[oc].ravel()
        blocks.append((d, a_, idx, sub, G, off, keep, ns[keep]))
        off += keep.size * K
assert off == NC, (off, NC)
print(f"  {len(blocks):,} blocks, block-conditional moments built", flush=True)

# ---------------------------------------------------------------- the test
m = A["size"].astype(np.float64)
k = A["obs_k"].astype(np.float64)
E = A["obs_e"].astype(np.float64)
T = k - E
muT = m * rbar
varT = np.where(npop > 1, m * (npop - m) / np.maximum(npop - 1, 1) * sig2, 0.0)
TESTABLE = A["valid"] & (varT > 0) & (m >= MIN_CLADE) & (npop > m)
zst = np.where(TESTABLE, (T - muT) / np.sqrt(np.maximum(varT, 1e-300)), -np.inf)
p_norm = np.where(TESTABLE, stats.norm.sf(zst), 1.0)
# exact hypergeometric tail: P(at least k missing among m drawn from npop with kpop)
p_hyp = np.ones(NC)
tm = TESTABLE
p_hyp[tm] = stats.hypergeom.sf(k[tm] - 1, npop[tm], kpop[tm], m[tm])
p_use = np.maximum(p_norm, p_hyp)                    # the conservative one
NT = int(TESTABLE.sum())
PCUT = EFALSE / max(NT, 1)
print(f"\n  testable {NT:,} of {int(A['valid'].sum()):,} valid "
      f"({100*NT/max(A['valid'].sum(),1):.1f}%)")
print(f"  Bonferroni at E[false] <= {EFALSE:g}: p <= {PCUT:.3e}")
for nm, pv in (("normal", p_norm), ("hypergeom", p_hyp), ("max (used)", p_use)):
    print(f"    {nm:>12}: {int((TESTABLE & (pv <= PCUT)).sum()):>9,} combos clear")

# ---------------------------------------------------------------- the plane
# ⚑ The old plane (63/66) had the MARGIN on y, which this route retired.  The
# natural axis now is -log10(p): the Bonferroni cut is a single horizontal line,
# and because that cut is EFALSE/N_test it differs per arm -- so the arm goes on
# the facet and the line in each facet IS the rule, exactly as clade size did
# before.  Only TESTABLE combos are binned; the rest have no p.
PI_ED = np.linspace(0.0, 1.0, 51)
NLP_ED = np.concatenate([np.linspace(0.0, 60.0, 121), [np.inf]])
SIZE_ED = np.array([4, 6, 10, 20, 50, 200, 10 ** 9])
SIZE_LAB = ["4-5", "6-9", "10-19", "20-49", "50-199", "200+"]
strat = np.clip(np.searchsorted(SIZE_ED, m, side="right") - 1, 0, 5)
pi_all = np.where(m > 0, k / np.maximum(m, 1), np.nan)
nlp = np.where(TESTABLE, -np.log10(np.maximum(p_use, 1e-300)), np.nan)
H = np.zeros((6, len(PI_ED) - 1, len(NLP_ED) - 1), np.int64)
for i in range(6):
    sel_ = TESTABLE & (strat == i)
    if sel_.any():
        H[i] = np.histogram2d(pi_all[sel_], nlp[sel_], bins=[PI_ED, NLP_ED])[0].astype(np.int64)
np.savez_compressed(RES / f"exact_{tag}_plane.npz", pi_edges=PI_ED, nlp_edges=NLP_ED,
                    size_labels=np.array(SIZE_LAB), H=H, n_testable=NT, n_valid=int(A["valid"].sum()),
                    nlp_threshold=-np.log10(PCUT), efalse=EFALSE, arm=np.array(arm))
print(f"  plane: {H.sum():,} testable combos binned; threshold at "
      f"-log10(p) = {-np.log10(PCUT):.2f}; max -log10(p) = {np.nanmax(nlp):.1f}")

# ---------------------------------------------------------------- validation
out_val = {}
if VALIDATE:
    mu_e = A["s1"] / B
    sd_e = np.sqrt(np.maximum(A["s2"] / B - mu_e ** 2, 0))
    Vm = A["obs_h"] * 0  # placeholder, not used
    gd = TESTABLE & (sd_e > 1e-9)
    # empirical permutation SD of z_V  vs  closed form / sqrt(V) ; compare on z_std scale
    p_emp = (1.0 + A["n_ge"]) / (B + 1.0)
    res = gd & (A["n_ge"] >= 1) & (A["n_ge"] <= 20)
    r1 = float(np.median(p_norm[res] / p_emp[res])) if res.sum() > 100 else None
    r2 = float(np.median(p_hyp[res] / p_emp[res])) if res.sum() > 100 else None
    print(f"\n  VALIDATION against the stored {B} permutations "
          f"({int(res.sum()):,} combos where they resolve):")
    print(f"    median p_norm  / p_emp = {r1}")
    print(f"    median p_hyper / p_emp = {r2}    (>1 = conservative)")
    out_val = dict(n_compare=int(res.sum()), ratio_norm=r1, ratio_hyper=r2)

# ---------------------------------------------------------------- collapse
PASS = TESTABLE & (p_use <= PCUT) & (T > 0)
hits = []
for d, a_, idx, sub, G, off_, keep, ns in blocks:
    sl = slice(off_, off_ + keep.size * K)
    gg, zz = np.nonzero(PASS[sl].reshape(keep.size, K))
    if not gg.size:
        continue
    o2 = np.argsort(sub, kind="stable"); cb = np.concatenate([[0], np.cumsum(np.bincount(sub, minlength=G))])
    owner = np.zeros(G, np.int64); owner[sub] = g[idx]
    for gi, z_ in zip(gg, zz):
        gk = int(keep[gi]); cells = idx[o2[cb[gk]:cb[gk + 1]]]
        s_ = off_ + gi * K + z_
        hits.append(dict(clone=int(owner[gk]), depth=d, anchor=int(a_), tape=int(z_),
                         size=int(m[s_]), nmiss=int(k[s_]), pi=float(k[s_] / m[s_]),
                         eb=float(E[s_] / m[s_]), z=float(zst[s_]),
                         pn=float(p_norm[s_]), ph=float(p_hyp[s_]), pu=float(p_use[s_]),
                         lam_h=float(A["obs_h"][s_]), npop=int(npop[s_]),
                         cells=cells))
kept, seen = [], {}
for h in sorted(hits, key=lambda h: -h["lam_h"]):
    key = (h["clone"], h["tape"]); ss = set(h["cells"].tolist())
    if any(ss & pv for pv in seen.get(key, [])):
        continue
    seen.setdefault(key, []).append(ss); kept.append(h)
tot = int(miss.sum()); cxt = sum(h["nmiss"] for h in kept)
print(f"\n  {len(hits):,} combos pass -> {len(kept):,} events "
      f"({len({(h['clone'],h['tape']) for h in kept})} distinct clone x tape)")
print(f"  cell x tape entries inside them: {cxt:,} of {tot:,} = {100*cxt/max(tot,1):.2f}%")
if kept:
    pi_ = np.array([h["pi"] for h in kept]); sz_ = np.array([h["size"] for h in kept])
    print(f"  pi_hat median {np.median(pi_):.3f} (expected {np.median([h['eb'] for h in kept]):.3f}); "
          f">=0.99 {100*(pi_>=0.99).mean():.1f}%  <0.90 {100*(pi_<0.90).mean():.1f}%")
    print(f"  smallest called clade {int(sz_.min())} cells; "
          f"max p used {max(h['pu'] for h in kept):.2e}")
    print(f"\n  {'clade':>10} {'events':>7} {'med pi':>8} {'med z':>7} {'med p_used':>12}")
    for lo, hi in [(4,6),(6,10),(10,20),(20,50),(50,200),(200,10**9)]:
        s_ = [h for h in kept if lo <= h["size"] < hi]
        if s_:
            print(f"  {lo:>4}-{hi-1:<5} {len(s_):>7} {np.median([h['pi'] for h in s_]):>8.3f} "
                  f"{np.median([h['z'] for h in s_]):>7.2f} "
                  f"{np.median([h['pu'] for h in s_]):>12.2e}")

out = dict(arm=arm, max_depth=MAXD, efalse=EFALSE, p_threshold=PCUT,
           n_slots=int(NC), n_valid=int(A["valid"].sum()), n_testable=NT,
           n_combos_pass=len(hits), n_events=len(kept),
           cellxtape_in_events=int(cxt), total_missing=tot,
           frac_missing_in_events=cxt / max(tot, 1),
           method="score standardised by exact finite-population permutation moments; "
                  "p = max(normal, hypergeometric); Bonferroni",
           validation=out_val)
if kept:
    out.update(pi_median=float(np.median(pi_)), pi_ge_099=float((pi_ >= 0.99).mean()),
               pi_lt_090=float((pi_ < 0.90).mean()), min_clade=int(sz_.min()),
               max_p_used=float(max(h["pu"] for h in kept)))
(RES / f"exact_{tag}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"exact_{tag}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\tdepth\tanchor\ttape\tclade_cells\tn_pop\tn_missing\t"
             "pi_hat\texpected_rate\tz_std\tp_normal\tp_hyper\tp_used\tlambda_hard\n")
    for h in sorted(kept, key=lambda h: h["pu"]):
        fh.write(f"{h['clone']}\t{cl_names[h['clone']]}\t{h['depth']}\t{h['anchor']}\t"
                 f"{h['tape']}\t{h['size']}\t{h['npop']}\t{h['nmiss']}\t{h['pi']:.4f}\t"
                 f"{h['eb']:.4f}\t{h['z']:.3f}\t{h['pn']:.3e}\t{h['ph']:.3e}\t"
                 f"{h['pu']:.3e}\t{h['lam_h']:.3f}\n")
print(f"\nwrote results/exact_{tag}.json / .tsv.gz")
