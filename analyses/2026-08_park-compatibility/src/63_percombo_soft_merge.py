#!/usr/bin/env python3
r"""
================================================================================
 Merge 62's parts: calibrate the score margin, call events, characterise them
================================================================================

Input : results/percombo_soft/percombo_soft_{arm}_d{D}_p*.npz   (accumulators)
        results/percombo_soft/percombo_soft_{arm}_d{D}_cal.npz  (calib draws)

THE MARGIN, per (clade, tape) combo, for a statistic T in {z, L_soft, L_hard}:

    m_obs    = T_obs   - M,     M = max over permutations 0..B-1
    m_cal[j] = T_{B+j} - M      (fresh permutations, j = 0..NCAL-1)

Under H0 the calibration draws are exchangeable with the observed value given M
(swapping "observed" for "draw B+j" leaves M untouched, since M depends only on
0..B-1), so  expected false(t) = mean_j #{ m_cal[j] >= t }  is unbiased, and
    FDR(t) = expected false(t) / #{ m_obs >= t },  q(t) = min_{t'<=t} FDR(t').

⚠ NCAL IS THE ONLY SOURCE OF PRECISION on the null count: the B accumulated
permutations build M, they are not null draws of the margin.  SE ~ sd/sqrt(NCAL).
Raise it by rerunning `62 --calib N` alone -- nothing else need be redone.

⚠ THRESHOLD PER CLADE-SIZE STRATUM at an ABSOLUTE BUDGET of expected-false
candidates, not a rate: the null conditions on clade size but a GLOBAL threshold
does not, and a rate computed on CANDIDATES cannot bound the error on EVENTS
after the overlap collapse.

--------------------------------------------------------------------------------
 DIVISION OF LABOUR -- three jobs, three statistics, one scan
--------------------------------------------------------------------------------
 DETECT      z = (k-E)/sqrt(V), the score test of delta=0 in the nested family
             X_c ~ Bern(sigma(eta_c + delta)).  Agnostic to completeness, so
             pi_hat can then be REPORTED rather than assumed.  One-sidedness is
             automatic (only T_obs above every permutation earns a positive
             margin); a called event additionally needs z_obs > 0.
 ATTRIBUTE   Lambda_hard.  For pi_hat ~ 1 it grows ~ m*(log(1-eps) - log p~),
             linearly in m, while stepping up to a parent clade with present
             cells costs ~4 nats each -- so it is maximised at the LARGEST clade
             that is still COMPLETE, which is where a Dollo loss on a stem
             belongs.  Hence the overlap collapse is ordered by Lambda_hard, not
             by the detection margin.  (⚠ 63's first version collapsed by the
             margin -- the detection statistic -- which is exactly the
             over-attribution the whole design is meant to avoid.)
 REPORT      pi_hat = k/m and ebar = E/m; plus, on the CALLED shortlist only,
             the exact LRT 2[l(delta_hat) - l(0)] and delta_hat by Newton.
             ⚠ delta_hat DIVERGES for a complete loss (complete separation), so
             it is the effect size for GRADED shifts and pi_hat is the readout
             for complete ones -- the divergence is the signature of a Dollo
             loss, not a numerical failure.

All three detectors are calibrated here on the SAME permutations and reported
side by side, so the correction to the Lambda_soft version is visible in the
output rather than only in the notes.

Usage: 63_percombo_soft_merge.py <arm> [--maxd 6] [--budget 2] [--pimin 0]
Outputs: results/percombo_score_{arm}_d{D}.json        (aggregate, committable)
         results/percombo_score_{arm}_d{D}.tsv.gz      (events; GITIGNORED)
         results/percombo_score_{arm}_d{D}_plane.npz   (the figure's 2-D density)
================================================================================
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_CLADE = 4


def arg(flag, default, cast=str):
    return cast(sys.argv[sys.argv.index(flag) + 1]) if flag in sys.argv else default


arm = sys.argv[1]
MAXD = arg("--maxd", 6, int)
BUDGET = arg("--budget", 2.0, float)
PIMIN = arg("--pimin", 0.0, float)
tag = f"{arm}_d{MAXD}"

parts = sorted((RES / "percombo_soft").glob(f"percombo_soft_{tag}_p*.npz"))
assert parts, f"no accumulator parts for {tag}"
A = {}
B = 0
for f in parts:
    z = np.load(f, allow_pickle=False)
    if not A:
        for key in ("obs_z", "obs_s", "obs_h", "obs_k", "obs_e", "obs_v", "size", "valid"):
            A[key] = z[key]
        for key in ("mx_z", "mx_s", "mx_h"):
            A[key] = z[key].copy()
        A["n_ge_z"] = z["n_ge_z"].astype(np.int32).copy()
        A["s1"] = z["s1"].copy(); A["s2"] = z["s2"].copy()
    else:
        assert np.array_equal(z["obs_z"], A["obs_z"]), f"{f.name}: observed differs"
        assert np.array_equal(z["valid"], A["valid"]), f"{f.name}: layout differs"
        for key in ("mx_z", "mx_s", "mx_h"):
            np.maximum(A[key], z[key], out=A[key])
        A["n_ge_z"] += z["n_ge_z"]; A["s1"] += z["s1"]; A["s2"] += z["s2"]
    B += int(z["meta"][4]); EPS = float(z["meta"][7])
    del z
calf = RES / "percombo_soft" / f"percombo_soft_{tag}_cal.npz"
assert calf.exists(), f"no calibration draws: run 62 --calib N on {arm}"
zc = np.load(calf, allow_pickle=False)
NCAL = zc["z"].shape[0]
assert int(zc["meta"][0]) == B, (
    f"calib draws assumed B={int(zc['meta'][0])} but {B} were accumulated -- the "
    f"calibration permutations must be DISJOINT from the accumulated set")
size, valid = A["size"], A["valid"]
print(f"{tag}: {len(parts)} parts, B={B:,} accumulated, {NCAL} calibration draws, "
      f"{size.size:,} slots ({valid.sum():,} scorable); "
      f"comparison draws: soft {zc['lam_s'].shape[0]}, hard {zc['lam_h'].shape[0]}", flush=True)

SIZE_EDGES = np.array([4, 6, 10, 20, 50, 200, 10 ** 9])
SIZE_LABEL = ["4-5", "6-9", "10-19", "20-49", "50-199", "200+"]
NB = len(SIZE_LABEL)
strat = np.clip(np.searchsorted(SIZE_EDGES, size, side="right") - 1, 0, NB - 1)
pi = np.where(size > 0, A["obs_k"].astype(np.float64) / np.maximum(size, 1), np.nan)
eb = np.where(size > 0, A["obs_e"].astype(np.float64) / np.maximum(size, 1), np.nan)

GRID = {"z": np.unique(np.round(np.concatenate(
            [np.arange(0.0, 3.0, 0.02), np.geomspace(3.0, 200.0, 300)]), 4)),
        "soft": np.unique(np.round(np.concatenate(
            [np.arange(0.0, 5.0, 0.05), np.geomspace(5.0, 2000.0, 300)]), 4))}
GRID["hard"] = GRID["soft"]
CALKEY = {"z": "z", "soft": "lam_s", "hard": "lam_h"}
OBSKEY = {"z": "obs_z", "soft": "obs_s", "hard": "obs_h"}
MXKEY = {"z": "mx_z", "soft": "mx_s", "hard": "mx_h"}


def margins(stat):
    """⚠ z carries all NCAL calibration draws; the two comparison statistics carry
    only the first NCMP (memory -- see 62).  Use whatever is stored for each."""
    ok = valid & np.isfinite(A[OBSKEY[stat]]) & (A[MXKEY[stat]] > -1e29)
    mo = np.where(ok, A[OBSKEY[stat]] - A[MXKEY[stat]], -np.inf)
    nd = zc[CALKEY[stat]].shape[0]
    mc = [np.where(ok, zc[CALKEY[stat]][j] - A[MXKEY[stat]], -np.inf) for j in range(nd)]
    return ok, mo, mc


def thresholds(stat, mo, mc, ok, verbose=False):
    g = GRID[stat]; thr = np.full(NB, np.inf); js = np.full(NB, -1, int)
    ef = np.zeros(NB)
    for i in range(NB):
        m = ok & (strat == i)
        if not m.any():
            continue
        o = np.array([(mo[m] >= t).sum() for t in g], float)
        u = np.array([np.mean([(v[m] >= t).sum() for v in mc]) for t in g])
        c = np.flatnonzero(u <= BUDGET)
        if c.size and o[int(c[0])] > 0:
            js[i] = int(c[0]); thr[i] = g[js[i]]; ef[i] = u[js[i]]
            if verbose:
                print(f"   {SIZE_LABEL[i]:>8} {int(m.sum()):>12,} {thr[i]:>10.3f} "
                      f"{o[js[i]]:>10,.0f} {u[js[i]]:>10.2f}")
        elif verbose:
            print(f"   {SIZE_LABEL[i]:>8} {int(m.sum()):>12,} {'NONE':>10} "
                  f"{'--':>10} {'--':>10}")
    return thr, ef


print(f"\n  per-stratum threshold on the SCORE margin, expected false <= {BUDGET:g}:")
print(f"   {'clade':>8} {'combos':>12} {'thr':>10} {'above':>10} {'exp.false':>10}")
okz, mo_z, mc_z = margins("z")
THRS, EFS = thresholds("z", mo_z, mc_z, okz, verbose=True)
tie_o = float((mo_z[okz] >= 0).mean())
tie_u = float(np.mean([(v[okz] >= 0).mean() for v in mc_z]))
print(f"  tie diagnostic: P(margin>=0) observed {tie_o:.4f} vs calibration {tie_u:.4f} "
      f"(a continuous null would give {1/(B+1):.5f})")

# ---- all three detectors, same permutations, slot level
print(f"\n  {'detector':>10} {'above thr':>11} {'exp.false':>10} {'med pi':>8} "
      f"{'pi>=0.99':>9} {'clades>=50':>11} {'of those z<0.5':>15}")
CMP = {}
for stat in ("z", "soft", "hard"):
    ok_, mo_, mc_ = margins(stat)
    thr_, ef_ = thresholds(stat, mo_, mc_, ok_)
    hit = ok_ & (mo_ >= thr_[strat])
    if stat == "z":
        hit &= A["obs_z"] > 0
    nb = int(hit.sum()); big = hit & (size >= 50)
    art = int((big & (A["obs_z"] < 0.5)).sum())
    CMP[stat] = dict(above=nb, exp_false=float(ef_.sum()),
                     med_pi=float(np.median(pi[hit])) if nb else None,
                     pi99=float((pi[hit] >= 0.99).mean()) if nb else None,
                     n_big=int(big.sum()), n_big_lowz=art,
                     thresholds={SIZE_LABEL[i]: (float(thr_[i]) if np.isfinite(thr_[i]) else None)
                                 for i in range(NB)})
    print(f"  {stat:>10} {nb:>11,} {ef_.sum():>10.1f} "
          f"{(np.median(pi[hit]) if nb else float('nan')):>8.3f} "
          f"{(100*(pi[hit]>=0.99).mean() if nb else float('nan')):>8.1f}% "
          f"{int(big.sum()):>11,} {art:>15,}")

# ---- the plane: 2-D density of (pi_hat, score margin), observed and null
PI_ED = np.linspace(0.0, 1.0, 51)
MG_ED = np.concatenate([[-1e30], np.linspace(-8.0, 40.0, 97), [1e30]])
ND = zc["k"].shape[0]
H_obs = np.zeros((NB, len(PI_ED) - 1, len(MG_ED) - 1), np.int64)
H_cal = np.zeros_like(H_obs)
for i in range(NB):
    m = okz & (strat == i)
    if not m.any():
        continue
    H_obs[i] = np.histogram2d(pi[m], mo_z[m], bins=[PI_ED, MG_ED])[0].astype(np.int64)
    for j in range(ND):
        pj = zc["k"][j][m].astype(np.float64) / np.maximum(size[m], 1)
        H_cal[i] += np.histogram2d(pj, mc_z[j][m], bins=[PI_ED, MG_ED])[0].astype(np.int64)
np.savez_compressed(RES / f"percombo_score_{tag}_plane.npz",
                    pi_edges=PI_ED, margin_edges=MG_ED,
                    size_labels=np.array(SIZE_LABEL), H_obs=H_obs, H_cal=H_cal,
                    n_cal_draws=ND, thresholds=THRS, grid=GRID["z"])
print(f"\n  wrote the plane ({H_obs.sum():,} combos binned, {ND} null draws)")

# ---- rebuild the layout: slot -> (clone, depth, anchor, tape, cells)
z0 = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
Y, clone = z0["recovered"], z0["clone"]
s_ = (Y.sum(1) >= THR[arm]) & (clone >= 0)
Y, clone = Y[s_], clone[s_]
codes = np.load(RES / (f"prefix_codes6_{arm}.npz" if MAXD > 4
                       else f"prefix_codes_{arm}.npz"), allow_pickle=False)["codes"]
n, K = Y.shape
cl_names, g_clone = np.unique(clone, return_inverse=True)
miss = ~Y
THRM = THRS[strat]
PASS = okz & (mo_z >= THRM) & (A["obs_z"] > 0)
hits, off = [], 0
for d in range(1, MAXD + 1):
    for a_ in range(K):
        cd = codes[:, a_, d - 1]; ok = cd >= 0
        if ok.sum() < MIN_CLADE:
            continue
        idx = np.flatnonzero(ok)
        key = g_clone[idx].astype(np.int64) * (int(cd[idx].max()) + 2) + cd[idx]
        _, sub = np.unique(key, return_inverse=True)
        G = int(sub.max() + 1); ns = np.bincount(sub, minlength=G)
        keep = np.flatnonzero(ns >= MIN_CLADE)
        if keep.size == 0:
            continue
        sl0 = slice(off, off + keep.size * K)
        gg, zz = np.nonzero(PASS[sl0].reshape(keep.size, K))
        if gg.size:
            o2 = np.argsort(sub, kind="stable"); cb = np.concatenate([[0], np.cumsum(ns)])
            owner = np.zeros(G, np.int64); owner[sub] = g_clone[idx]
            for gi, z_ in zip(gg, zz):
                gk = int(keep[gi]); cells = idx[o2[cb[gk]:cb[gk + 1]]]
                sl = off + gi * K + z_
                hits.append(dict(clone=int(owner[gk]), depth=d, anchor=int(a_),
                                 tape=int(z_), size=int(ns[gk]),
                                 z=float(A["obs_z"][sl]), margin=float(mo_z[sl]),
                                 lam_h=float(A["obs_h"][sl]), lam_s=float(A["obs_s"][sl]),
                                 nullmax=float(A["mx_z"][sl]), n_ge=int(A["n_ge_z"][sl]),
                                 pi=float(pi[sl]), eb=float(eb[sl]), V=float(A["obs_v"][sl]),
                                 strat=SIZE_LABEL[strat[sl]], n_missing=int(A["obs_k"][sl]),
                                 cells=cells))
        off += keep.size * K
assert off == size.size, (off, size.size)

# ---- overlap collapse ordered by Lambda_hard: the largest COMPLETE clade wins
def collapse(hs, keyfn):
    kept, seen = [], {}
    for h in sorted(hs, key=keyfn):
        k_ = (h["clone"], h["tape"]); ss = set(h["cells"].tolist())
        if any(ss & prev for prev in seen.get(k_, [])):
            continue
        seen.setdefault(k_, []).append(ss); kept.append(h)
    return kept


kept = collapse(hits, lambda h: -h["lam_h"])
alt = collapse(hits, lambda h: -h["margin"])
sel = [h for h in kept if h["pi"] >= PIMIN]
tot = int(miss.sum()); cxt = sum(h["n_missing"] for h in sel)
print(f"\n  {len(hits):,} combos pass -> {len(kept):,} events collapsed by Lambda_hard "
      f"({len(alt):,} if collapsed by the margin instead)")
print(f"  cell x tape entries inside them: {cxt:,} of {tot:,} = {100*cxt/max(tot,1):.2f}%")

# ---- exact LRT and delta_hat on the shortlist only
out_ev = []
if sel:
    def fit_twoway(M, iters=40):
        cm = np.clip(M.mean(0), 1e-3, 1 - 1e-3)
        be = np.log(cm / (1 - cm)); al = np.zeros(M.shape[0])
        for _ in range(iters):
            for ax in (0, 1):
                p = 1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :])))
                if ax == 0:
                    g = (M - p).sum(1); h = (p * (1 - p)).sum(1)
                    al = np.clip(al + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
                    al -= al.mean()
                else:
                    g = (M - p).sum(0); h = (p * (1 - p)).sum(0)
                    be = np.clip(be + np.clip(g / np.maximum(h, 1e-9), -2, 2), -12, 12)
        return np.clip(1.0 / (1.0 + np.exp(-(al[:, None] + be[None, :]))), 1e-6, 1 - 1e-6)

    P = fit_twoway(miss.astype(float))
    eta = np.log(P / (1 - P))
    for c in range(g_clone.max() + 1):
        ix = np.flatnonzero(g_clone == c)
        e = eta[ix]; tgt = miss[ix].sum(0).astype(float); gam = np.zeros(K)
        for _ in range(60):
            p = 1.0 / (1.0 + np.exp(-(e + gam[None, :])))
            gam -= np.clip((p.sum(0) - tgt) / np.maximum((p * (1 - p)).sum(0), 1e-9), -3, 3)
            gam = np.clip(gam, -15, 15)
        P[ix] = np.clip(1.0 / (1.0 + np.exp(-(e + gam[None, :]))), 1e-6, 1 - 1e-6)
    ETA = np.log(P / (1 - P))
    for h in sel:
        cc = h["cells"]; t = h["tape"]
        X = miss[cc, t].astype(float); et = ETA[cc, t]
        d = 0.0
        for _ in range(200):
            q = 1.0 / (1.0 + np.exp(-(et + d)))
            g = (X - q).sum(); hh = (q * (1 - q)).sum()
            if hh < 1e-12:
                break
            st = g / hh; d += float(np.clip(st, -5, 5))
            if abs(st) < 1e-10:
                break
        q = 1.0 / (1.0 + np.exp(-(et + d)))
        l1 = float((X * np.log(np.clip(q, 1e-300, 1)) +
                    (1 - X) * np.log(np.clip(1 - q, 1e-300, 1))).sum())
        p0 = P[cc, t]
        l0 = float((X * np.log(p0) + (1 - X) * np.log1p(-p0)).sum())
        h["delta"] = d; h["LRT"] = 2 * (l1 - l0)
        out_ev.append(h)
    dd = np.array([h["delta"] for h in sel]); LR = np.array([h["LRT"] for h in sel])
    p_ = np.array([h["pi"] for h in sel]); sz_ = np.array([h["size"] for h in sel])
    zz_ = np.array([h["z"] for h in sel]); dp_ = np.array([h["depth"] for h in sel])
    print(f"  z: median {np.median(zz_):.2f}, min {zz_.min():.2f}; "
          f"exact LRT median {np.median(LR):.2f}; delta_hat median {np.median(dd):.2f} "
          f"({100*(dd>10).mean():.0f}% diverged, i.e. complete separation)")
    print(f"  pi_hat: median {np.median(p_):.3f} (expected {np.median([h['eb'] for h in sel]):.3f}); "
          f">=0.99 {100*(p_>=0.99).mean():.1f}%  <0.90 {100*(p_<0.90).mean():.1f}%")
    print(f"  smallest called clade {sz_.min()} cells; all beat every permutation: "
          f"{all(h['n_ge'] == 0 for h in sel)}")
    print(f"\n  {'clade size':>12} {'n':>5} {'med pi':>8} {'med z':>7} {'med LRT':>9} "
          f"{'med L_hard':>11}")
    for lo, hi in [(4, 6), (6, 10), (10, 20), (20, 50), (50, 200), (200, 10**9)]:
        m = (sz_ >= lo) & (sz_ < hi)
        if m.sum():
            print(f"   {lo:>5}-{hi-1:<6} {int(m.sum()):>5} {np.median(p_[m]):>8.3f} "
                  f"{np.median(zz_[m]):>7.2f} {np.median(LR[m]):>9.2f} "
                  f"{np.median([h['lam_h'] for h, k_ in zip(sel, m) if k_]):>11.1f}")

out = dict(arm=arm, max_depth=MAXD, B=B, n_calib=NCAL, eps=EPS, budget=BUDGET,
           pimin=PIMIN, detector="score z=(k-E)/sqrt(V)",
           attributor="Lambda_hard (largest complete clade)",
           n_combo_slots=int(size.size), n_scorable=int(okz.sum()),
           tie_obs=tie_o, tie_calib=tie_u, continuous_null=1 / (B + 1),
           size_labels=SIZE_LABEL,
           strat_thresholds={SIZE_LABEL[i]: (float(THRS[i]) if np.isfinite(THRS[i]) else None)
                             for i in range(NB)},
           exp_false_total=float(EFS.sum()), detector_comparison=CMP,
           n_combos_above=len(hits), n_events=len(kept),
           n_events_margin_collapse=len(alt), n_events_pimin=len(sel),
           cellxtape_in_events=int(cxt), total_missing_entries=tot,
           frac_missing_in_events=cxt / max(tot, 1))
if sel:
    out.update(z_median=float(np.median(zz_)), z_min=float(zz_.min()),
               LRT_median=float(np.median(LR)), delta_median=float(np.median(dd)),
               frac_delta_diverged=float((dd > 10).mean()),
               pi_median=float(np.median(p_)), pi_ge_099=float((p_ >= 0.99).mean()),
               pi_lt_090=float((p_ < 0.90).mean()), min_called_clade=int(sz_.min()),
               all_beat_every_permutation=bool(all(h["n_ge"] == 0 for h in sel)),
               pi_by_depth={int(d): dict(n=int((dp_ == d).sum()),
                                         median_pi=float(np.median(p_[dp_ == d])),
                                         median_z=float(np.median(zz_[dp_ == d])))
                            for d in range(1, MAXD + 1) if (dp_ == d).sum()})
(RES / f"percombo_score_{tag}.json").write_text(json.dumps(out, indent=1))
with gzip.open(RES / f"percombo_score_{tag}.tsv.gz", "wt") as fh:
    fh.write("clone\tclone_bc\tdepth\tanchor\ttape\tclade_cells\tn_missing\tpi_hat\t"
             "expected_rate\tz\tV\tmargin\tnull_max\tdelta_hat\tLRT\tlambda_hard\t"
             "lambda_soft\tn_ge\tsize_stratum\n")
    for h in sorted(out_ev, key=lambda h: -h["z"]):
        fh.write(f"{h['clone']}\t{cl_names[h['clone']]}\t{h['depth']}\t{h['anchor']}\t"
                 f"{h['tape']}\t{h['size']}\t{h['n_missing']}\t{h['pi']:.4f}\t"
                 f"{h['eb']:.4f}\t{h['z']:.3f}\t{h['V']:.3f}\t{h['margin']:.3f}\t"
                 f"{h['nullmax']:.3f}\t{h['delta']:.3f}\t{h['LRT']:.3f}\t"
                 f"{h['lam_h']:.3f}\t{h['lam_s']:.3f}\t{h['n_ge']}\t{h['strat']}\n")
print(f"\nwrote results/percombo_score_{tag}.json / .tsv.gz / _plane.npz")
