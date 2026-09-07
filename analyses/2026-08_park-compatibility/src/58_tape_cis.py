#!/usr/bin/env python3
r"""
================================================================================
 58 -- PREMISE CHECK, done on writes: is symbol composition TAPE-dependent?
================================================================================

⚠⚠ THIS SUPERSEDES THE per-tape BLOCK OF 56, WHICH WAS CONFOUNDED.  56 counted
raw (cell, tape, site) entries and took its floor from a random split of those
entries.  Both halves of that are wrong for the same reason: a symbol written at
tape y in a clone founder reappears in every descendant, so entries are not
independent draws.  The resulting "excess" tracked CLONAL CONCENTRATION, not the
tapes -- median TV 0.040 in Pre-TX (many small clones) against 0.563 in Mouse 2
(one clone holding 98.9% of the weight) -- and the per-SITE statistic moved the
same way, which settles it, since a site cannot be cis to a tape.

WHY IT MATTERS.  The co-integrated-symbol test assumes pegRNAs act in TRANS, so
xi must be identical at all 166 tapes.  If any CIS component existed, then --
because tape losses are correlated with each other (rho_tape = 0.25, fig 3b) --
a clade that lost tape z would carry a non-random tape mixture and so a shifted
symbol composition, with no silencing at all.  That is the very signal we intend
to claim, manufactured.  So the premise has to be measured, on the right unit.

--------------------------------------------------------------------------------
 THE DESIGN
--------------------------------------------------------------------------------
UNIT: independent post-MRCA writes within a clone (57_writes.extract_writes).
Ancestral content is monomorphic in the clone and contributes nothing, so
identity by descent is gone by construction.

CONTRAST, clone-stratified: for tape t, compare t's write composition against the
composition of all OTHER tapes IN THE SAME CLONES, each clone weighted by its
share of tape t's writes.  Clone founder effects therefore cancel in the contrast
rather than being assumed absent.

    p_ts   = SUM_c n_cts / SUM_c n_ct.
    ref_ts = SUM_c (n_ct. / N_t.) * (M_cs - n_cts) / (M_c. - n_ct.)
    TV_t   = 0.5 * SUM_s |p_ts - ref_ts|

NULL: within each (clone, site) block, permute the TAPE LABEL among writes.  This
preserves every clone's composition, every site's composition, and every
(clone, tape) write count, and destroys only which tape a symbol landed on --
exactly the null "composition is tape-independent".  No re-extraction is needed,
because the writes themselves are defined by prefix structure, not by tape
identity.  ⚠ The site block matters: site 6 has elevated q (README), so a null
that mixed sites would be too wide.

REPORTED: TV_t against its own null mean and its null 95th percentile; the number
of tapes exceeding; and a single global statistic SUM_t N_t * TV_t with a
permutation p, so the question gets one answer and not 166.

--------------------------------------------------------------------------------
 ⭐ THE PER-SYMBOL DECOMPOSITION -- why the residual is not just a nuisance
--------------------------------------------------------------------------------
Measured 2026-09-07: the pseudoreplication was ~98% of the apparent effect, but a
small residual survives at p = 1/201 in EVERY arm, and its relative size orders by
statistical power (Mouse 3 +2.6% -> Pre-TX ~+21%).  Two readings, distinguishable
only by the SHAPE of the excess, which the TV summary throws away:

 (1) IT IS THE MECHANISM, at cell level.  Tape t's writes come only from cells
     that RECOVERED tape t -- i.e. cells where integration t is expressed.  In
     exactly those cells pegRNA t is expressed too, so symbol s(t) is over-
     supplied.  Trans action writes it to every tape, but the comparison draws on
     a partly different cell population.  So a genuine trans co-integration
     mechanism PREDICTS this residual: a small, positive, tape-specific
     enrichment of s(t) at tape t.  If so this is not a premise check at all, it
     is the test, run clade-free on every write in the arm -- and the recovered
     argmax IS the z -> s(z) map, checkable for five-way cross-arm concordance.
     ⇒ CONCENTRATED on one symbol, POSITIVE.

 (2) IT IS HOMOPLASY-COLLAPSE BIAS.  The dedup merges two independent writes of
     the same symbol under the same parent, so a tape with a bushier trie
     collapses more, and collapse hits FREQUENT symbols hardest.  The null
     preserves each (clone, tape) write count but reassigns symbols, so it does
     not reproduce tape-specific collapse.
     ⇒ DIFFUSE, over the frequent-symbol end, and CORRELATED WITH BUSHINESS.

So this script now records, per tape, the signed per-symbol excess in null-sd
units, the concentration (top |z| against the rms of the rest), and each tape's
writes-per-parent as the bushiness covariate.  Shape settles which reading holds;
corr(excess, bushiness) is the direct check on (2).

Usage: 58_tape_cis.py <arm> [--nperm 200] [--poly codes|codes_or_none]
Output: results/tape_cis_{arm}.json
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_w = __import__("57_writes")
load_arm, extract_writes = _w.load_arm, _w.extract_writes

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

arm = sys.argv[1]
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 200
POLY = sys.argv[sys.argv.index("--poly") + 1] if "--poly" in sys.argv else "codes"
MIN_W = 50                      # clones contributing fewer writes carry no signal
rng = np.random.default_rng(20260907)

Y, clone, codes, S, syms, tapes = load_arm(arm)
# ⚠ composition is measured on DESIGN-CONFORMING symbols only -- see
# 57_writes.conforming for why junk dominated the first run of this scan.
_cmask, _remap, syms = _w.conforming(syms)
K, A = len(tapes), len(syms)
print(f"  alphabet: {int(_cmask.sum())} conforming of {len(_cmask)} symbols", flush=True)

# ---- collect writes, clone by clone ---------------------------------------
names, sizes = np.unique(clone, return_counts=True)
w_c, w_t, w_j, w_s, w_p = [], [], [], [], []
for ci, cl in enumerate(names):
    if sizes[ci] < 2:
        continue
    g = clone == cl
    t, j, s, pk = extract_writes(codes[g], S[g], Y[g], poly=POLY, parents=True)
    _k = _remap[s] >= 0
    t, j, s, pk = t[_k], j[_k], _remap[s[_k]], pk[_k]
    if len(t) < MIN_W:
        continue
    w_c.append(np.full(len(t), ci, dtype=np.int64)); w_t.append(t)
    w_j.append(j); w_s.append(s)
    w_p.append(pk * (len(names) + 1) + ci)      # parent keys are clone-local
assert w_c, f"{arm}: no clone reached {MIN_W} writes"
w_c, w_t, w_j, w_s, w_p = (np.concatenate(x) for x in (w_c, w_t, w_j, w_s, w_p))
uc, w_c = np.unique(w_c, return_inverse=True)
NC, NW = len(uc), len(w_t)
print(f"{arm}: {NW:,} independent writes over {NC:,} clones "
      f"({sizes[np.isin(np.arange(len(names)), uc)].sum():,} cells), "
      f"{A} symbols, poly={POLY}", flush=True)

M_cs = np.bincount(w_c * A + w_s, minlength=NC * A).reshape(NC, A).astype(float)
M_c = M_cs.sum(1)


def tv_all(t_lab, diff=None):
    """TV_t and N_t for every tape, clone-stratified as documented above.

    If `diff` is an array it is filled with the SIGNED per-symbol excess
    p_ts - ref_ts, which is what the shape diagnostic reads."""
    tv = np.zeros(K); Nt = np.zeros(K)
    order = np.argsort(t_lab, kind="stable")
    ts, cs, ss = t_lab[order], w_c[order], w_s[order]
    bnd = np.searchsorted(ts, np.arange(K + 1))
    for t in range(K):
        lo, hi = bnd[t], bnd[t + 1]
        if hi - lo < 20:
            continue
        c, s = cs[lo:hi], ss[lo:hi]
        uct, cinv = np.unique(c, return_inverse=True)   # only clones using tape t
        nct = len(uct)
        n_cts = np.bincount(cinv * A + s, minlength=nct * A).reshape(nct, A).astype(float)
        n_ct = n_cts.sum(1)
        oth = M_cs[uct] - n_cts                  # other tapes, same clones
        den = np.maximum(M_c[uct] - n_ct, 1.0)
        act = n_ct.sum()
        p = n_cts.sum(0) / act
        ref = (n_ct[:, None] / act * (oth / den[:, None])).sum(0)
        tv[t] = 0.5 * np.abs(p - ref).sum(); Nt[t] = act
        if diff is not None:
            diff[t] = p - ref
    return tv, Nt


D_obs = np.zeros((K, A))
tv_obs, Nt = tv_all(w_t, diff=D_obs)
G_obs = float((Nt * tv_obs).sum())

# ---- null: permute the tape label within (clone, site) ---------------------
# ⚠ vectorised: a per-block python loop is ~15k blocks x 200 perms on Pre-TX.
# lexsort by (random, block) lays the labels of each block down in random order,
# and writing them back at the block-grouped positions IS a within-block shuffle.
blk = w_c * 8 + w_j
pos = np.argsort(blk, kind="stable")        # positions grouped by block
tv_null = np.zeros((NPERM, K)); G_null = np.zeros(NPERM)
Dsum = np.zeros((K, A)); Dsq = np.zeros((K, A)); Dbuf = np.zeros((K, A))
for b in range(NPERM):
    src = np.lexsort((rng.random(NW), blk))
    perm = np.empty_like(w_t); perm[pos] = w_t[src]
    Dbuf[:] = 0.0
    tv_null[b], _ = tv_all(perm, diff=Dbuf)
    Dsum += Dbuf; Dsq += Dbuf * Dbuf
    G_null[b] = float((Nt * tv_null[b]).sum())
    if b % 25 == 0:
        print(f"  perm {b}/{NPERM}  G_null={G_null[b]:,.0f}  (obs {G_obs:,.0f})", flush=True)

# ⚠ z is against the null's own mean and sd for that (tape, symbol) cell, so the
# collapse bias that the null DOES reproduce is already divided out; what is left
# is the part the null cannot make.
Dm = Dsum / NPERM
Dsd = np.sqrt(np.maximum(Dsq / NPERM - Dm * Dm, 0.0))
Z = np.where(Dsd > 0, (D_obs - Dm) / np.maximum(Dsd, 1e-12), 0.0)

scored = Nt > 0
nm, n95 = tv_null.mean(0), np.percentile(tv_null, 95, axis=0)
exc = tv_obs - nm
p_glob = (1 + int((G_null >= G_obs).sum())) / (NPERM + 1)
# concentration: does ONE symbol carry the tape's excess (co-integration), or is
# it spread over the frequent-symbol end (collapse bias)?
top_i = np.argmax(Z, axis=1)
top_z = Z[np.arange(K), top_i]
rest = Z.copy(); rest[np.arange(K), top_i] = 0.0
rms = np.sqrt((rest ** 2).sum(1) / max(A - 1, 1))
conc = np.where(rms > 0, top_z / np.maximum(rms, 1e-12), 0.0)
# bushiness: writes per distinct parent, the covariate collapse bias must track
bush = np.zeros(K)
for t_ in range(K):
    m_ = w_t == t_
    if m_.sum():
        bush[t_] = m_.sum() / len(np.unique(w_p[m_]))
sc = scored & (Nt >= 100)
r_bush = float(np.corrcoef(bush[sc], (tv_obs - nm)[sc])[0, 1]) if sc.sum() > 3 else float("nan")
r_conc = float(np.corrcoef(bush[sc], conc[sc])[0, 1]) if sc.sum() > 3 else float("nan")
xi_p = np.bincount(w_s, minlength=A) / NW
r_freq = float(np.corrcoef(xi_p, np.abs(Z[sc]).mean(0))[0, 1])
np.savez_compressed(RES / f"tape_cis_Z_{arm}.npz", Z=Z, D_obs=D_obs, Dm=Dm, Dsd=Dsd,
                    Nt=Nt, bush=bush, symbols=syms, tapes=tapes, xi=xi_p)

res = {
    "arm": arm, "poly": POLY, "n_writes": int(NW), "n_clones": int(NC),
    "shape": {
        "top_z_median": float(np.median(top_z[sc])),
        "top_z_max": float(top_z[sc].max()),
        "concentration_median": float(np.median(conc[sc])),
        "n_tapes_top_z_gt4": int((top_z[sc] > 4).sum()),
        "n_tapes_conc_gt3": int((conc[sc] > 3).sum()),
        "n_tapes_scored": int(sc.sum()),
        "corr_bushiness_excess": r_bush,
        "corr_bushiness_concentration": r_conc,
        "corr_symbol_freq_meanabsZ": r_freq,
    },
    "map": [{"tape": str(tapes[t]), "idx": int(t), "n_writes": int(Nt[t]),
             "symbol": str(syms[top_i[t]]), "z": float(top_z[t]),
             "concentration": float(conc[t]), "bushiness": float(bush[t])}
            for t in np.where(sc)[0]],
    "n_symbols": int(A), "nperm": NPERM,
    "tapes_scored": int(scored.sum()),
    "tv_obs_median": float(np.median(tv_obs[scored])),
    "tv_null_median": float(np.median(nm[scored])),
    "excess_median": float(np.median(exc[scored])),
    "excess_max": float(exc[scored].max()),
    "n_tapes_above_null_p95": int((tv_obs[scored] > n95[scored]).sum()),
    "G_obs": G_obs, "G_null_mean": float(G_null.mean()),
    "G_null_max": float(G_null.max()), "p_global": p_glob,
    "per_tape": [{"tape": str(tapes[t]), "idx": int(t), "n_writes": int(Nt[t]),
                  "tv": float(tv_obs[t]), "tv_null_mean": float(nm[t]),
                  "tv_null_p95": float(n95[t]), "excess": float(exc[t])}
                 for t in np.argsort(-np.where(scored, exc, -np.inf))[:int(scored.sum())]],
}
(RES / f"tape_cis_{arm}.json").write_text(json.dumps(res, indent=1))
print(f"{arm}: TV median obs {res['tv_obs_median']:.4f} vs null {res['tv_null_median']:.4f} "
      f"(excess {res['excess_median']:+.4f}); "
      f"{res['n_tapes_above_null_p95']}/{res['tapes_scored']} tapes above null p95", flush=True)
print(f"  global  G_obs {G_obs:,.0f}  null mean {G_null.mean():,.0f} "
      f"max {G_null.max():,.0f}   p = {p_glob:.4f}", flush=True)
sh = res["shape"]
print(f"  SHAPE   top z: median {sh['top_z_median']:.2f} max {sh['top_z_max']:.2f}; "
      f"{sh['n_tapes_top_z_gt4']}/{sh['n_tapes_scored']} tapes with top z > 4; "
      f"concentration median {sh['concentration_median']:.2f}, "
      f"{sh['n_tapes_conc_gt3']} tapes > 3", flush=True)
print(f"  COLLAPSE corr(bushiness, excess) = {sh['corr_bushiness_excess']:+.3f}; "
      f"corr(bushiness, concentration) = {sh['corr_bushiness_concentration']:+.3f}; "
      f"corr(symbol freq, mean|z|) = {sh['corr_symbol_freq_meanabsZ']:+.3f}", flush=True)
print("  top 8 tapes by z: " + "  ".join(
    f"{tapes[t]}->{syms[top_i[t]]}(z={top_z[t]:.1f},c={conc[t]:.1f})"
    for t in np.where(sc)[0][np.argsort(-top_z[sc])[:8]]), flush=True)
