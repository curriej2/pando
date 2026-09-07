#!/usr/bin/env python3
r"""
================================================================================
 59 -- DOES THE CO-INTEGRATED SYMBOL VANISH WHEN A TAPE IS SILENCED?
       Hand-inspectable examples inside large clades.        (README direction 2)
================================================================================

THE FIRST ORTHOGONAL TEST OF THE MECHANISM.  Everything in the event catalogue
infers silencing from MISSINGNESS, which is also what technical dropout looks
like.  The cassette PB-U6-pegRNA-NNNNGGA-EF1a-mRFP-TAPE-TargetBC carries a pegRNA
and a tape on ONE integration, and pegRNAs act in TRANS, so integration z's
symbol s(z) is written into EVERY tape in the cell.  If silencing kills the whole
locus then

    losing tape z from the readout  =>  symbol s(z) stops appearing AT ALL OTHER
                                        TAPES in the same cells.

That is a different channel -- the symbols written into other tapes -- which
transcript capture cannot reach.

⚠ The map z -> s(z) is UNKNOWN (TargetBC is 10-nt, NNNN is 4-nt, never linked).
So do not test a known pair: ask of each event WHICH symbol is depleted, and let
the map's own structure be the evidence -- five-way cross-arm concordance (all
166 TargetBCs are identical across arms), near-injectivity at the collision rate
predicted for 166 draws from 256, and corr(beta_z, xi_s(z)) > 0.

--------------------------------------------------------------------------------
 DESIGN -- the four decisions, taken 2026-09-07
--------------------------------------------------------------------------------
 UNIT      independent post-MRCA writes (57_writes).  ⚠ NOT raw entries: a symbol
           in the clade founder recurs in every descendant, and the pseudo-
           replication that creates is what invalidated 56's first premise check.
 LAYER     sub-clone clades in the two large clones.  Step 0 measured
           W = 214,076 writes inside Subclone clone 6's depth-1 clade and 28,933
           inside Mouse 2 clone 76's -- so ~1,238 and ~216 expected copies of a
           median-xi symbol.  This layer, not the clone-wide one, because the
           rest of the clone is then an INTERNAL control no clone-level or batch
           effect can reach.  (The README expected clone-wide to be needed for
           power; measurement says otherwise.)
 CONTRAST  rest of clone, SITE-MATCHED.  Site matters: TV vs pooled at site 6 is
           ~2x sites 1-5, and clade writes sit at deeper sites than clone writes,
           so unmatched, depth alone would masquerade as a compositional shift.
           ⚠ Cells that lost the ANCHOR or never reached `depth` are put in
           NEITHER group -- the only evidence they lie outside the clade is that
           we could not read their anchor (script 10's D-vs-S distinction).
 EXCLUDED  the anchor tape (the clade is DEFINED by its symbols) and the lost
           tape z (the cross-tape rule that keeps Lambda honest).

--------------------------------------------------------------------------------
 THE STATISTIC -- a site-stratified common odds ratio, reported in nats
--------------------------------------------------------------------------------
For symbol s, site j: inside a_j of A_j writes are s; outside b_j of B_j are.

    inside odds = exp(theta_s) * outside odds,   one theta_s, per-site baselines mu_j

    Lambda_s = l(theta_hat_s) - l(0)   in NATS, signed by sign(theta_hat_s).

Two reasons for the odds ratio rather than a raw frequency difference:
 * RENORMALISATION IS DIVIDED OUT.  Removing s inflates every other symbol by
   xi_s/(1-xi_s) ~ 1%; in the odds of s against the rest that cancels to second
   order, so the other ~100 symbols do not all drift in one direction.
 * Under complete removal theta -> -inf, under partial removal of a symbol
   supplied by j integrations theta -> log(1 - 1/j).  ⚠ So MAGNITUDE IS NOT THE
   SIGNAL: 166 integrations map onto ~100-106 observed symbols (mean multiplicity
   ~1.6) and xi spans 570x, so how deep the depletion goes says little.  The
   RECOVERED IDENTITY s(z), and its reproducibility, is the evidence.

⚠⚠ THE STATISTIC DOES NOT DEPEND ON z EXCEPT THROUGH THE CLADE.  A clade shows
whatever composition shift it has; if it lost several tapes we should see several
depleted symbols.  So a "control tape" is not a control here.  The real test is
whether the depleted SET matches the lost SET, consistently, across events -- and
the honest per-example control is a clade that lost nothing, which must be flat.

--------------------------------------------------------------------------------
 NULL
--------------------------------------------------------------------------------
The project's rule: permute the label being tested, blocked by the label above
it.  Whole cell ROWS are shuffled WITHIN CLONE, so clade membership is destroyed
while every cell's own tape and symbol content travels with it.  The writes are
then RE-EXTRACTED, because which writes are post-MRCA depends on who is in the
clade -- that is the expensive part and it is not optional.
Reported: the null distribution of min_s theta_s and of max_s Lambda_s, i.e. of
the extreme over ~100 symbols, so the multiplicity is calibrated rather than
Bonferroni-corrected.

Usage: 59_depletion.py <arm> [--top N] [--nperm 200] [--poly codes|codes_or_none]
       [--control]   also score size-matched clades in the same clone with NO
                     called loss, which must come out flat
Outputs results/depletion_{arm}.json           (aggregate, committable)
        results/depletion_{arm}.tsv.gz         (per symbol per event; gitignored)
================================================================================
"""
import gzip, json, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
_w = __import__("57_writes")
load_arm, extract_writes, clade_of, _events = (
    _w.load_arm, _w.extract_writes, _w.clade_of, _w._events)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
NSITE = 6

arm = sys.argv[1]
TOP = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 4
NPERM = int(sys.argv[sys.argv.index("--nperm") + 1]) if "--nperm" in sys.argv else 200
POLY = sys.argv[sys.argv.index("--poly") + 1] if "--poly" in sys.argv else "codes"
CONTROL = "--control" in sys.argv
rng = np.random.default_rng(20260907)

Y, clone, codes, S, syms, tapes = load_arm(arm)
# ⚠ conforming symbols only: junk is a per-amplicon parsing artefact and so is
# intrinsically tape-specific (57_writes.conforming).
_cmask, _remap, syms = _w.conforming(syms)
A = len(syms)
print(f"{arm}: alphabet {A} conforming of {len(_cmask)} symbols", flush=True)
THETA = np.linspace(-8.0, 8.0, 321)


def tabulate(t, j, s):
    """(J, A) counts of symbol s at site j, over the writes handed in."""
    n = np.bincount(j * A + s, minlength=NSITE * A)
    return n.reshape(NSITE, A).astype(float)


def fit_theta(a, b):
    """Site-stratified common log-OR for every symbol at once.

    a, b : (J, A) inside / outside counts.  Returns theta_hat (A,), Lam (A,).
    For each theta on the grid the per-site baselines mu_j are profiled out by
    Newton on   A_j*sigma(mu+theta) + B_j*sigma(mu) = a_j + b_j,  vectorised
    over (site, symbol, theta).  Lambda is then max_theta l - l(0), in nats.
    """
    A_j, B_j = a.sum(1), b.sum(1)                       # (J,)
    Aj = A_j[:, None, None]; Bj = B_j[:, None, None]
    tot = (a + b)[:, :, None]                           # (J, A, 1)
    th = THETA[None, None, :]
    mu = np.zeros((NSITE, A, len(THETA)))
    for _ in range(60):
        p1 = 1.0 / (1.0 + np.exp(-(mu + th)))
        p0 = 1.0 / (1.0 + np.exp(-mu))
        g = Aj * p1 + Bj * p0 - tot
        h = Aj * p1 * (1 - p1) + Bj * p0 * (1 - p0)
        mu -= np.clip(g / np.maximum(h, 1e-12), -2.0, 2.0)
    p1 = np.clip(1.0 / (1.0 + np.exp(-(mu + th))), 1e-15, 1 - 1e-15)
    p0 = np.clip(1.0 / (1.0 + np.exp(-mu)), 1e-15, 1 - 1e-15)
    ll = (a[:, :, None] * np.log(p1) + (Aj - a[:, :, None]) * np.log1p(-p1)
          + b[:, :, None] * np.log(p0) + (Bj - b[:, :, None]) * np.log1p(-p0)).sum(0)
    k = np.argmax(ll, axis=1)
    i0 = int(np.argmin(np.abs(THETA)))
    return THETA[k], ll[np.arange(A), k] - ll[:, i0]


def score(ins, out, drop):
    ti, ji, si = _w.filter_writes(
        *extract_writes(codes[ins], S[ins], Y[ins], drop_tapes=drop, poly=POLY), _remap)
    to, jo, so = _w.filter_writes(
        *extract_writes(codes[out], S[out], Y[out], drop_tapes=drop, poly=POLY), _remap)
    if len(ti) < 50 or len(to) < 50:
        return None
    a, b = tabulate(ti, ji, si), tabulate(to, jo, so)
    th, lam = fit_theta(a, b)
    return a, b, th, lam


rows, summary = [], []
cvals = _w.clone_values(clone)          # events `clone` column is an INDEX, not the id
ev = sorted(_events(arm), key=lambda r: -int(r["clade_cells"]))[:TOP]
for r in ev:
    cl, an = int(cvals[int(r["clone"])]), int(r["anchor_tape"])
    d, tp, nc = int(r["depth"]), int(r["tape"]), int(r["clade_cells"])
    t0 = time.time()
    ins, out = clade_of(codes, clone, cl, an, d, want_cells=nc)
    res = score(ins, out, (an, tp))
    if res is None:
        print(f"  clone {cl} anchor {an} d{d} tape {tp}: too few writes", flush=True)
        continue
    a, b, th, lam = res
    signed = np.where(th < 0, lam, -lam)                # + = depletion evidence
    ordr = np.argsort(-signed)
    # --- null: re-permute clade membership within the clone, re-extract -------
    g = np.where(clone == cl)[0]
    nin = int(ins.sum())
    null_min_th, null_max_lam = [], []
    for bi in range(NPERM):
        pick = rng.permutation(len(g))[:nin]
        m1 = np.zeros(len(clone), dtype=bool); m1[g[pick]] = True
        m0 = np.zeros(len(clone), dtype=bool); m0[g] = True; m0 &= ~m1
        rp = score(m1, m0, (an, tp))
        if rp is None:
            continue
        _, _, th_p, lam_p = rp
        sp = np.where(th_p < 0, lam_p, -lam_p)
        null_min_th.append(float(th_p.min())); null_max_lam.append(float(sp.max()))
    nml = np.array(null_max_lam)
    top = ordr[0]
    p_ext = (1 + int((nml >= signed[top]).sum())) / (len(nml) + 1)
    print(f"  clone {cl} anchor {an:3d} d{d} tape {tp:3d}  n_in={nin:5d} "
          f"W_in={int(a.sum()):7,d} W_out={int(b.sum()):7,d}  "
          f"top={syms[top]} theta={th[top]:+.3f} Lambda={signed[top]:8.1f} "
          f"| null max Lambda mean {nml.mean():6.1f} max {nml.max():7.1f} "
          f"p={p_ext:.4f}   [{time.time()-t0:.0f}s]", flush=True)
    print("     next 4: " + "  ".join(
        f"{syms[i]} {th[i]:+.2f}/{signed[i]:.0f}" for i in ordr[1:5]), flush=True)
    summary.append({
        "clone": cl, "anchor": an, "depth": d, "tape": tp, "n_in": nin,
        "n_out": int(out.sum()), "W_in": int(a.sum()), "W_out": int(b.sum()),
        "lambda_event_nats": float(r["lambda_nats"]),
        "top_symbol": str(syms[top]), "top_theta": float(th[top]),
        "top_lambda": float(signed[top]), "p_permutation": p_ext,
        "null_max_lambda_mean": float(nml.mean()), "null_max_lambda_max": float(nml.max()),
        "nperm_used": int(len(nml)),
        "ranked": [{"symbol": str(syms[i]), "theta": float(th[i]),
                    "lambda": float(signed[i]), "n_in": int(a[:, i].sum()),
                    "n_out": int(b[:, i].sum())} for i in ordr[:8]],
    })
    for i in range(A):
        rows.append((cl, an, d, tp, str(syms[i]), th[i], signed[i],
                     int(a[:, i].sum()), int(b[:, i].sum())))

out_json = {"arm": arm, "poly": POLY, "nperm": NPERM, "n_symbols": int(A),
            "events": summary}
(RES / f"depletion_{arm}.json").write_text(json.dumps(out_json, indent=1))
with gzip.open(RES / f"depletion_{arm}.tsv.gz", "wt") as fh:
    fh.write("clone\tanchor\tdepth\ttape\tsymbol\ttheta\tlambda_signed\tn_in\tn_out\n")
    for r_ in rows:
        fh.write("\t".join(str(x) for x in r_) + "\n")
print(f"{arm}: {len(summary)} events scored -> depletion_{arm}.json", flush=True)
