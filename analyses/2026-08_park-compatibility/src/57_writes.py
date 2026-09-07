#!/usr/bin/env python3
r"""
================================================================================
 57 -- INDEPENDENT POST-MRCA WRITES: the sample unit for the depletion test
================================================================================

⚠⚠ WHY THIS MODULE EXISTS, AND THE MISTAKE IT FIXES.  The first version of the
premise check (56) counted raw (cell, tape, site) entries and took its noise
floor from a random split of those entries.  That is pseudoreplication: a symbol
written at tape y in a clone founder reappears in EVERY descendant, so the
effective sample size is a small fraction of the entry count and the floor was
far too low.  It is the same lesson as "pooled rho weights clones by n_C(n_C-1)"
(README, estimator rule 1) in a different guise.  Every composition statistic in
this direction of work must be computed on DEDUPLICATED WRITES.

--------------------------------------------------------------------------------
 THE UNIT
--------------------------------------------------------------------------------
A tape is an append-only record, so a symbol already present when an integration
was silenced stays there.  Only insertions laid down AFTER the loss can show a
depletion.  For a cell group G (a clone, or a clade), tape y and site j:

    parent  = the depth-(j-1) prefix code (for j = 1, a single root)
    children = the depth-j prefix codes of G's cells carrying that parent

  A parent is POLYMORPHIC in G at site j iff its children take >= 2 distinct
  values.  Then the MRCA of G cannot have had site j written -- if it had, every
  descendant would carry that one symbol -- so every distinct child code under it
  is one write that happened INSIDE G, i.e. after any loss on G's stem.

Two properties make this the right unit:

 1. IDENTITY BY DESCENT IS DEDUPLICATED EXACTLY.  Cells sharing a write share the
    child code, so they contribute it once.  No effective-sample-size guesswork.
 2. THE CLADE'S OWN ANCESTRY IS EXCLUDED.  Ancestral content is by construction
    monomorphic in G and so contributes nothing.

⚠ Cost of the dedup, stated: two INDEPENDENT writes of the same symbol under the
same parent collapse to one child code, so frequent symbols are deflated.  This
is the homoplasy collapse the m-estimation work (README "Is m recoverable from
s?") already characterises.  It biases inside and outside groups identically, so
a CONTRAST survives it; an absolute xi estimate would not.

--------------------------------------------------------------------------------
 POLYMORPHISM CONVENTION -- and why the default is the conservative one
--------------------------------------------------------------------------------
"NONE" (site unedited) is also evidence the MRCA had site j blank, so counting it
as a child value finds more polymorphic parents.  But a tape truncated by a JUNK
symbol reads as NONE from that point too (README: 6.88% of (cell, tape) instances
contain junk), and that is missing data, not an unedited site.

    --poly codes        (DEFAULT) require >= 2 distinct REAL child codes.
                        Immune to junk truncation.  Fewer writes.
    --poly codes_or_none  count NONE as a child value.  More writes, and admits
                        junk truncation as false polymorphism.

Both are reported by the CLI so the gap is visible rather than assumed.

--------------------------------------------------------------------------------
 CLI -- the Step 0 power table
--------------------------------------------------------------------------------
    57_writes.py <arm> [--top N] [--clades] [--poly codes|codes_or_none]

reports W per clone (and per called clade with --clades) and the expected copies
of a median-xi symbol, which is the power calculation for the depletion test:
seeing a symbol vanish needs W * xi_s to be comfortably above a few nats.
================================================================================
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}


def load_arm(arm):
    """Y (recovered), clone, codes6, S, symbols, tapes -- all row-aligned."""
    z0 = np.load(RES / f"dropout_matrix_{arm}.npz")
    Y, clone = z0["recovered"], z0["clone"]
    sel = (Y.sum(1) >= THR[arm]) & (clone >= 0)
    Y, clone = Y[sel], clone[sel]
    codes = np.load(RES / f"prefix_codes6_{arm}.npz")["codes"]
    zs = np.load(RES / f"symbols_{arm}.npz")
    S, syms, tapes = zs["S"], zs["symbols"], zs["tapes"]
    assert codes.shape[0] == Y.shape[0] == S.shape[0], "row misalignment"
    assert ((codes >= 0) == (S >= 0)).all(), "codes/symbol determination differ"
    return Y, clone, codes, S, syms, tapes


def extract_writes(codes_g, S_g, rec_g, drop_tapes=(), poly="codes", parents=False):
    """Independent post-MRCA writes for one cell group.

    codes_g (m,K,D) int32, S_g (m,K,D) int16, rec_g (m,K) bool.
    Returns (tape, site, symbol) int arrays, one entry per write; with
    parents=True also the parent key of each write, which is what the COLLAPSE
    diagnostic needs -- writes per distinct parent measures how bushy a tape's
    trie is, and the dedup collapses hardest exactly where it is bushiest.
    `drop_tapes` are excluded from the pool -- ALWAYS pass the anchor tape (the
    clade is defined by its symbols) and the tape under test (cross-tape rule).
    """
    m, K, D = codes_g.shape
    keep = np.ones(K, dtype=bool)
    for t in drop_tapes:
        keep[t] = False
    ti = np.where(keep)[0]
    if m < 2 or len(ti) == 0:
        e = np.zeros(0, dtype=np.int64)
        return e, e, e
    C = int(codes_g.max()) + 2                       # child alphabet incl. NONE
    o_t, o_j, o_s, o_p = [], [], [], []
    for j in range(D):
        ch = codes_g[:, ti, j].astype(np.int64)
        sy = S_g[:, ti, j].astype(np.int64)
        if j == 0:
            par = np.where(rec_g[:, ti], 0, -1).astype(np.int64)
        else:
            par = codes_g[:, ti, j - 1].astype(np.int64)
        ok = par >= 0
        if not ok.any():
            continue
        tp = np.broadcast_to(ti[None, :], (m, len(ti)))
        # one key per (tape, parent, child-or-NONE); NONE encoded as 0
        key = (tp[ok].astype(np.int64) * (int(par.max()) + 2) + par[ok]) * C + (ch[ok] + 1)
        uk, first = np.unique(key, return_index=True)
        pkey, cval = uk // C, uk % C
        # distinct child values per (tape, parent), then which parents are poly
        upk, inv = np.unique(pkey, return_inverse=True)
        if poly == "codes":
            nval = np.bincount(inv, weights=(cval > 0).astype(float), minlength=len(upk))
        elif poly == "codes_or_none":
            nval = np.bincount(inv, minlength=len(upk))
        else:
            raise ValueError(poly)
        good = (nval[inv] >= 2) & (cval > 0)         # real writes under poly parents
        if not good.any():
            continue
        o_t.append(tp[ok][first[good]])
        o_j.append(np.full(int(good.sum()), j, dtype=np.int64))
        o_s.append(sy[ok][first[good]])
        if parents:
            o_p.append(pkey[good] * 8 + j)           # parent identity, site-tagged
    if not o_t:
        e = np.zeros(0, dtype=np.int64)
        return (e, e, e, e) if parents else (e, e, e)
    out = (np.concatenate(o_t), np.concatenate(o_j), np.concatenate(o_s))
    return out + (np.concatenate(o_p),) if parents else out


def clade_of(codes, clone, cl, anchor, depth, want_cells=None):
    """The cell mask of one catalogue clade, plus the 'determined not in it' mask.

    ⚠ `cl` is the RAW clone value, not the events table's `clone` column -- that
    column holds the g_clone INDEX into np.unique(clone), and the two differ
    wherever a clone id was dropped by the tape filter (Mouse 3 index 110 is raw
    clone 111, a 210-cell clone, while raw clone 110 has 2 cells).  Callers
    reading events_*.tsv.gz must map through `clone_values(clone)` first.

    Cells excluded from the clade because they lost the anchor or never reached
    `depth` are NOT put in the comparison group: the only evidence they lie
    outside is that we could not read their anchor.  That is the same D-vs-S
    distinction script 10 makes -- missing must not be read as absent.
    """
    g = clone == cl
    code = codes[:, anchor, depth - 1]
    det = g & (code >= 0)
    vals, cnt = np.unique(code[det], return_counts=True)
    if want_cells is not None:
        hit = vals[cnt == want_cells]
        assert len(hit) == 1, f"clade size {want_cells} matches {len(hit)} codes"
        best = hit[0]
    else:
        best = vals[np.argmax(cnt)]
    inside = det & (code == best)
    return inside, det & ~inside


def conforming(syms):
    """Mask + dense remap for DESIGN-CONFORMING insert symbols (NNNNGGA).

    ⚠⚠ The pegRNA insert is a 4-nt NNNN followed by the fixed GGA.  Of the ~104
    symbols that clear the count threshold, 6 do not conform (Mouse 2: TGGGA,
    GATGGA, CACGGA, GCGCGG, GTCG and a 91-nt read-through), carrying 2.2% of
    edits -- and measured 2026-09-07 they DOMINATE the per-tape signal: Mouse 2's
    strongest tape effect was TGGGA at z = 67, Subclone's GAATGGATGAT at z = 99.
    That is expected once stated: junk is a PARSING artefact of a particular
    amplicon, so it is intrinsically tape-specific, and it says nothing about the
    pegRNA writing pool.

    ⚠ This does NOT contradict the README correction "junk values ARE heritable
    characters".  Junk is a perfectly good LINEAGE character -- it is inherited.
    It is not evidence about insert COMPOSITION, which is the only thing the
    depletion test reads.  The two uses need different alphabets.

    Returns (mask over the old ids, remap old id -> dense new id or -1, new names).
    """
    m = np.array([len(s) == 7 and s.endswith("GGA")
                  and set(s[:4]) <= set("ACGT") for s in syms], dtype=bool)
    remap = np.full(len(syms), -1, dtype=np.int64)
    remap[m] = np.arange(int(m.sum()))
    return m, remap, np.asarray(syms)[m]


def filter_writes(t_, j_, s_, remap):
    """Drop writes whose symbol is non-conforming and renumber the rest."""
    ns = remap[s_]
    k = ns >= 0
    return t_[k], j_[k], ns[k]


def clone_values(clone):
    """np.unique(clone) -- index i of the events table's `clone` column -> raw id."""
    return np.unique(clone)


def _events(arm, tag="_lam2_b2"):
    rows = []
    with gzip.open(RES / f"events_{arm}{tag}.tsv.gz", "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            rows.append(dict(zip(hdr, ln.rstrip("\n").split("\t"))))
    return rows


if __name__ == "__main__":
    arm = sys.argv[1]
    TOP = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 5
    POLY = sys.argv[sys.argv.index("--poly") + 1] if "--poly" in sys.argv else None
    Y, clone, codes, S, syms, tapes = load_arm(arm)
    xi = json.loads((RES / "xi_vectors.json").read_text())[arm]["xi"]
    xmed = float(np.median(list(xi.values())))
    print(f"{arm}: {len(clone):,} cells, {len(syms)} symbols, median xi {xmed:.5f}", flush=True)
    modes = [POLY] if POLY else ["codes", "codes_or_none"]

    names, sizes = np.unique(clone, return_counts=True)
    for i in np.argsort(-sizes)[:TOP]:
        g = clone == names[i]
        line = f"  clone {names[i]:>6} n={sizes[i]:5d}"
        for md in modes:
            t, j, s = extract_writes(codes[g], S[g], Y[g], poly=md)
            line += f" | {md}: W={len(t):8,d} E[med xi]={len(t)*xmed:7.0f}"
        print(line, flush=True)

    if "--clades" in sys.argv:
        cvals = clone_values(clone)
        ev = sorted(_events(arm), key=lambda r: -int(r["clade_cells"]))[:TOP]
        for r in ev:
            cl, an = int(cvals[int(r["clone"])]), int(r["anchor_tape"])
            d, tp, nc = int(r["depth"]), int(r["tape"]), int(r["clade_cells"])
            ins, out = clade_of(codes, clone, cl, an, d, want_cells=nc)
            for lab, msk in (("clade", ins), ("out", out)):
                line = (f"  clone {cl:>6} anchor {an:3d} d{d} tape {tp:3d} "
                        f"{lab:5s} n={int(msk.sum()):5d}")
                for md in modes:
                    t, j, s = extract_writes(codes[msk], S[msk], Y[msk],
                                             drop_tapes=(an, tp), poly=md)
                    line += f" | {md}: W={len(t):8,d}"
                print(line, flush=True)
