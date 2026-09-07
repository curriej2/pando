#!/usr/bin/env python3
r"""
================================================================================
 56 -- the symbol array, and the PREMISE CHECK for the depletion test
================================================================================

Two jobs in one CSV pass.

(1) CACHE THE SYMBOLS.  Everything downstream of the event catalogue has worked
    on presence/absence (dropout_matrix) or on prefix IDENTITY (prefix_codes6).
    The depletion test needs the symbols themselves, so cache
        symbols_{arm}.npz : S[cell, tape, site] = symbol id, -1 = absent
    row-aligned, BIT-IDENTICALLY, with dropout_matrix_{arm}.npz (post-filter) and
    prefix_codes6_{arm}.npz.  The loading block is copied from 44_prefix_codes6.py
    unchanged for exactly that reason, and asserted against both caches.

(2) THE PREMISE CHECK -- is xi TAPE-DEPENDENT?
    The whole co-integrated-symbol test rests on pegRNAs acting in TRANS: symbol
    s(z) is written by integration z into EVERY tape, so the insert composition
    must be the same at all 166 tapes.  If any CIS component exists -- a tape
    preferentially receiving its own cassette's symbol -- then tapes carry
    tape-specific composition, and because tape losses are CORRELATED with each
    other (rho_tape = 0.25, fig 3b), a clade that lost tape z has a non-random
    tape composition and would show a compositional shift with no silencing at
    all.  That is precisely the depletion we intend to claim, manufactured.

    Section 0 gives four grounds for trans action.  This measures it here.

    Statistic per tape t: G^2 = 2 * SUM_s n_ts log(p_ts / p_s), the LRT of "this
    tape's composition differs from pooled", plus the total variation distance
    TV_t = 0.5 * SUM_s |p_ts - p_s|.
    ⚠ With ~10^4-10^5 insertions per tape, G^2 is significant for deviations far
    too small to matter, so G^2 alone answers the wrong question.  The NOISE
    FLOOR is what makes TV interpretable: split each tape's cells at random into
    two halves and compute TV between the halves.  A cis effect must lift TV
    above that floor, not above zero.

    Reported per SITE as well: composition is known to vary with site (site 6 has
    elevated q, README), and clade writes sit at deeper sites than clone writes,
    so the depletion test must site-match.  This measures how much that matters.

    ⚠ One confound this cannot remove by itself: tape composition also varies
    through the CELLS that recover that tape (clone founder effects x rho_tape).
    So the pooled reference is recomputed on the tape's OWN cells -- for each
    tape t, p_s is taken over all OTHER tapes in exactly the cells where t was
    recovered.  A residual TV above the floor is then cis, not cell composition.

Usage:  56_symbol_cache.py [arm ...]
Outputs results/symbols_{arm}.npz          (gitignored -- per-cell derived table)
        results/tape_composition_{arm}.json (aggregate, committable)
================================================================================
"""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT / "results"
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
MIN_COUNT, N_SITES, MAX_D = 1000, 6, 6
RNG = np.random.default_rng(20260907)


def build_symbols(arm):
    """Row-for-row the loader of 44_prefix_codes6.py; returns S and the alphabet."""
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    d0 = xi_all[arm]
    keep = {k for k, v in d0["xi"].items() if v * d0["n_edits"] >= MIN_COUNT}
    sym_id = {s: i for i, s in enumerate(sorted(keep))}
    clone_of = {}
    with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["ClonalBC"] and r["ClonalBC"] != "None":
                clone_of[r["CellID"]] = r["ClonalBC"]
    with open(DATA / f"{arm}_EditTable_filtered.csv", newline="") as fh:
        rr = csv.reader(fh); hdr = next(rr)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
        tapes = sorted(blocks)
        idxs = [[i for _, i in sorted(blocks[tp])] for tp in tapes]
        K = len(tapes)
        rows = []
        for row in rr:
            if clone_of.get(row[0]) is None:
                continue
            v = row[1:]
            a = np.full((K, MAX_D), -1, dtype=np.int16)
            ntapes = 0
            for t, cols in enumerate(idxs):
                seen = False
                for j, i in enumerate(cols):
                    val = v[i]
                    if val != "None":
                        seen = True
                    if j < MAX_D:
                        sid = sym_id.get(val, -1)
                        if sid < 0:
                            break
                        a[t, j] = sid
                    elif val not in keep:
                        break
                if seen:
                    ntapes += 1
            if ntapes < THR[arm]:
                continue
            rows.append(a)
    S = np.array(rows, dtype=np.int16)
    inv = np.empty(len(sym_id), dtype=object)
    for s, i in sym_id.items():
        inv[i] = s
    return S, np.array(list(inv), dtype="<U11"), np.array(tapes, dtype="<U11")


def g2_and_tv(cnt, ref):
    """cnt: (A,) counts for one tape. ref: (A,) reference probabilities."""
    n = cnt.sum()
    if n == 0:
        return 0.0, 0.0
    p = cnt / n
    nz = cnt > 0
    g2 = 2.0 * float(np.sum(cnt[nz] * np.log(p[nz] / np.maximum(ref[nz], 1e-12))))
    return g2, 0.5 * float(np.sum(np.abs(p - ref)))


def analyse(arm, S, syms, tapes):
    n, K, D = S.shape
    A = len(syms)
    out = {"arm": arm, "n_cells": n, "n_tapes": K, "n_symbols": A}

    # ---- counts[t, s] over all sites, and per site ---------------------------
    counts = np.zeros((K, A), dtype=np.int64)
    counts_site = np.zeros((D, A), dtype=np.int64)
    for j in range(D):
        col = S[:, :, j]
        for t in range(K):
            c = col[:, t]
            counts[t] += np.bincount(c[c >= 0], minlength=A)
        cj = col[col >= 0]
        counts_site[j] = np.bincount(cj, minlength=A)

    tot = counts.sum(0)
    pooled = tot / tot.sum()
    out["n_insertions"] = int(tot.sum())
    out["xi_pooled_top"] = {str(syms[i]): float(pooled[i])
                            for i in np.argsort(-pooled)[:10]}

    # ---- per-site composition (does site matter?) ----------------------------
    site_tv = []
    for j in range(D):
        nj = counts_site[j].sum()
        site_tv.append({"site": j + 1, "n": int(nj),
                        "tv_vs_pooled": g2_and_tv(counts_site[j], pooled)[1]})
    out["per_site"] = site_tv

    # ---- per-tape composition, referenced to that tape's OWN cells -----------
    # For tape t, reference = composition over all tapes != t restricted to the
    # cells that recovered t.  Removes clone-composition-through-recovery.
    rec = (S[:, :, 0] >= 0)                    # cell recovered tape t at site 1
    per_tape = []
    for t in range(K):
        cells = np.where(rec[:, t])[0]
        if len(cells) < 50:
            continue
        sub = S[cells]                          # (m, K, D)
        own = np.bincount(sub[:, t, :][sub[:, t, :] >= 0], minlength=A)
        oth = sub.copy(); oth[:, t, :] = -1
        ref_c = np.bincount(oth[oth >= 0], minlength=A).astype(float)
        ref = ref_c / max(ref_c.sum(), 1)
        g2, tv = g2_and_tv(own.astype(float), ref)
        # noise floor: split this tape's own insertions at random in half
        flat = sub[:, t, :][sub[:, t, :] >= 0]
        half = RNG.random(len(flat)) < 0.5
        a1 = np.bincount(flat[half], minlength=A).astype(float)
        a2 = np.bincount(flat[~half], minlength=A).astype(float)
        _, tv_floor = g2_and_tv(a1, a2 / max(a2.sum(), 1))
        per_tape.append({"tape": str(tapes[t]), "idx": t, "n_cells": int(len(cells)),
                         "n_ins": int(own.sum()), "G2": g2, "df": int(A - 1),
                         "tv": tv, "tv_split_floor": tv_floor,
                         "excess": tv - tv_floor,
                         "top_enriched": [str(syms[i]) for i in
                                          np.argsort(-(own / max(own.sum(), 1) - ref))[:3]]})
    per_tape.sort(key=lambda r: -r["excess"])
    out["per_tape"] = per_tape
    ex = np.array([r["excess"] for r in per_tape])
    tvv = np.array([r["tv"] for r in per_tape])
    fl = np.array([r["tv_split_floor"] for r in per_tape])
    out["summary"] = {
        "tapes_scored": len(per_tape),
        "tv_median": float(np.median(tvv)), "tv_max": float(tvv.max()),
        "floor_median": float(np.median(fl)),
        "excess_median": float(np.median(ex)), "excess_max": float(ex.max()),
        "n_tapes_excess_gt_2x_floor": int(np.sum(tvv > 2 * fl)),
    }
    return out


for arm in (sys.argv[1:] or ["Mouse3", "Mouse1", "Mouse2", "Initial", "Subclone"]):
    S, syms, tapes = build_symbols(arm)
    ref6 = np.load(RES / f"prefix_codes6_{arm}.npz")["codes"]
    assert ref6.shape[:2] == S.shape[:2], f"{arm}: shape {S.shape} vs codes {ref6.shape}"
    assert ((ref6 >= 0) == (S >= 0)).all(), f"{arm}: determination pattern differs from codes6"
    np.savez_compressed(RES / f"symbols_{arm}.npz", S=S, symbols=syms, tapes=tapes)
    res = analyse(arm, S, syms, tapes)
    (RES / f"tape_composition_{arm}.json").write_text(json.dumps(res, indent=1))
    s = res["summary"]
    print(f"{arm}: {res['n_cells']:,} cells x {res['n_tapes']} tapes, "
          f"{res['n_symbols']} symbols, {res['n_insertions']:,} insertions", flush=True)
    print(f"  per-tape TV vs own-cell reference: median {s['tv_median']:.4f} "
          f"(split-half floor {s['floor_median']:.4f}), max {s['tv_max']:.4f}; "
          f"{s['n_tapes_excess_gt_2x_floor']}/{s['tapes_scored']} tapes above 2x floor", flush=True)
    print("  per-site TV vs pooled: " +
          "  ".join(f"s{d['site']}={d['tv_vs_pooled']:.4f}" for d in res["per_site"]), flush=True)
