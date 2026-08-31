#!/usr/bin/env python3
r"""
================================================================================
 D.4b step 4-5: the cross-tape compatibility test
================================================================================

Two characters are COMPATIBLE iff their clades are nested or disjoint.  Partial
overlap means no single tree carries both.  Counting three witnesses:

    n11 = |S1 & S2|      n10 = |S1 \ S2|      n01 = |S2 \ S1|

    INCOMPATIBLE  iff  all three are non-zero.

Only CROSS-TAPE pairs are tested: within a tape the prefixes form a trie and are
nested by construction (D.4 fact 1), verified on real data by script 10.

--------------------------------------------------------------------------------
 HOW IT IS COMPUTED -- cross-tabulation, not pairwise set operations
--------------------------------------------------------------------------------
At a fixed level L, the characters on tape z PARTITION the cells: a cell has
exactly one level-L prefix, or none.  So represent them as a LABEL VECTOR (one
entry per cell) rather than as sets.  Then for a tape pair (z,z') and a level pair
(L,L'), one cross-tabulation of the two label vectors gives every character pair
at those levels at once:

        n11 = tab[p][p']        n10 = rowsum[p] - n11        n01 = colsum[p'] - n11

because the row IS one character's clade and the column IS the other's, so the
interior entry is exactly their intersection and the margins supply the rest.

That replaces C-choose-2 pairwise tests (19.9e9 over the whole dataset, 53% of them
in one clone) with 166-choose-2 x 6 x 6 ~ 493k table builds per clone, each a single
pass over that clone's cells.  Disjoint pairs are the ZERO entries: never
enumerated, counted by subtraction.  The identity is verified against brute force
in __main__ before anything else runs.

--------------------------------------------------------------------------------
 THE THREE-VALUED PART (D.4b Procedure step 4)
--------------------------------------------------------------------------------
A label is one of: a real prefix code / DETERMINED-ABSENT / UNDETERMINED.
DETERMINED-ABSENT must be its own category, not dropped: a cell carrying p on tape
z whose tape z' was read but stopped short is genuinely in S1\S2 and must count
toward n10.

    missing-excluded  : build the table only over cells determined on BOTH tapes
    missing-as-absent : fold UNDETERMINED into DETERMINED-ABSENT

Both are run; the spread between them is the dropout sensitivity.  D.4b is explicit
that a single number is not a result.

Cells whose level-L prefix is unique (or shared by fewer than MIN_CARRIERS cells)
carry no stored character, and are folded into DETERMINED-ABSENT.  That is exact:
such a cell is genuinely determined NOT to carry any character we evaluate.

Outputs results/compatibility_{table}.json and per-character degrees.
================================================================================
"""
import json, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
MIN_CARRIERS = 3
ABSENT, UNDET = -1, -2


def load(tbl):
    z = np.load(RES / f"characters_{tbl}.npz", allow_pickle=True)
    key = z["dmask_key"]; lens = z["dmask_lens"]; flat = z["dmask"]
    off = np.concatenate([[0], np.cumsum(lens)])
    D = {(int(a), int(b), int(c)): np.unpackbits(flat[off[i]:off[i+1]]).astype(bool)
         for i, (a, b, c) in enumerate(key)}
    return z, D


def clone_labels(z, D, ci, n_cells, min_carriers):
    """labels[(tape, level)] -> int array over cells; and char id -> (tape, level, size)."""
    sel = np.flatnonzero((z["ch_clone"] == ci) & (z["ch_len"] >= min_carriers))
    lab = {}
    meta = {}
    for i in sel:
        t, L = int(z["ch_tape"][i]), int(z["ch_level"][i])
        k = (t, L)
        if k not in lab:
            d = D.get((ci, t, L))
            base = np.full(n_cells, UNDET, np.int32)
            if d is not None:
                base[d[:n_cells]] = ABSENT
            lab[k] = base
        car = z["carriers"][z["ch_start"][i]:z["ch_start"][i] + z["ch_len"][i]]
        lab[k][car] = i
        meta[i] = (t, L, int(z["ch_len"][i]))
    return lab, meta


def crosstab_pairs(a, b, excluded):
    """Yield (id_a, id_b, n11, n10, n01) for co-occurring pairs; plus (#pairs, #disjoint)."""
    if excluded:
        m = (a != UNDET) & (b != UNDET)
        a, b = a[m], b[m]
    else:
        a = np.where(a == UNDET, ABSENT, a)
        b = np.where(b == UNDET, ABSENT, b)
    ra = np.unique(a[a >= 0]); rb = np.unique(b[b >= 0])
    if ra.size == 0 or rb.size == 0:
        return [], 0, 0
    ia = np.searchsorted(ra, a); ib = np.searchsorted(b * 0 + rb.searchsorted(b, side="left") * 0 + rb, b) if False else np.searchsorted(rb, b)
    va = (a >= 0) & (ra[np.clip(ia, 0, ra.size - 1)] == a)
    vb = (b >= 0) & (rb[np.clip(ib, 0, rb.size - 1)] == b)
    rowsum = np.bincount(ia[va], minlength=ra.size)
    colsum = np.bincount(ib[vb], minlength=rb.size)
    both = va & vb
    tab = np.bincount(ia[both] * rb.size + ib[both], minlength=ra.size * rb.size)
    nz = np.flatnonzero(tab)
    out = []
    for k in nz:
        i, j = divmod(int(k), rb.size)
        n11 = int(tab[k]); n10 = int(rowsum[i]) - n11; n01 = int(colsum[j]) - n11
        out.append((int(ra[i]), int(rb[j]), n11, n10, n01))
    return out, ra.size * rb.size, ra.size * rb.size - nz.size


def run_clone(z, D, ci, n_cells, conv_excluded):
    lab, meta = clone_labels(z, D, ci, n_cells, MIN_CARRIERS)
    keys = sorted(lab)
    by_tape = defaultdict(list)
    for t, L in keys: by_tape[t].append(L)
    tapes = sorted(by_tape)
    tested = incompat = disjoint = 0
    deg = defaultdict(int)
    for x in range(len(tapes)):
        for y in range(x + 1, len(tapes)):
            t1, t2 = tapes[x], tapes[y]
            for L1 in by_tape[t1]:
                a = lab[(t1, L1)]
                for L2 in by_tape[t2]:
                    pairs, npair, ndis = crosstab_pairs(a, lab[(t2, L2)], conv_excluded)
                    tested += npair; disjoint += ndis
                    for i, j, n11, n10, n01 in pairs:
                        if n10 > 0 and n01 > 0:
                            incompat += 1; deg[i] += 1; deg[j] += 1
    return dict(characters=len(meta), pairs=tested, disjoint=disjoint,
                incompatible=incompat,
                compatible=tested - incompat), deg


def brute_force(z, D, ci, n_cells, conv_excluded):
    sel = np.flatnonzero((z["ch_clone"] == ci) & (z["ch_len"] >= MIN_CARRIERS))
    S, Dm, TL = {}, {}, {}
    for i in sel:
        t, L = int(z["ch_tape"][i]), int(z["ch_level"][i])
        s = np.zeros(n_cells, bool)
        s[z["carriers"][z["ch_start"][i]:z["ch_start"][i]+z["ch_len"][i]]] = True
        S[i] = s; TL[i] = (t, L)
        d = D.get((ci, t, L)); Dm[i] = d[:n_cells] if d is not None else np.zeros(n_cells, bool)
    ids = list(S); inc = tot = 0
    for x in range(len(ids)):
        for y in range(x+1, len(ids)):
            i, j = ids[x], ids[y]
            if TL[i][0] == TL[j][0]: continue
            tot += 1
            e = (Dm[i] & Dm[j]) if conv_excluded else np.ones(n_cells, bool)
            s1, s2 = S[i] & e, S[j] & e
            if (s1 & s2).any() and (s1 & ~s2).any() and (s2 & ~s1).any(): inc += 1
    return tot, inc


if __name__ == "__main__":
    tbl = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"
    z, D = load(tbl)
    sizes = z["clone_sizes"]
    nchar = np.bincount(z["ch_clone"][z["ch_len"] >= MIN_CARRIERS], minlength=len(sizes))

    # ---- validate the cross-tab engine against brute force on a small real clone
    cand = [i for i in range(len(sizes)) if 30 <= nchar[i] <= 400]
    print(f"{tbl}: validating cross-tab against brute force on {len(cand[:3])} real clones")
    for ci in cand[:3]:
        for exc in (False, True):
            fast, _ = run_clone(z, D, ci, int(sizes[ci]), exc)
            tot, inc = brute_force(z, D, ci, int(sizes[ci]), exc)
            tag = "excluded" if exc else "as-absent"
            ok = (fast["pairs"] == tot) and (fast["incompatible"] == inc)
            print(f"  clone {ci:>4} ({sizes[ci]:>4} cells, {nchar[ci]:>4} chars) {tag:>10}: "
                  f"pairs {fast['pairs']:>7,}/{tot:>7,}  incompat {fast['incompatible']:>6,}/{inc:>6,}"
                  f"  {'PASS' if ok else 'FAIL'}")
            assert ok, "cross-tab does not reproduce brute force"
    print("  cross-tab engine verified\n")

    # ---- full run, clones in increasing character count
    order = sorted([i for i in range(len(sizes)) if nchar[i] > 0], key=lambda i: nchar[i])
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 1800.0
    res, degs, t0 = [], {}, time.time()
    for ci in order:
        if time.time() - t0 > budget:
            print(f"  [time budget {budget:.0f}s reached; {len(res)}/{len(order)} clones done]")
            break
        row = dict(clone=int(ci), cells=int(sizes[ci]), characters=int(nchar[ci]))
        for exc, tag in ((False, "as_absent"), (True, "excluded")):
            st, dg = run_clone(z, D, ci, int(sizes[ci]), exc)
            row[tag] = st
            if exc: degs[ci] = dg
        res.append(row)
    out = RES / f"compatibility_{tbl}.json"
    out.write_text(json.dumps(dict(table=tbl, min_carriers=MIN_CARRIERS, clones=res), indent=1))

    def agg(k):
        P = sum(r[k]["pairs"] for r in res); I = sum(r[k]["incompatible"] for r in res)
        return P, I, 100*(P-I)/P if P else float("nan")
    print(f"{tbl}: {len(res):,} clones, {sum(r['characters'] for r in res):,} characters")
    for k, lab in (("as_absent","missing-as-absent"), ("excluded","missing-excluded")):
        P, I, f = agg(k)
        print(f"  {lab:>18}: {P:>14,} pairs | {I:>12,} incompatible | "
              f"compatibility {f:6.2f}%")
    P0,I0,f0 = agg("as_absent"); P1,I1,f1 = agg("excluded")
    print(f"  {'spread':>18}: {f1-f0:+.2f} percentage points")
    # conflict-graph concentration, missing-excluded
    alld = []
    for ci, dg in degs.items(): alld.extend(dg.values())
    if alld:
        a = np.sort(np.array(alld))[::-1]; c = np.cumsum(a)/a.sum()
        nchars = sum(r["characters"] for r in res)
        print(f"  conflict graph (excluded): {len(a):,} of {nchars:,} characters have degree>0 "
              f"({100*len(a)/nchars:.1f}%)")
        for pct in (1,5,10):
            k = max(int(len(a)*pct/100)-1, 0)
            print(f"    top {pct:>2}% of conflicting characters carry {100*c[k]:5.1f}% of all conflict edges")
    print(f"  -> {out.name}   [{time.time()-t0:.0f}s]")
