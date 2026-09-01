#!/usr/bin/env python3
r"""Compatibility test, sparse-matrix formulation.  Replaces script 13's engine.

Script 13's cross-tab is CORRECT (verified against brute force) but slow: it loops
over 166-choose-2 tape pairs x 36 level pairs ~ 493k numpy calls per clone, and on
the tiny arrays of a small clone the call overhead dominates -- 14.5 s per clone,
which is ~15 hours for the dataset and far worse for the big clones.

The whole computation is two sparse products.  Let

    M  = characters x cells, 1 where the cell is in S_i      (membership)
    N  = characters x cells, 1 where the cell is in D_i      (determined)

Then, for the missing-excluded convention, every witness count is an entry of a
matrix product:

    n11_ij = |S_i & S_j|          = (M M^T)_ij
    n10_ij = |S_i & D_j| - n11    = (M N^T)_ij   - (M M^T)_ij
    n01_ij = |S_j & D_i| - n11    = (M N^T)_ji   - (M M^T)_ij

    INCOMPATIBLE iff  n11 > 0  and  n10 > 0  and  n01 > 0

Since incompatibility requires n11 > 0, only the NON-ZERO entries of M M^T are ever
materialised -- i.e. exactly the co-occurring pairs.  Disjoint pairs (the vast
majority) are counted by subtraction from the analytic total and never touched.

For missing-as-absent, D is all cells, so N is all-ones and the products collapse:
n10 = |S_i| - n11, n01 = |S_j| - n11.  Only M M^T is needed.

Within-tape pairs are masked out afterwards (they are nested by construction).
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
MIN_CARRIERS = 3


def load(tbl):
    z = np.load(RES / f"characters_{tbl}.npz", allow_pickle=True)
    key = z["dmask_key"]; lens = z["dmask_lens"]; flat = z["dmask"]
    off = np.concatenate([[0], np.cumsum(lens)])
    D = {(int(a), int(b), int(c)): np.unpackbits(flat[off[i]:off[i+1]]).astype(bool)
         for i, (a, b, c) in enumerate(key)}
    return z, D


def clone_matrices(z, D, ci, n):
    """M = characters x cells membership; Ntilde = DISTINCT determination masks.

    D depends only on (tape, level), never on the prefix, so there are at most
    166*6 = 996 distinct masks per clone rather than one per character. Building
    N as characters x cells made M @ N.T a near-dense C x C product (241M entries
    for Mouse2's largest clone) and the subsequent B[r,c] fancy-index blew past
    128 G. Indexing a C x 996 matrix by mask id instead is ~4 orders smaller.
    """
    sel = np.flatnonzero((z["ch_clone"] == ci) & (z["ch_len"] >= MIN_CARRIERS))
    if sel.size == 0: return None
    tape = z["ch_tape"][sel].astype(np.int32)
    lev = z["ch_level"][sel].astype(np.int32)

    keys = {}
    mask_of = np.empty(sel.size, np.int32)
    rows_n, cols_n = [], []
    for r in range(sel.size):
        k = (int(tape[r]), int(lev[r]))
        if k not in keys:
            keys[k] = len(keys)
            d = D.get((ci, k[0], k[1]))
            dc = np.flatnonzero(d[:n]) if d is not None else np.array([], np.int32)
            rows_n.append(np.full(dc.size, keys[k], np.int32)); cols_n.append(dc)
        mask_of[r] = keys[k]

    rows_m, cols_m = [], []
    for r, i in enumerate(sel):
        car = z["carriers"][z["ch_start"][i]:z["ch_start"][i] + z["ch_len"][i]]
        rows_m.append(np.full(car.size, r, np.int32)); cols_m.append(car)

    C = sel.size
    M = sparse.csr_matrix((np.ones(sum(x.size for x in cols_m), np.int32),
                           (np.concatenate(rows_m), np.concatenate(cols_m))), shape=(C, n))
    Nt = sparse.csr_matrix((np.ones(sum(x.size for x in cols_n), np.int32),
                            (np.concatenate(rows_n) if rows_n else np.array([], np.int32),
                             np.concatenate(cols_n) if cols_n else np.array([], np.int32))),
                           shape=(len(keys), n))
    return M, Nt, mask_of, tape, np.asarray(M.sum(axis=1)).ravel()


def analyse(M, Nt, mask_of, tape, sizes, excluded, budget=4e7, edges=None):
    """Blocked so peak memory is independent of clone size.

    M @ M.T is near-dense for a large clone (15,533^2 ~ 2.4e8 non-zeros for Mouse2's
    biggest), and materialising it whole is what OOM-killed the large clones even
    after the compact-mask fix. Processing BLOCK rows at a time caps the live
    product at BLOCK x C instead of C x C; `budget` is the target entry count, so
    memory is tunable and no longer scales with the square of the clone.
    """
    C = M.shape[0]
    _, cnt = np.unique(tape, return_counts=True)
    total = (C * C - int((cnt.astype(np.int64) ** 2).sum())) // 2
    P = np.asarray((M @ Nt.T).todense()) if excluded else None
    block = max(256, min(C, int(budget // max(C, 1))))
    deg = np.zeros(C, np.int64)
    inc_tot = co_tot = 0
    for s in range(0, C, block):
        e = min(s + block, C)
        A = (M[s:e] @ M.T).tocoo()
        r = A.row.astype(np.int64) + s          # global row index
        c = A.col.astype(np.int64)
        n11 = A.data
        del A
        keep = (r < c) & (tape[r] != tape[c]) & (n11 > 0)
        r, c, n11 = r[keep], c[keep], n11[keep]
        del keep
        co_tot += r.size
        if excluded:
            n10 = P[r, mask_of[c]] - n11
            n01 = P[c, mask_of[r]] - n11
        else:
            n10 = sizes[r] - n11
            n01 = sizes[c] - n11
        inc = (n10 > 0) & (n01 > 0)
        inc_tot += int(inc.sum())
        if inc.any():
            deg += np.bincount(np.concatenate([r[inc], c[inc]]), minlength=C)
            if edges is not None:
                edges.append(np.stack([r[inc], c[inc]]).astype(np.int32))
        del r, c, n11, n10, n01, inc
    return dict(characters=C, pairs=int(total), cooccurring=int(co_tot),
                incompatible=int(inc_tot),
                compatible=int(total - inc_tot)), deg


if __name__ == "__main__":
    tbl = sys.argv[1] if len(sys.argv) > 1 else "Mouse3"
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 1200.0
    z, D = load(tbl); sizes_cl = z["clone_sizes"]
    nchar = np.bincount(z["ch_clone"][z["ch_len"] >= MIN_CARRIERS], minlength=len(sizes_cl))
    order = sorted([i for i in range(len(sizes_cl)) if nchar[i] > 0], key=lambda i: nchar[i])

    # cross-check against script 13's verified engine on a few clones
    sys.path.insert(0, str(Path(__file__).parent))
    import importlib.util
    spec = importlib.util.spec_from_file_location("s13", Path(__file__).parent / "13_compatibility.py")
    s13 = importlib.util.module_from_spec(spec)
    sys.argv = ["x"]; spec.loader.exec_module(s13)
    print(f"{tbl}: cross-checking sparse engine against script 13\n")
    for ci in order[:3]:
        mats = clone_matrices(z, D, ci, int(sizes_cl[ci]))
        for exc, tag in ((False, "as-absent"), (True, "excluded")):
            fast, _ = analyse(*mats, exc)
            ref, _ = s13.run_clone(z, D, ci, int(sizes_cl[ci]), exc)
            ok = fast["pairs"] == ref["pairs"] and fast["incompatible"] == ref["incompatible"]
            print(f"  clone {ci:>4} {tag:>10}: pairs {fast['pairs']:>8,}/{ref['pairs']:>8,}"
                  f"  incompat {fast['incompatible']:>7,}/{ref['incompatible']:>7,}  "
                  f"{'PASS' if ok else 'FAIL'}")
            assert ok
    print("  sparse engine verified\n")

    res, degs_exc, degs_abs, edge_store, t0 = [], {}, {}, {}, time.time()
    for ci in order:
        if time.time() - t0 > budget:
            print(f"  [budget reached: {len(res)}/{len(order)} clones]"); break
        mats = clone_matrices(z, D, ci, int(sizes_cl[ci]))
        if mats is None: continue
        row = dict(clone=int(ci), cells=int(sizes_cl[ci]), characters=int(nchar[ci]))
        for exc, tag in ((False, "as_absent"), (True, "excluded")):
            ed = [] if not exc else None          # as-absent edges: the constructible graph
            st, dg = analyse(*mats, exc, edges=ed)
            row[tag] = st
            (degs_exc if exc else degs_abs)[ci] = dg
            if ed and sum(e.shape[1] for e in ed) <= 60_000_000:
                edge_store[ci] = np.concatenate(ed, axis=1)
        res.append(row)
        if len(res) % 25 == 0 or nchar[ci] > 5000:
            print(f"    [{len(res):>5}/{len(order)}] clone {ci} "
                  f"({int(sizes_cl[ci]):,} cells, {int(nchar[ci]):,} chars) "
                  f"{time.time()-t0:.0f}s", flush=True)
    (RES / f"compatibility_{tbl}.json").write_text(
        json.dumps(dict(table=tbl, min_carriers=MIN_CARRIERS, clones=res), indent=1))
    # persist the conflict-graph degrees: the degree DISTRIBUTION is the D.4b step-5
    # deliverable, and recomputing it means rerunning the whole table.
    for tag, dd in (("excluded", degs_exc), ("as_absent", degs_abs)):
        if dd:
            np.savez_compressed(RES / f"conflict_degrees_{tag}_{tbl}.npz",
                clone=np.concatenate([np.full(d.size, ci, np.int32) for ci, d in dd.items()]),
                char=np.concatenate([np.arange(d.size, dtype=np.int32) for d in dd.values()]),
                degree=np.concatenate([d.astype(np.int32) for d in dd.values()]))
    if edge_store:
        np.savez_compressed(RES / f"conflict_edges_as_absent_{tbl}.npz",
            clone=np.concatenate([np.full(e.shape[1], ci, np.int32) for ci, e in edge_store.items()]),
            a=np.concatenate([e[0] for e in edge_store.values()]),
            b=np.concatenate([e[1] for e in edge_store.values()]))

    def agg(k):
        P = sum(r[k]["pairs"] for r in res); I = sum(r[k]["incompatible"] for r in res)
        return P, I, 100*(P-I)/P if P else float("nan")
    print(f"{tbl}: {len(res):,}/{len(order):,} clones, "
          f"{sum(r['characters'] for r in res):,} characters, {time.time()-t0:.0f}s")
    for k, lab in (("as_absent","missing-as-absent"), ("excluded","missing-excluded")):
        P, I, f = agg(k)
        print(f"  {lab:>18}: {P:>15,} pairs | {I:>13,} incompatible | compatibility {f:6.2f}%")
    print(f"  {'spread':>18}: {agg('excluded')[2]-agg('as_absent')[2]:+.2f} points")
    alld = np.concatenate([d[d > 0] for d in degs_exc.values() if (d > 0).any()]) if degs_exc else np.array([])
    if alld.size:
        a = np.sort(alld)[::-1]; cum = np.cumsum(a)/a.sum()
        nch = sum(r["characters"] for r in res)
        print(f"  conflict graph: {a.size:,}/{nch:,} characters conflict ({100*a.size/nch:.1f}%)")
        for p in (1,5,10,25):
            print(f"    top {p:>2}% of conflicting characters carry "
                  f"{100*cum[max(int(a.size*p/100)-1,0)]:5.1f}% of edges")
