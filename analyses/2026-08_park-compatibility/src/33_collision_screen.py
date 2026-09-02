#!/usr/bin/env python3
r"""
Park's cross-clone collision screen, vectorised, cached as per-cell flags.

WHY IT IS NEEDED HERE.  The lineage test asks whether cells of one clone lose the
same tapes.  A cell assigned to the wrong clone is unrelated to its clone-mates,
so it can only DILUTE within-clone agreement -- it biases the test toward the
null.  Reporting the test with and without the screen therefore brackets the
answer rather than picking a side.

THE RULE (paper methods, as implemented in 18_subclone_groundtruth.py): a cell is
flagged when its sequential match to its OWN clone's consensus falls below 0.50
while another clone of the same arm leads it by >= 0.20, on >= 5 co-observed tapes.

  consensus prefix per (clone, tape) = longest prefix carried by a MAJORITY of that
      clone's cells that observe the tape (majority, so a minority of contaminants
      cannot set it)
  match(cell, clone) = mean over co-observed tapes of
      lcp(cell prefix, consensus prefix) / len(consensus prefix)

⚠ DEVIATION, forced by the clone-size distribution.  The paper screens against
"all other clones of the same mouse"; with a median clone size of 2-4 cells most
"consensus" sequences here would be one cell's prefixes.  Competitor clones are
therefore restricted to those with >= MIN_COMP cells.  A cell's OWN clone is always
used, whatever its size.

Output: results/collision_flags_{arm}.npz  (flag, clone, n_cells)
"""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT / "results"
MIN_COUNT, N_SITES, MIN_COMP = 1000, 6, 5
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}


def load_clones():
    clone = {}
    with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["ClonalBC"] and r["ClonalBC"] != "None":
                clone[r["CellID"]] = r["ClonalBC"]
    return clone


def run(arm, clone_of, xi_all):
    d = xi_all[arm]
    keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
    sym_id = {s: i for i, s in enumerate(sorted(keep))}
    A = len(sym_id)
    thr = THR[arm]

    with open(DATA / f"{arm}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh); hdr = next(r)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
        tapes = sorted(blocks)
        idxs = [[i for _, i in sorted(blocks[tp])] for tp in tapes]
        K = len(tapes)

        rows, labels = [], []
        for row in r:
            cl = clone_of.get(row[0])
            if cl is None:
                continue
            v = row[1:]
            a = np.full((K, N_SITES), -1, dtype=np.int16)
            ntapes = 0
            for t, cols in enumerate(idxs):
                seen = False
                for j, i in enumerate(cols):
                    val = v[i]
                    if val != "None":
                        seen = True
                    sid = sym_id.get(val, -1)
                    if sid < 0:
                        break
                    a[t, j] = sid
                if seen:
                    ntapes += 1
            if ntapes < thr:
                continue
            rows.append(a); labels.append(cl)

    sym = np.array(rows, dtype=np.int16)          # (n, K, 6)
    n = sym.shape[0]
    names, inv = np.unique(np.array(labels), return_inverse=True)
    G = names.size
    sizes = np.bincount(inv, minlength=G)
    print(f"{arm}: {n:,} cells, {G:,} clones (median {np.median(sizes):.0f}, "
          f"max {sizes.max():,}); alphabet {A}", flush=True)

    # ---- consensus per (clone, tape): cost is linear in n, not n x G
    cons = np.full((G, K, N_SITES), -1, dtype=np.int16)
    order = np.argsort(inv, kind="stable")
    starts = np.concatenate([[0], np.cumsum(sizes)])
    tape_ix = np.arange(K)
    for g in range(G):
        idx = order[starts[g]:starts[g + 1]]
        S = sym[idx]                                # (m, K, 6)
        obs = S[:, :, 0] >= 0                       # observes the tape at all
        nobs = obs.sum(0)
        act = obs.copy()                            # still matching consensus
        for dpt in range(N_SITES):
            has = act & (S[:, :, dpt] >= 0)
            if not has.any():
                break
            code = (tape_ix[None, :] * A + S[:, :, dpt])[has]
            cnt = np.bincount(code, minlength=K * A).reshape(K, A)
            best = cnt.argmax(1); bc = cnt.max(1)
            ok = bc * 2 >= nobs                     # majority of observers
            ok &= bc > 0
            cons[g, ok, dpt] = best[ok].astype(np.int16)
            act &= (S[:, :, dpt] == cons[g, :, dpt][None, :]) & ok[None, :]
    conslen = (cons >= 0).sum(2)                    # (G, K)

    # ---- match every cell against every competitor consensus
    def score_all(g):
        c = cons[g]                                  # (K,6)
        L = conslen[g]
        eq = (sym == c[None]) & (c[None] >= 0)
        lcp = np.minimum.accumulate(eq, axis=2).sum(2)      # leading matches
        co = (sym[:, :, 0] >= 0) & (L > 0)[None, :]
        num = np.where(co, lcp / np.maximum(L, 1)[None, :], 0.0).sum(1)
        k = co.sum(1)
        return num / np.maximum(k, 1), k

    own = np.zeros(n); kown = np.zeros(n, int)
    for g in range(G):
        m = inv == g
        if not m.any():
            continue
        s, k = score_all(g)
        own[m] = s[m]; kown[m] = k[m]
    comp = np.where(sizes >= MIN_COMP)[0]
    best = np.zeros(n)
    for g in comp:
        s, _ = score_all(g)
        s = np.where(inv == g, -1.0, s)             # never compete with own clone
        best = np.maximum(best, s)
    flag = (kown >= 5) & (own < 0.50) & (best >= own + 0.20)

    print(f"  competitors: {comp.size} clones with >= {MIN_COMP} cells")
    print(f"  flagged {int(flag.sum()):,}/{n:,} cells ({100*flag.mean():.2f}%)",
          flush=True)
    np.savez_compressed(RES / f"collision_flags_{arm}.npz",
                        flag=flag, clone=inv, clone_names=names,
                        own=own, best=best, k_co=kown)
    print(f"  -> collision_flags_{arm}.npz", flush=True)


if __name__ == "__main__":
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    clone_of = load_clones()
    for a in (sys.argv[1:] or ["Mouse2", "Subclone"]):
        run(a, clone_of, xi_all)
