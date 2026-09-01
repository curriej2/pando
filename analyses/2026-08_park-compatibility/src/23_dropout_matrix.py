#!/usr/bin/env python3
r"""
================================================================================
 Cache the (cell x tape) dropout / depth matrices -- one pass over the CSVs
================================================================================

Everything in the Fig-3 dropout characterisation reads this cache, so the 627 MB
of edit tables are parsed exactly once.

For every (cell, tape) instance we record three things:

  recovered  bool   at least one of the tape's 6 sites is not the string "None"
  depth      int8   number of LEADING sites carrying an alphabet symbol
                    (= `read_len` in 10_build_characters.py: the prefix length)
  term       int8   why the read stopped, same codes as script 10:
                      ABSENT   0  every site is None -> tape not recovered
                      UNEDITED 1  symbols then None  -> trailing unedited sites
                      JUNK     2  symbols then a sub-threshold value
                      COMPLETE 3  all 6 sites carry symbols

`recovered` is exactly `term != ABSENT`, kept explicitly for clarity.

ALPHABET.  Frequency rule only (README "junk values ARE heritable characters"):
a value is a symbol iff it occurs >= MIN_COUNT times in its own table.  No
sequence-pattern requirement.  Read off xi_vectors.json, as scripts 10/22 do.

NOT APPLIED HERE: the paper's cell filter (>= 100 recovered tapes for
Initial/Subclone, >= 20 for the mice).  The filter is a selection on precisely
the quantity Fig 3b measures, so the cache stores every cell and the threshold is
applied at analysis time, both ways.

CELL METADATA.  `sample` comes from the CellID suffix (no join needed); `clone`
is the ClonalBC join, -1 where the cell has no barcode.  Cell IDs themselves are
not stored -- nothing downstream needs them.

Output: results/dropout_matrix_{table}.npz   (gitignored: *.npz)
================================================================================
"""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES  = ROOT / "results"
MIN_COUNT = 1000
N_SITES   = 6
ARMS      = ["Mouse3", "Mouse2", "Mouse1", "Initial", "Subclone"]   # ascending size
ABSENT, UNEDITED, JUNK, COMPLETE = 0, 1, 2, 3


def load_clones():
    clone = {}
    with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["ClonalBC"] and r["ClonalBC"] != "None":
                clone[r["CellID"]] = r["ClonalBC"]
    return clone


def run(tbl, clone, xi_all):
    d = xi_all[tbl]
    keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}

    with open(DATA / f"{tbl}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh)
        hdr = next(r)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            tp, s = c.rsplit(".", 1)
            blocks[tp].append((int(s[4:]), i))
        tapes = sorted(blocks)
        idxs = [[i for _, i in sorted(blocks[tp])] for tp in tapes]
        K = len(tapes)

        depth_rows, term_rows, samples, clones = [], [], [], []
        for row in r:
            vals = row[1:]
            dep = np.zeros(K, dtype=np.int8)
            trm = np.zeros(K, dtype=np.int8)
            for t, cols in enumerate(idxs):
                L, code = 0, COMPLETE
                for i in cols:
                    v = vals[i]
                    if v in keep:
                        L += 1
                        continue
                    code = (ABSENT if (v == "None" and L == 0)
                            else UNEDITED if v == "None" else JUNK)
                    break
                dep[t] = L
                trm[t] = code
            depth_rows.append(dep)
            term_rows.append(trm)
            cid = row[0]
            samples.append(cid.split("-1_", 1)[1] if "-1_" in cid else "?")
            clones.append(clone.get(cid, ""))

    depth = np.array(depth_rows, dtype=np.int8)
    term  = np.array(term_rows,  dtype=np.int8)
    recovered = term != ABSENT

    snames = sorted(set(samples))
    sidx = np.array([snames.index(s) for s in samples], dtype=np.int16)
    cnames = sorted({c for c in clones if c})
    cmap = {c: i for i, c in enumerate(cnames)}
    cidx = np.array([cmap.get(c, -1) for c in clones], dtype=np.int32)

    out = RES / f"dropout_matrix_{tbl}.npz"
    np.savez_compressed(
        out, tapes=np.array(tapes), depth=depth, term=term, recovered=recovered,
        sample=sidx, sample_names=np.array(snames),
        clone=cidx, clone_names=np.array(cnames))

    n, K = recovered.shape
    Rc = recovered.sum(1)
    print(f"{tbl}: {n:,} cells x {K} tapes | recovered {100*recovered.mean():.2f}% "
          f"| R_c mean {Rc.mean():.1f} sd {Rc.std(ddof=1):.1f} "
          f"min {Rc.min()} max {Rc.max()} | no ClonalBC {(cidx < 0).sum():,} "
          f"| clones {len(cnames):,} | samples {snames}", flush=True)
    print(f"  -> {out.name} ({out.stat().st_size/1e6:.1f} MB)", flush=True)


if __name__ == "__main__":
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    clone = load_clones()
    print(f"clone assignments: {len(clone):,} cells", flush=True)
    for t in (sys.argv[1:] or ARMS):
        run(t, clone, xi_all)
