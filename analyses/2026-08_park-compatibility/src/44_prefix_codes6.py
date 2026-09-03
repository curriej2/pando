#!/usr/bin/env python3
r"""
Prefix codes to depth 6 -- the full array -- as prefix_codes6_{arm}.npz.

WHY.  The soft-variant catalogue (src/43) found pi_hat < 0.90 ("partial") for
57.2% of Subclone events against 8-32% in the other arms.  Two different things
are pooled under that label:

  (a) genuine graded silencing -- a heritable RATE shift, which would need a
      lineage-varying rate in the model;
  (b) a COMPLETE loss on a clade SMALLER than the one we tested -- an artefact of
      clade resolution, needing nothing but a finer clade.

Prefix clades stopped at depth 4 while Subclone's clones run 200-11,000 cells, so
(b) should dominate there.  Rebuilding to depth 6 -- the tapes have six sites, so
the resolution exists -- lets the test run at maximum resolution: if the partial
events are (b) they resolve into complete events at depth 5-6; if they persist
they are (a), and the graded component is real.

⚠ The cell filter and row order are kept BIT-IDENTICAL to 35_lineage_depth.py's
build_codes so the new cache aligns with dropout_matrix_{arm}.npz exactly (the
`seen` flag is set before the depth branch in both, so the tape count is
unchanged).  Written to a NEW filename so the depth-4 results stay reproducible.
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


def build(arm):
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
    n, K, _ = S.shape
    A = len(sym_id) + 1
    codes = np.full((n, K, MAX_D), -1, dtype=np.int32)
    prev = np.zeros((n, K), dtype=np.int64)
    ok = np.ones((n, K), dtype=bool)
    for d in range(MAX_D):
        ok &= S[:, :, d] >= 0
        prev = prev * A + (S[:, :, d].astype(np.int64) + 1)
        flat = np.where(ok, prev, -1).ravel()
        _, inv = np.unique(flat, return_inverse=True)
        codes[:, :, d] = np.where(ok, inv.reshape(n, K).astype(np.int32), -1)
        print(f"  {arm} depth {d+1}: {int(ok.sum()):,} determined (cell,tape) of {n*K:,} "
              f"({100*ok.mean():.1f}%), {len(np.unique(flat))-1:,} distinct prefixes", flush=True)
    ref = np.load(RES / f"prefix_codes_{arm}.npz", allow_pickle=False)["codes"]
    assert ref.shape[0] == n, f"row count differs from the depth-4 cache: {ref.shape[0]} vs {n}"
    assert (ref[:, :, :4] >= 0).sum() == (codes[:, :, :4] >= 0).sum(), \
        "determination pattern differs from the depth-4 cache"
    np.savez_compressed(RES / f"prefix_codes6_{arm}.npz", codes=codes)
    print(f"  -> prefix_codes6_{arm}.npz  (aligned with the depth-4 cache)", flush=True)


if __name__ == "__main__":
    for a in (sys.argv[1:] or ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]):
        build(a)
