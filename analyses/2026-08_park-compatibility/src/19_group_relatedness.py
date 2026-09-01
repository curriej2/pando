#!/usr/bin/env python3
r"""Are some ClonalBC groups actually ONE colony split by multi-integration barcoding?

The subclone colonies were single-cell sorted at day 10, ten days AFTER recording
began, so each founder already carried several edits per tape. Cells of one colony
therefore share a deep ancestral prefix on essentially every tape. Cells of two
DIFFERENT colonies descend from founders that diverged before day -11 and edited
independently from day 0, so they should share prefixes only by chance (~q=0.017).

So: build each group's majority consensus per tape and measure, for every pair of
groups, the mean shared-prefix depth across tapes. One colony split across two
barcodes => near-identical consensus. Two real colonies => near zero.

This decides whether the large "straddling" clades in script 18 are false clades or
a broken ground-truth label.
"""
import csv, json, sys
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = Path(__file__).resolve().parents[1] / "results"
MIN_COUNT, TAPE_FILTER = 1000, 100

d = json.loads((RES / "xi_vectors.json").read_text())["Subclone"]
keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
clone = {}
with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if r["ClonalBC"] and r["ClonalBC"] != "None": clone[r["CellID"]] = r["ClonalBC"]

by_group = defaultdict(list)
with open(DATA / "Subclone_EditTable_filtered.csv", newline="") as fh:
    rdr = csv.reader(fh); hdr = next(rdr)[1:]
    blocks = defaultdict(list)
    for i, c in enumerate(hdr):
        tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
    tapes = sorted(blocks); blocks = {t: [i for _, i in sorted(v)] for t, v in blocks.items()}
    for row in rdr:
        g = clone.get(row[0])
        if g is None: continue
        v = row[1:]
        if sum(1 for t in tapes if any(v[i] != "None" for i in blocks[t])) < TAPE_FILTER: continue
        rec = []
        for t in tapes:
            pr = ()
            for i in blocks[t]:
                if v[i] not in keep: break
                pr += (v[i],)
            rec.append(pr)
        by_group[g].append(rec)

groups = [g for g, v in sorted(by_group.items(), key=lambda kv: -len(kv[1])) if len(v) >= 50]
rng = np.random.default_rng(0)
cons = {}
for g in groups:
    v = by_group[g]
    idx = rng.choice(len(v), min(400, len(v)), replace=False)
    out = []
    for ti in range(len(tapes)):
        prefs = [v[i][ti] for i in idx if v[i][ti]]
        if not prefs: out.append(()); continue
        need, c = len(prefs) / 2.0, ()
        while True:
            nxt = Counter(p[len(c)] for p in prefs if len(p) > len(c) and p[:len(c)] == c)
            if not nxt: break
            s, n = nxt.most_common(1)[0]
            if n < need: break
            c += (s,)
        out.append(c)
    cons[g] = out

print(f"{len(groups)} ClonalBC groups; mean consensus depth per tape: " +
      ", ".join(f"{g[:6]}={np.mean([len(x) for x in cons[g]]):.2f}" for g in groups[:6]) + " ...\n")
print("mean SHARED PREFIX DEPTH between group consensuses (tapes where both non-empty)")
print("         " + "".join(f"{g[:6]:>8}" for g in groups))
M = np.zeros((len(groups), len(groups)))
for a, ga in enumerate(groups):
    row = ""
    for b, gb in enumerate(groups):
        sh, k = 0, 0
        for x, y in zip(cons[ga], cons[gb]):
            if not x or not y: continue
            l = 0
            while l < min(len(x), len(y)) and x[l] == y[l]: l += 1
            sh += l; k += 1
        M[a, b] = sh / k if k else 0
        row += f"{M[a,b]:>8.2f}"
    print(f"{ga[:8]:>8} {row}")
off = M[~np.eye(len(groups), dtype=bool)]
print(f"\n  off-diagonal: median {np.median(off):.3f}, max {off.max():.3f}")
hi = [(groups[a], groups[b], M[a, b]) for a in range(len(groups))
      for b in range(a + 1, len(groups)) if M[a, b] > 1.0]
print(f"  pairs sharing >1.0 symbols of prefix on average — i.e. plausibly ONE colony:")
for a, b, v in sorted(hi, key=lambda x: -x[2]): print(f"      {a} + {b}   {v:.2f}")
