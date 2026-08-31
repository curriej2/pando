#!/usr/bin/env python3
r"""Is a junk value a heritable molecular event, or a readout artifact?

Scripts 06-10 treat any non-NNNNGGA value as junk that TRUNCATES the readable
prefix.  That is right for a sequencing artifact and wrong for a real insertion:
pegRNA scaffold read-through genuinely writes scaffold sequence into the tape, and
once written it is irreversible and inherited like any other edit.

The two hypotheses make opposite predictions, and the test needs no tree.

    A value V observed at site k of tape z arose ONCE in some ancestor
      => every carrier inherited that ancestor's sites 1..k-1
      => the PRECEDING PREFIXES of the carriers should be concordant.

    A readout artifact lands on phylogenetically unrelated cells
      => preceding prefixes should look like the background population.

Statistic: for each (tape, site k>=2, value), take the carriers whose sites 1..k-1
are all readable symbols, and compute

    concordance = (carriers sharing the single most common preceding prefix)
                  / (carriers)

Real NNNNGGA symbols are the positive control.  Concordance rises trivially as the
carrier count falls, so junk and real are compared WITHIN carrier-count bins.
"""
import csv, re, json, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = Path(__file__).resolve().parents[1] / "results"
SYM = re.compile(r"^[ACGT]{4}GGA$")
TBL = sys.argv[1] if len(sys.argv) > 1 else "Mouse1"

d = json.loads((RES / "xi_vectors.json").read_text())[TBL]
keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= 1000}

# (tape, site, value) -> Counter of preceding prefixes
pre = defaultdict(Counter)
with open(DATA / f"{TBL}_EditTable_filtered.csv", newline="") as fh:
    r = csv.reader(fh); hdr = next(r)[1:]
    blocks = defaultdict(list)
    for i, c in enumerate(hdr):
        tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
    blocks = {tp: [i for _, i in sorted(v)] for tp, v in blocks.items()}
    for row in r:
        vals = row[1:]
        for tp, idxs in blocks.items():
            v = [vals[i] for i in idxs]
            prefix = []
            for k, x in enumerate(v):
                if k >= 1 and x != "None":                  # site k+1, needs a prefix
                    if len(prefix) == k:                    # sites 1..k all readable
                        pre[(tp, k + 1, x)][tuple(prefix)] += 1
                if x in keep: prefix.append(x)
                else: break

rows = []
for (tp, k, val), c in pre.items():
    n = sum(c.values())
    if n < 3: continue
    rows.append((SYM.match(val) is not None and val in keep,
                 n, c.most_common(1)[0][1] / n, val))

real = [r for r in rows if r[0]]
junk = [r for r in rows if not r[0]]
print(f"{TBL}: {len(real):,} real (tape,site,symbol) triples | {len(junk):,} junk triples\n")
print(f"{'carriers':>12} {'real n':>8} {'real conc':>10} {'junk n':>8} {'junk conc':>10} {'verdict':>22}")
for lo, hi in [(3,4),(5,9),(10,24),(25,99),(100,499),(500,10**9)]:
    R = [r[2] for r in real if lo <= r[1] <= hi]
    J = [r[2] for r in junk if lo <= r[1] <= hi]
    if len(R) < 5 or len(J) < 5: continue
    mr, mj = float(np.median(R)), float(np.median(J))
    v = "junk ~ real (heritable)" if mj > 0.7*mr else ("intermediate" if mj > 0.3*mr else "junk ~ background")
    lab = f"{lo}-{hi if hi<10**9 else '+'}"
    print(f"{lab:>12} {len(R):>8,} {mr:>10.3f} {len(J):>8,} {mj:>10.3f} {v:>22}")

print("\ntop junk values by carrier count, with their concordance:")
print(f"{'value':>34} {'triples':>8} {'median carriers':>16} {'median conc':>12}")
byval = defaultdict(list)
for isreal, n, conc, val in junk: byval[val].append((n, conc))
for val, v in sorted(byval.items(), key=lambda kv: -sum(x[0] for x in kv[1]))[:10]:
    ns = [x[0] for x in v]; cs = [x[1] for x in v]
    print(f"{val[:34]:>34} {len(v):>8,} {np.median(ns):>16.0f} {np.median(cs):>12.3f}")
