#!/usr/bin/env python3
"""Extract the full insert-probability vector xi for each table (and each site).

Step 0 saved only the summary q. The coupon-inversion estimator for m needs the
whole vector, since E[s|m] = sum_i (1 - (1-xi_i)^m) depends on the shape of xi,
not just on q. Stdlib only.
"""
import csv, re, json
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
OUT = Path(__file__).resolve().parents[1] / "results"
SYM = re.compile(r"^[ACGT]{4}GGA$")
TABLES = ["Initial", "Mouse1", "Mouse2", "Mouse3", "Subclone"]

out = {}
for t in TABLES:
    with open(DATA / f"{t}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh)
        site_of = [c.rsplit(".", 1)[1] for c in next(r)[1:]]
        glob = Counter(); per_site = defaultdict(Counter)
        for row in r:
            for v, s in zip(row[1:], site_of):
                if SYM.match(v):
                    glob[v] += 1; per_site[s][v] += 1
    n = sum(glob.values())
    xi = {k: c / n for k, c in glob.most_common()}
    q = sum(v * v for v in xi.values())
    out[t] = dict(
        n_edits=n, M_obs=len(xi), q=q,
        xi=xi,
        per_site={s: {"n": sum(c.values()), "M_obs": len(c),
                      "q": sum((v / sum(c.values()))**2 for v in c.values()),
                      "xi": {k: v / sum(c.values()) for k, v in c.most_common()}}
                  for s, c in per_site.items()},
    )
    print(f"{t:>10}  n={n:>10,}  M_obs={len(xi):>4}  q={q:.6f}")

OUT.mkdir(exist_ok=True)
(OUT / "xi_vectors.json").write_text(json.dumps(out))
print(f"\nwrote {OUT/'xi_vectors.json'}")
