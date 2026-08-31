#!/usr/bin/env python3
r"""Extract the full insert-probability vector xi for each table (and each site).

ALPHABET RULE, REVISED 2026-08-31 (see README "junk values ARE heritable"):
a value is a symbol iff it occurs >= MIN_COUNT times in its table. The old rule
also demanded the NNNNGGA design pattern; scripts 11 and 12 showed that
non-conforming values (scaffold read-through, length variants) are clone-restricted
FAR above a calibrated null -- more so than real symbols -- i.e. they are heritable
lineage events, not readout artifacts. Demanding conformance discarded real signal.
Frequency is the evidence of reality; matching the design pattern is not.

Emits two vectors:
  xi            over the kept alphabet (count >= MIN_COUNT), renormalised. Used by
                every downstream script.
  xi_untrimmed  over every value seen >= 2 times, renormalised. Only script 05
                needs it, to show that the rare tail fakes estimator resolution.

⚠ Scripts 01-03 predate this rule and still describe the NNNNGGA-only alphabet;
they are kept as the Step-0 record. Stdlib only.
"""
import csv, re, json
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
OUT = Path(__file__).resolve().parents[1] / "results"
SYM = re.compile(r"^[ACGT]{4}GGA$")      # only to REPORT how many were promoted
MIN_COUNT = 1000
TABLES = ["Initial", "Mouse1", "Mouse2", "Mouse3", "Subclone"]

out = {}
for t in TABLES:
    with open(DATA / f"{t}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh)
        site_of = [c.rsplit(".", 1)[1] for c in next(r)[1:]]
        glob = Counter(); per_site = defaultdict(Counter)
        for row in r:
            for v, s in zip(row[1:], site_of):
                if v != "None":
                    glob[v] += 1; per_site[s][v] += 1

    n_obs = sum(glob.values())
    keep = {k for k, c in glob.items() if c >= MIN_COUNT}
    kept = {k: glob[k] for k in keep}
    n = sum(kept.values())
    xi = {k: c / n for k, c in sorted(kept.items(), key=lambda kv: -kv[1])}
    q = sum(v * v for v in xi.values())
    promoted = sorted((k for k in keep if not SYM.match(k)), key=lambda k: -glob[k])

    tail = {k: c for k, c in glob.items() if c >= 2}
    nt = sum(tail.values())

    def site_block(c):
        sub = {k: v for k, v in c.items() if k in keep}
        tt = sum(sub.values())
        return {"n": tt, "M_obs": len(sub),
                "q": sum((v / tt) ** 2 for v in sub.values()),
                "xi": {k: v / tt for k, v in sorted(sub.items(), key=lambda kv: -kv[1])}}

    out[t] = dict(
        n_edits=n, n_observed=n_obs, M_obs=len(xi), q=q, xi=xi,
        promoted=promoted, n_promoted_edits=sum(glob[k] for k in promoted),
        xi_untrimmed={k: c / nt for k, c in sorted(tail.items(), key=lambda kv: -kv[1])},
        per_site={s: site_block(c) for s, c in per_site.items()},
    )
    print(f"{t:>10}  n={n:>10,}  M={len(xi):>4}  q={q:.6f}  "
          f"promoted {len(promoted):>2} non-NNNNGGA symbols "
          f"({100*sum(glob[k] for k in promoted)/n:.2f}% of kept edits)")

OUT.mkdir(exist_ok=True)
(OUT / "xi_vectors.json").write_text(json.dumps(out))
print(f"\nwrote {OUT/'xi_vectors.json'}")
