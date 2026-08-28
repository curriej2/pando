#!/usr/bin/env python3
"""Step 0 of the D.4b protocol: first contact with the Park edit tables.

Streams each *_EditTable_filtered.csv and reports, per table:
  - matrix shape (cells x tapes x sites) and the tape/site vocabulary
  - the missingness rate (the literal string 'None')
  - the insertion-symbol alphabet, split into NNNNGGA-conforming and not
  - q = sum_i xi_i^2 over conforming symbols  -- the quantity that controls
    homoplasy, and the only channel by which composition affects topology.

Standard library only: no project env needed, and streaming keeps peak RSS
to a few MB regardless of the 268 MB inputs.
"""
import csv, re, sys, json
from collections import Counter
from pathlib import Path

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
OUT = Path(__file__).resolve().parents[1] / "results"
SYMBOL = re.compile(r"^[ACGT]{4}GGA$")   # Park's NNNNGGA insertion design
MISSING = "None"
TABLES = ["Initial", "Mouse1", "Mouse2", "Mouse3", "Subclone"]


def profile(path):
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        assert header[0] == "CellID", header[0]
        cols = header[1:]
        tapes, sites = set(), set()
        for c in cols:
            tape, site = c.rsplit(".", 1)
            tapes.add(tape); sites.add(site)

        good, bad = Counter(), Counter()
        n_missing = n_cells = 0
        for row in rdr:
            n_cells += 1
            for v in row[1:]:
                if v == MISSING:
                    n_missing += 1
                elif SYMBOL.match(v):
                    good[v] += 1
                else:
                    bad[v] += 1

    n_good, n_bad = sum(good.values()), sum(bad.values())
    n_obs = n_good + n_bad
    total = n_obs + n_missing
    q = sum((c / n_good) ** 2 for c in good.values()) if n_good else float("nan")
    return dict(
        table=path.stem.replace("_EditTable_filtered", ""),
        cells=n_cells, tapes=len(tapes), sites=len(sites), columns=len(cols),
        entries=total, missing=n_missing, missing_frac=n_missing / total,
        symbols_observed=len(good), symbols_possible=4 ** 4,
        edits_conforming=n_good, edits_nonconforming=n_bad,
        nonconforming_frac_of_observed=n_bad / n_obs if n_obs else 0.0,
        nonconforming_distinct=len(bad),
        q=q, effective_alphabet=1 / q, q_flat_floor=1 / len(good),
        top_symbol=good.most_common(1)[0][0],
        top_symbol_frac=good.most_common(1)[0][1] / n_good,
    )


def main():
    rows = [profile(DATA / f"{t}_EditTable_filtered.csv") for t in TABLES]
    OUT.mkdir(exist_ok=True)
    (OUT / "symbol_composition.json").write_text(json.dumps(rows, indent=2))

    hdr = ("table", "cells", "tapes", "sites", "missing%", "M_obs",
           "edits", "junk%", "q", "1/q")
    print("".join(f"{h:>12}" for h in hdr))
    for r in rows:
        print(f"{r['table']:>12}{r['cells']:>12,}{r['tapes']:>12}{r['sites']:>12}"
              f"{100*r['missing_frac']:>11.1f}%{r['symbols_observed']:>12}"
              f"{r['edits_conforming']:>12,}{100*r['nonconforming_frac_of_observed']:>11.2f}%"
              f"{r['q']:>12.6f}{r['effective_alphabet']:>12.1f}")
    print(f"\nnotes predicted q ~= 0.004 (flat over 256); pooled measured q "
          f"~= {sum(r['q'] for r in rows)/len(rows):.4f}")


if __name__ == "__main__":
    main()
