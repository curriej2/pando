"""Is `None` tape-dropout, or is it 'site not yet edited'?

Sequential filling means an observed tape should read: symbols at sites 1..j,
then nothing at j+1..6. If `None` marks unedited sites, missingness is a SUFFIX
within every tape. If `None` marks a lost amplicon, whole tapes go None at once.
An `None` with a symbol AFTER it is neither -- that would be real per-site loss.
"""
import csv, re
from collections import Counter
from pathlib import Path

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
SYMBOL = re.compile(r"^[ACGT]{4}GGA$")

f = DATA / "Mouse3_EditTable_filtered.csv"
with open(f, newline="") as fh:
    rdr = csv.reader(fh)
    hdr = next(rdr)[1:]
    tapes = sorted({c.rsplit(".", 1)[0] for c in hdr})
    idx = {t: [i for i, c in enumerate(hdr) if c.rsplit(".", 1)[0] == t] for t in tapes}

    pat = Counter()
    fill = Counter()
    n_tape_instances = 0
    for row in rdr:
        vals = row[1:]
        for t in tapes:
            v = [vals[i] for i in idx[t]]
            n_tape_instances += 1
            present = [x != "None" for x in v]
            k = sum(present)
            if k == 0:
                pat["all None (tape absent)"] += 1
                continue
            # is the non-None set exactly a prefix 1..k ?
            if present[:k] == [True]*k and not any(present[k:]):
                pat["prefix filled, None suffix"] += 1
                fill[k] += 1
            else:
                pat["INTERNAL None (gap)"] += 1
                fill[k] += 1

print(f"tape-instances examined (cells x 166 tapes): {n_tape_instances:,}\n")
for k, v in pat.most_common():
    print(f"  {k:<32} {v:>12,}  {100*v/n_tape_instances:6.2f}%")

print("\nnumber of filled sites per OBSERVED tape (edits laid down so far):")
tot = sum(fill.values())
for k in sorted(fill):
    print(f"  {k} site(s): {fill[k]:>10,}  {100*fill[k]/tot:6.2f}%")
