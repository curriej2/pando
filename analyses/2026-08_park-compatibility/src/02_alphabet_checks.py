#!/usr/bin/env python3
"""Two checks on the Step-0 alphabet result.

(a) RAREFACTION. M_obs ranges 167-244 across tables. Is that biology, or just
    that the tables differ 17-fold in edit count? Subsample every table down to
    the smallest table's edit count and re-count distinct symbols.

(b) PER-SITE xi. The Step-0 estimator pools over cells, tapes AND sites. That is
    unbiased for xi only under row A3 (insert probabilities shared across tapes,
    sites and time). Sites are ordered in time, so per-site q is a direct
    empirical test of A3 on real data.
"""
import csv, re, random, json
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
OUT = Path(__file__).resolve().parents[1] / "results"
SYMBOL = re.compile(r"^[ACGT]{4}GGA$")
TABLES = ["Initial", "Mouse1", "Mouse2", "Mouse3", "Subclone"]
RAREFY_TO = 1_322_277        # Mouse3's conforming-edit count, the smallest
SEED = 0


def collect(path):
    """Return (global Counter, {site: Counter})."""
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        site_of = [c.rsplit(".", 1)[1] for c in next(rdr)[1:]]
        glob = Counter()
        per_site = defaultdict(Counter)
        for row in rdr:
            for v, s in zip(row[1:], site_of):
                if SYMBOL.match(v):
                    glob[v] += 1
                    per_site[s][v] += 1
    return glob, per_site


def q_of(counter):
    n = sum(counter.values())
    return sum((c / n) ** 2 for c in counter.values()), n


def rarefy_distinct(counter, target, seed=SEED):
    """Distinct symbols expected in a sample of `target` edits, without replacement."""
    total = sum(counter.values())
    if total <= target:
        return len(counter), False
    rng = random.Random(seed)
    syms, weights = zip(*counter.items())
    # sample `target` edits without replacement from the multiset
    pool = []
    for s, w in zip(syms, weights):
        pool.append((s, w))
    # hypergeometric-style: draw sequentially by shuffling an index space is too
    # big; instead sample indices and map back via cumulative counts.
    picks = rng.sample(range(total), target)
    picks.sort()
    seen, cum, j = set(), 0, 0
    for s, w in pool:
        hi = cum + w
        while j < len(picks) and picks[j] < hi:
            seen.add(s); j += 1
            # once this symbol is seen, skip its remaining picks
            while j < len(picks) and picks[j] < hi:
                j += 1
        cum = hi
    return len(seen), True


def main():
    rows, site_rows = [], []
    for t in TABLES:
        glob, per_site = collect(DATA / f"{t}_EditTable_filtered.csv")
        q, n = q_of(glob)
        m_rare, did = rarefy_distinct(glob, RAREFY_TO)
        rows.append(dict(table=t, M_obs=len(glob), edits=n, q=q,
                         M_rarefied=m_rare, rarefied=did))
        for s in sorted(per_site):
            qs, ns = q_of(per_site[s])
            site_rows.append(dict(table=t, site=s, M_obs=len(per_site[s]),
                                  edits=ns, q=qs))

    print("(a) RAREFACTION — all tables cut to %s edits\n" % f"{RAREFY_TO:,}")
    print(f"{'table':>10}{'edits':>13}{'M_obs':>8}{'M@rarefied':>13}")
    for r in rows:
        tag = "" if r["rarefied"] else "  (already smallest)"
        print(f"{r['table']:>10}{r['edits']:>13,}{r['M_obs']:>8}{r['M_rarefied']:>13}{tag}")

    print("\n(b) PER-SITE q  — tests row A3 (xi constant across sites/time)\n")
    print(f"{'table':>10}" + "".join(f"{s:>11}" for s in ["Site1","Site2","Site3","Site4","Site5","Site6"]))
    for t in TABLES:
        qs = {r["site"]: r["q"] for r in site_rows if r["table"] == t}
        print(f"{t:>10}" + "".join(f"{qs.get(s, float('nan')):>11.5f}" for s in
              ["Site1","Site2","Site3","Site4","Site5","Site6"]))
    print(f"\n{'table':>10}" + "".join(f"{s:>11}" for s in ["n1","n2","n3","n4","n5","n6"]))
    for t in TABLES:
        ns = {r["site"]: r["edits"] for r in site_rows if r["table"] == t}
        print(f"{t:>10}" + "".join(f"{ns.get(s,0):>11,}" for s in
              ["Site1","Site2","Site3","Site4","Site5","Site6"]))

    OUT.mkdir(exist_ok=True)
    (OUT / "alphabet_checks.json").write_text(
        json.dumps(dict(rarefaction=rows, per_site=site_rows), indent=2))


if __name__ == "__main__":
    main()
