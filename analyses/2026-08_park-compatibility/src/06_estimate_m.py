#!/usr/bin/env python3
"""First pass: estimate m per (clone, tape, prefix) from the distinct-symbol count s.

For every prefix p' on every tape within every clone, the carriers of p' show some
set of symbols at the next site. s = |that set|. The m independent write events that
produced them are unobservable, but E[s|m] = sum_i (1-(1-xi_i)^m) inverts against the
measured xi, and (m_hat - s) is then the estimated number of homoplastic recurrences.

ALPHABET RULE (one rule, used consistently for both xi and s): a value counts as a
symbol iff it matches NNNNGGA *and* occurs >=MIN_COUNT times in its table. Everything
else -- scaffold read-through, indel variants, the ultra-rare artifact tail -- is not
a symbol and TRUNCATES the readable prefix at that site. Costs 0.05% of edits.

Outputs results/m_estimates_{table}.tsv.gz and results/m_summary.json
"""
import csv, re, json, gzip, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT / "results"
SYM = re.compile(r"^[ACGT]{4}GGA$")
MIN_COUNT = 1000          # same cut as the identifiability analysis
MIN_CARRIERS = 3          # D.4b Procedure step 2: clades of >=3 cells
TABLES = sys.argv[1:] or ["Mouse1", "Mouse2", "Mouse3"]

# ---- clone assignment -------------------------------------------------------
clone = {}
with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
    for row in csv.DictReader(fh):
        bc = row["ClonalBC"]
        if bc and bc != "None":
            clone[row["CellID"]] = bc
print(f"clone assignments loaded: {len(clone):,} cells", flush=True)

xi_all = json.loads((RES / "xi_vectors.json").read_text())

# ---- E[s|m] inversion table, per site ---------------------------------------
M_GRID = np.unique(np.round(np.logspace(0, 5, 3000)).astype(int)).astype(float)

def build_inverter(xi_map):
    xi = np.array(list(xi_map.values())); xi = xi / xi.sum()
    exp_s = (1.0 - np.power(1.0 - xi[None, :], M_GRID[:, None])).sum(axis=1)
    q = float((xi ** 2).sum())
    def invert(s_vals):
        return np.interp(np.asarray(s_vals, dtype=float), exp_s, M_GRID)
    return invert, q, len(xi), exp_s[-1]

summary = {}
for tbl in TABLES:
    d = xi_all[tbl]
    keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
    inv = {}
    for site, sd in d["per_site"].items():
        sub = {k: v for k, v in sd["xi"].items() if k in keep}
        inv[site] = build_inverter(sub)
    print(f"\n{tbl}: alphabet = {len(keep)} symbols "
          f"({100*sum(d['xi'][k] for k in keep):.3f}% of edits)", flush=True)

    # ---- walk the prefix trie -----------------------------------------------
    nodes = defaultdict(Counter)          # (clone, tape, prefix) -> Counter(next symbol)
    carriers = Counter()                  # (clone, tape, prefix) -> cells carrying p'
    n_cells = n_skipped = 0
    with open(DATA / f"{tbl}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh)
        hdr = next(r)[1:]
        tape_of = [c.rsplit(".", 1)[0] for c in hdr]
        # column index blocks per tape, in site order
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            t, s = c.rsplit(".", 1)
            blocks[t].append((int(s[4:]), i))
        blocks = {t: [i for _, i in sorted(v)] for t, v in blocks.items()}

        for row in r:
            cid = row[0]
            cl = clone.get(cid)
            if cl is None:
                n_skipped += 1
                continue
            n_cells += 1
            vals = row[1:]
            for tape, idxs in blocks.items():
                pref = ()
                for j, i in enumerate(idxs):
                    v = vals[i]
                    key = (cl, tape, pref)
                    carriers[key] += 1
                    if v not in keep:          # None, junk, or ultra-rare -> stop
                        break
                    nodes[key][v] += 1
                    pref = pref + (v,)
    print(f"  cells used {n_cells:,} (skipped {n_skipped:,} without a clone barcode)"
          f" | trie nodes {len(nodes):,}", flush=True)

    # ---- per-node estimates -------------------------------------------------
    out = RES / f"m_estimates_{tbl}.tsv.gz"
    kept = 0
    rows_by_level = defaultdict(list)
    with gzip.open(out, "wt") as fo:
        fo.write("clone\ttape\tlevel\tclade_size\tn_next\ts\tm_hat\trecurrences\n")
        for key, cnt in nodes.items():
            cl, tape, pref = key
            n_next = sum(cnt.values())
            if n_next < MIN_CARRIERS:
                continue
            s = len(cnt)
            site = f"Site{len(pref)+1}"
            m_hat = float(inv[site][0]([s])[0])
            m_hat = max(m_hat, s)                       # m >= s always
            m_hat = min(m_hat, n_next)                  # m <= n_next always
            rec = m_hat - s
            kept += 1
            rows_by_level[len(pref)].append((n_next, s, m_hat, rec))
            fo.write(f"{cl}\t{tape}\t{len(pref)}\t{carriers[key]}\t{n_next}\t{s}"
                     f"\t{m_hat:.3f}\t{rec:.3f}\n")
    print(f"  nodes with >={MIN_CARRIERS} carriers: {kept:,}  -> {out.name}", flush=True)

    q0 = inv["Site1"][1]; birthday = float(np.sqrt(2 / q0))
    lv = {}
    for L, rows in sorted(rows_by_level.items()):
        a = np.array(rows)
        lv[L] = dict(n_nodes=len(a), median_n_next=float(np.median(a[:, 0])),
                     median_s=float(np.median(a[:, 1])),
                     median_m=float(np.median(a[:, 2])),
                     mean_recurrences=float(a[:, 3].mean()),
                     frac_with_recurrence=float((a[:, 3] >= 1).mean()),
                     frac_above_birthday=float((a[:, 2] >= birthday).mean()))
    summary[tbl] = dict(alphabet=len(keep), q_site1=q0, birthday=birthday,
                        cells=n_cells, nodes=kept, by_level=lv)
    print(f"  birthday threshold m* = {birthday:.1f}")
    print(f"  {'L':>2} {'nodes':>9} {'med n_next':>11} {'med s':>7} {'med m':>7}"
          f" {'mean recur':>11} {'% w/ recur':>11} {'% m>=m*':>9}")
    for L, v in lv.items():
        print(f"  {L:>2} {v['n_nodes']:>9,} {v['median_n_next']:>11.1f} {v['median_s']:>7.1f}"
              f" {v['median_m']:>7.1f} {v['mean_recurrences']:>11.3f}"
              f" {100*v['frac_with_recurrence']:>10.1f}% {100*v['frac_above_birthday']:>8.1f}%")

prev = json.loads((RES / "m_summary.json").read_text()) if (RES / "m_summary.json").exists() else {}
prev.update(summary)          # merge: a partial re-run must not drop other tables
(RES / "m_summary.json").write_text(json.dumps(prev, indent=2))
print("\nwrote results/m_summary.json")
