#!/usr/bin/env python3
r"""
================================================================================
 Ground-truth test: does the compatible skeleton respect true colony boundaries?
================================================================================

THE WORRY. Under missing-as-absent, dropout SHRINKS clades. If true clade
{a,b,c,d} loses b and d, we observe {a,c} -- which is not a clade of the true tree
at all. Compatibility filtering removes characters that CONFLICT, but a set of
false clades can be mutually consistent. Nothing in the procedure guarantees the
skeleton is CORRECT, only that it is SELF-CONSISTENT. Script 17 showed the
obvious defence fails: 84% of distinct clades are asserted by a single tape, so
cross-tape support cannot filter artifacts.

THE TEST. The Subclone arm is ground truth. Ten wells were single-cell sorted and
grown as monoclonal colonies, so cells in different colonies share no ancestry
after day 10. Park validated their own trees this way (Figure 1e,f). We pool cells
across colonies, build characters BLIND to the barcodes, build the skeleton, and
ask whether its clades respect colony boundaries.

⚠ ClonalBC is a good but INEXACT proxy for colony. At MOI 0.3 about 14% of
founding cells carry >=2 barcode integrations, so one colony can appear as two
ClonalBC groups (the pipeline assigns each cell its highest-UMI barcode). We see
11 substantial groups for 8 usable colonies, and groups of 21 and 14 cells matching
the paper's two failed wells ("14 and 22 cells"). So three outcomes are distinct:

    WITHIN    clade sits inside one group          -> fine, within-colony structure
    UNION     clade is exactly a union of whole groups
                                                   -> a split colony correctly
                                                      rejoined (expect ~3), OR two
                                                      colonies wrongly merged
    STRADDLE  clade takes part of one group and part/all of another
                                                   -> a genuine FALSE CLADE, which
                                                      is the dropout failure mode

The straddle rate, measured against ground truth, is the number we want.
================================================================================
"""
import csv, json, re, sys, importlib.util
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT_DIR / "results"
MIN_COUNT, MIN_CARRIERS, TAPE_FILTER = 1000, 3, 100
PER_GROUP = int(sys.argv[1]) if len(sys.argv) > 1 else 200
RESTARTS = int(sys.argv[2]) if len(sys.argv) > 2 else 500
SEED = 0
SCREEN = (sys.argv[3].lower() not in ('0','off','no')) if len(sys.argv) > 3 else True

spec = importlib.util.spec_from_file_location("s17", Path(__file__).parent / "17_skeleton_linear.py")
s17 = importlib.util.module_from_spec(spec)
_argv = sys.argv; sys.argv = ["x"]; spec.loader.exec_module(s17); sys.argv = _argv

d = json.loads((RES / "xi_vectors.json").read_text())["Subclone"]
keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}

clone = {}
with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if r["ClonalBC"] and r["ClonalBC"] != "None": clone[r["CellID"]] = r["ClonalBC"]

# ---- read, filter, group ----------------------------------------------------
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
        vals = row[1:]
        if sum(1 for t in tapes if any(vals[i] != "None" for i in blocks[t])) < TAPE_FILTER: continue
        by_group[g].append([tuple(vals[i] for i in blocks[t]) for t in tapes])

rng = np.random.default_rng(SEED)
groups = [g for g, v in sorted(by_group.items(), key=lambda kv: -len(kv[1])) if len(v) >= 50]
cells, labels = [], []
for g in groups:
    v = by_group[g]
    pick = rng.choice(len(v), min(PER_GROUP, len(v)), replace=False)
    for i in pick: cells.append(v[i]); labels.append(g)
labels = np.array(labels); n = len(cells)
print(f"pooled {n:,} cells from {len(groups)} ClonalBC groups "
      f"(<= {PER_GROUP} each), {len(tapes)} tapes\n")
for g in groups: print(f"    {g:>12}  {len(by_group[g]):>6,} available  "
                       f"{int((labels==g).sum()):>4} sampled")


# ---- Park's cross-clone collision screen ------------------------------------
# Methods: "every cell of a clone was screened against the consensus of all other
# clones of the same mouse. A cell was flagged ... when its sequential match to its
# own clone's consensus fell below 0.50 while another clone led it by at least 0.20,
# on at least five co-observed tapes, and flagged cells were pruned."
#
# ⚠ The paper does not define "sequential match" numerically, so we pick one and
# state it: consensus prefix per (group, tape) is the longest prefix carried by a
# MAJORITY of that group's cells observing the tape (majority, so a minority of
# contaminants cannot set it). Then
#
#     match(cell, group) = mean over co-observed tapes of
#                          lcp(cell prefix, consensus prefix) / len(consensus prefix)
#
# A cell in the RIGHT group shares its founder's edits, so match is high. A cell
# from a DIFFERENT founding clone edited independently from day 0 and matches only
# by chance, at roughly q = 0.017 per site.
def consensus(recs, ntapes):
    """Per tape: longest prefix held by >=50% of the cells that observe that tape."""
    out = []
    for ti in range(ntapes):
        prefs = []
        for r in recs:
            pr = ()
            for v in r[ti]:
                if v not in keep: break
                pr += (v,)
            if pr: prefs.append(pr)
        if not prefs: out.append(()); continue
        need = len(prefs) / 2.0
        cons = ()
        while True:
            nxt = Counter(p[len(cons)] for p in prefs if len(p) > len(cons)
                          and p[:len(cons)] == cons)
            if not nxt: break
            sym, cnt = nxt.most_common(1)[0]
            if cnt < need: break
            cons = cons + (sym,)
        out.append(cons)
    return out


def match(rec, cons, ntapes):
    tot, k = 0.0, 0
    for ti in range(ntapes):
        c = cons[ti]
        if not c: continue
        pr = ()
        for v in rec[ti]:
            if v not in keep: break
            pr += (v,)
        if not pr: continue
        l = 0
        while l < min(len(pr), len(c)) and pr[l] == c[l]: l += 1
        tot += l / len(c); k += 1
    return (tot / k if k else 0.0), k


if SCREEN:
    cons = {g: consensus([cells[i] for i in np.flatnonzero(labels == g)], len(tapes))
            for g in groups}
    flag = np.zeros(n, bool)
    for i in range(n):
        own, k = match(cells[i], cons[labels[i]], len(tapes))
        if k < 5: continue
        best = max((match(cells[i], cons[g], len(tapes))[0]
                    for g in groups if g != labels[i]), default=0.0)
        if own < 0.50 and best >= own + 0.20: flag[i] = True
    print(f"\\n  collision screen: flagged {int(flag.sum()):,}/{n:,} cells "
          f"({100*flag.mean():.2f}%) — pruned")
    for g in groups:
        m = (labels == g)
        if m.sum(): print(f"      {g:>12} {int((m & flag).sum()):>4}/{int(m.sum()):>4} flagged")
    cells = [c for i, c in enumerate(cells) if not flag[i]]
    labels = labels[~flag]; n = len(cells)
    print(f"  {n:,} cells remain\\n")

# ---- build characters over the POOLED set, blind to labels -------------------
trie = defaultdict(list)
for ci, rec in enumerate(cells):
    for ti, t in enumerate(tapes):
        pref = ()
        for v in rec[ti]:
            if v not in keep: break
            pref += (v,)
            trie[(ti, pref)].append(ci)
clades = [(len(m), 1, np.array(sorted(m), np.int32)) for m in trie.values()
          if MIN_CARRIERS <= len(m) < n]
# deduplicate, tracking distinct tapes
dedup = {}
for (ti, pref), m in trie.items():
    if not (MIN_CARRIERS <= len(m) < n): continue
    k = np.array(sorted(m), np.int32).tobytes()
    e = dedup.setdefault(k, [set(), np.array(sorted(m), np.int32), 99])
    e[0].add(ti); e[2] = min(e[2], len(pref))
clades = [(c.size, len(tp), c) for tp, c, _ in dedup.values()]
depth  = [dep for _, _, dep in dedup.values()]
print(f"\n  characters {sum(1 for (ti,p),m in trie.items() if MIN_CARRIERS<=len(m)<n):,}"
      f" -> distinct clades {len(clades):,}")

acc, par, owner = s17.best_skeleton(clades, n, RESTARTS)
print(f"  skeleton: C = {len(acc):,} of {n-1:,} internal nodes = {len(acc)/(n-1):.3f}"
      f"  ({RESTARTS} restarts)")

# ---- score against ground truth ---------------------------------------------
gsize = {g: int((labels == g).sum()) for g in groups}
tally, examples = Counter(), defaultdict(list)
by_depth = defaultdict(Counter)
purity = []
detail = []
for k in acc:
    S = clades[k][2]
    lab = Counter(labels[S].tolist())
    full = [g for g, c in lab.items() if c == gsize[g]]
    part = [g for g, c in lab.items() if c < gsize[g]]
    if len(lab) == 1:
        cls = "WITHIN"
    elif not part:
        cls = f"UNION({len(full)})"
    else:
        cls = "STRADDLE"
    tally[cls] += 1
    by_depth[min(depth[k], 4)][cls] += 1
    top = max(lab.values())
    purity.append((S.size, top / S.size, len(lab), S.size - top))
    detail.append((int(S.size), float(top / S.size), cls, int(depth[k]),
                   int(clades[k][1]), sorted(((g, c) for g, c in lab.items()
                                              if c >= 0.5 * gsize[g]), key=lambda x: -x[1])))
    if cls == "STRADDLE" and len(examples["s"]) < 5:
        examples["s"].append((int(S.size), dict(lab)))

print(f"\n  {'class':>12} {'clades':>8} {'%':>7}")
for c, v in sorted(tally.items(), key=lambda kv: -kv[1]):
    print(f"  {c:>12} {v:>8,} {100*v/len(acc):>6.1f}%")
straddle = tally["STRADDLE"]
print(f"\n  ⇒ FALSE-CLADE RATE (straddling true colony boundaries): "
      f"{straddle:,}/{len(acc):,} = {100*straddle/len(acc):.2f}%")
print(f"\n  straddle rate by the SHALLOWEST prefix depth asserting the clade:")
print(f"  {'depth':>7} {'clades':>8} {'within':>8} {'straddle':>9} {'rate':>8}")
for dep in sorted(by_depth):
    c = by_depth[dep]; tot = sum(c.values())
    lab = f"{dep}+" if dep == 4 else str(dep)
    print(f"  {lab:>7} {tot:>8,} {c['WITHIN']:>8,} {c['STRADDLE']:>9,} "
          f"{100*c['STRADDLE']/tot:>7.1f}%")
P = np.array([[a, b, c, d] for a, b, c, d in purity])
print(f"\n  PURITY = largest single-colony share of a clade's cells")
print(f"  {'purity':>12} {'clades':>8} {'%':>7} {'median size':>12} {'median foreign cells':>21}")
for lo, hi in [(0.0,0.5),(0.5,0.8),(0.8,0.95),(0.95,0.99),(0.99,1.01)]:
    m = (P[:,1] >= lo) & (P[:,1] < hi)
    if not m.any(): continue
    lab = f">={lo:.2f}" if hi > 1 else f"{lo:.2f}-{hi:.2f}"
    print(f"  {lab:>12} {int(m.sum()):>8,} {100*m.mean():>6.1f}% "
          f"{np.median(P[m,0]):>12.0f} {np.median(P[m,3]):>21.0f}")
print(f"\n  median purity {np.median(P[:,1]):.3f};  "
      f"{100*(P[:,1]>=0.95).mean():.1f}% of clades are >=95% one colony")
print(f"  clades needing <=5 cells removed to become pure: "
      f"{100*(P[:,3]<=5).mean():.1f}%;  <=20 cells: {100*(P[:,3]<=20).mean():.1f}%")
# ---- are the LARGE merges shallow (homoplasy) or deep (something else)? ------
print(f"\n  DEPTH & SUPPORT by clade size and class")
print(f"  {'size':>12} {'class':>9} {'n':>6} {'med depth':>10} {'med support':>12} {'med purity':>11}")
for lo, hi in [(3,9),(10,49),(50,199),(200,100000)]:
    for cls in ("WITHIN", "STRADDLE"):
        sub = [d for d in detail if lo <= d[0] <= hi and d[2] == cls]
        if not sub: continue
        print(f"  {f'{lo}-{hi}':>12} {cls:>9} {len(sub):>6,} "
              f"{np.median([d[3] for d in sub]):>10.1f} "
              f"{np.median([d[4] for d in sub]):>12.1f} "
              f"{np.median([d[1] for d in sub]):>11.3f}")

big = sorted([d for d in detail if d[0] >= 50 and d[2] == "STRADDLE"], key=lambda d: -d[0])
print(f"\n  the {len(big)} straddling clades with >=50 cells "
      f"(colonies contributing >=50% of their sampled cells):")
print(f"  {'size':>7} {'depth':>6} {'support':>8} {'purity':>7}  whole colonies merged")
for d in big[:14]:
    print(f"  {d[0]:>7,} {d[3]:>6} {d[4]:>8} {d[1]:>7.3f}  "
          + ", ".join(f"{g}({c})" for g, c in d[5]))
bw = sorted([d for d in detail if d[0] >= 50 and d[2] == "WITHIN"], key=lambda d: -d[0])
print(f"\n  for contrast, the largest {min(len(bw),6)} WITHIN-colony clades >=50 cells:")
for d in bw[:6]:
    print(f"  {d[0]:>7,} {d[3]:>6} {d[4]:>8} {d[1]:>7.3f}  "
          + ", ".join(f"{g}({c})" for g, c in d[5]))
for sz, lab in examples["s"]:
    print(f"      e.g. clade of {sz} cells spanning {lab}")
(RES / "subclone_groundtruth.json").write_text(json.dumps(
    dict(cells=n, groups={g: gsize[g] for g in groups}, C=len(acc),
         frac=len(acc)/(n-1), tally=dict(tally),
         false_clade_rate=straddle/len(acc)), indent=1))
