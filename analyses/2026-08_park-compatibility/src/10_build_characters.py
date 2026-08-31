#!/usr/bin/env python3
r"""
================================================================================
 Character-set construction for the D.4b compatibility check   (steps 1-3)
================================================================================

A CHARACTER is a pair (tape z, prefix p), |p| >= 1.  Its clade S_{z,p} is the set
of cells whose tape z carries p as a prefix.  Two characters are COMPATIBLE iff
their clades are nested or disjoint; partial overlap means no single tree carries
both.  Only CROSS-TAPE pairs can ever conflict -- within a tape the prefixes form
a trie and are nested by construction (D.4 fact 1) -- which is what makes the
whole check cheap.

This script produces the characters.  The pairwise test is step 4, separately.

--------------------------------------------------------------------------------
 S AND D:  who carries it, versus whose membership we know
--------------------------------------------------------------------------------
Every cell sits in one of three buckets for a given character:

    in S            carries p on tape z
    in D, not in S  DETERMINED ABSENT -- we know it does not carry p
    outside D       UNKNOWN

S is always a subset of D: a cell undetermined for (z,p) is never in that
character's own clade.  The danger is CROSS-character -- a cell can be a perfectly
good member of S1 while being undetermined for character 2, on a different tape:

    cell | tape 1        in S1 | tape 2          in S2
    -----+---------------------+-----------------------------
      a  | reads X...    yes   | reads Y...      yes
      b  | reads X...    yes   | ENTIRELY None   UNDETERMINED
      c  | reads Z...    no    | reads Y...      yes
      d  | reads Z...    no    | reads W...      no

Naive set difference on S alone gives S1&S2={a}, S1\S2={b}, S2\S1={c}: all three
witnesses non-empty, so the pair is called INCOMPATIBLE.  But the only evidence
that b lies outside S2 is that its tape 2 could not be read.  The truth may well
be S1 subset of S2 -- nested, compatible.  With D, b cannot witness S1\S2, which
becomes (S1 & D2) \ S2 = empty, and the pair is correctly called nested.

At 24% tape-absence among analysed cells this mechanism would fire constantly.
D is what stops dropout from being read as evidence of absence.

    missing-excluded  : restrict every test to D1 & D2
    missing-as-absent : set D = all cells

Both conventions therefore fall out of one data structure (D.4b Procedure step 4
asks for both).

MEASURED COST OF TRUNCATING AT JUNK (Mouse1): 6.88% of (cell, tape) instances
contain a junk value; in 45% of those, readable symbols follow it, so truncation
discards 176,176 symbols = 3.03% of conforming edits.  Junk does not stop the tape
recording -- it stops us PARSING, and since a character is an ordered prefix from
site 1, one unreadable position makes every longer prefix unknown.
POSSIBLE REFINEMENT, NOT IMPLEMENTED: a post-junk symbol can never confirm
membership but can refute it -- a cell reading (s1, s2, JUNK, s4) is determined
ABSENT from any character whose 4th symbol is not s4.  That would enlarge D with
negative determinations only.  Revisit if the missing-excluded run is underpowered.

--------------------------------------------------------------------------------
 DETERMINATION LOGIC  -- richer than D.4b assumes, and in our favour
--------------------------------------------------------------------------------
`None` in the delivered tables conflates two different things, and script 03
showed they are separable: 44% of (cell, tape) instances are ENTIRELY None (the
tape was not recovered), while 9% of entries are trailing Nones inside a tape
that WAS recovered (those sites are simply not yet edited).

For a cell whose tape z was recovered and reads L symbols before terminating:

    |p| <= L                     -> determined (in if the prefix matches, else out)
    |p| >  L, stopped by `None`  -> determined ABSENT.  Irreversibility plus
                                    sequential filling mean the cell genuinely had
                                    not reached that site.  This is informative
                                    data, not missing data.
    |p| >  L, stopped by junk    -> undetermined; we cannot read past the junk.
    tape entirely None           -> undetermined.  (We cannot distinguish "not
                                    recovered" from "recovered but wholly
                                    unedited"; since site 1 is edited in 98.4% of
                                    recovered tapes, nearly all of these are
                                    dropout. Treating them as undetermined is the
                                    conservative call, and the missing-as-absent
                                    run covers the other reading.)

--------------------------------------------------------------------------------
 CELL FILTER (step 1)
--------------------------------------------------------------------------------
The paper keeps a cell for lineage reconstruction if it recovered >= 100 TargetBCs
(= tapes) for Initial/Subclone, or >= 20 for Mouse1-3.  The delivered
*_EditTable_filtered.csv do NOT have this applied (minimum observed is 0-1 tapes),
so we apply it here.  See notes S D.4c.

Alphabet rule, unchanged from scripts 06/08: a value is a symbol iff it matches
NNNNGGA *and* occurs >= MIN_COUNT times in its table; anything else terminates the
readable prefix as junk.

Output: results/characters_{table}.npz  (gitignored -- carries clone barcodes)
================================================================================
"""
import csv, re, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES  = ROOT / "results"
SYM  = re.compile(r"^[ACGT]{4}GGA$")
MIN_COUNT    = 1000     # symbol must be this common in its table
MIN_CARRIERS = 2        # store >=2; the >=3 filter is applied at analysis time
TAPE_FILTER  = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
N_SITES = 6

# termination codes
ABSENT, UNEDITED, JUNK, COMPLETE = 0, 1, 2, 3


def load_clones():
    clone = {}
    with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["ClonalBC"] and r["ClonalBC"] != "None":
                clone[r["CellID"]] = r["ClonalBC"]
    return clone


def parse_cell(vals, idxs, keep):
    """Read one tape for one cell -> (prefix tuple, read_len, termination code)."""
    pref = []
    for j, i in enumerate(idxs):
        v = vals[i]
        if v in keep:
            pref.append(v); continue
        # not a symbol: figure out why we stopped
        if v == "None":
            return tuple(pref), len(pref), (ABSENT if not pref else UNEDITED)
        return tuple(pref), len(pref), JUNK
    return tuple(pref), len(pref), COMPLETE


def run(tbl, clone, xi_all):
    d = xi_all[tbl]
    keep = {k for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
    thr = TAPE_FILTER[tbl]

    with open(DATA / f"{tbl}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh); hdr = next(r)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            tp, s = c.rsplit(".", 1); blocks[tp].append((int(s[4:]), i))
        tapes = sorted(blocks)
        tape_idx = {tp: k for k, tp in enumerate(tapes)}
        blocks = {tp: [i for _, i in sorted(v)] for tp, v in blocks.items()}

        by_clone = defaultdict(list)          # clone -> list of parsed cells
        n_read = n_nobc = n_lowtape = 0
        for row in r:
            n_read += 1
            cl = clone.get(row[0])
            if cl is None:
                n_nobc += 1; continue
            vals = row[1:]
            # step 1: the paper's tape-recovery filter
            n_tapes = sum(1 for tp in tapes if any(vals[i] != "None" for i in blocks[tp]))
            if n_tapes < thr:
                n_lowtape += 1; continue
            rec = [parse_cell(vals, blocks[tp], keep) for tp in tapes]
            by_clone[cl].append(rec)

    print(f"\n{tbl}: {n_read:,} cells | no ClonalBC {n_nobc:,} | "
          f"<{thr} tapes {n_lowtape:,} | kept {sum(len(v) for v in by_clone.values()):,} "
          f"in {len(by_clone):,} clones", flush=True)

    # step 3: characters per clone
    clone_names, clone_sizes = [], []
    ch_clone, ch_tape, ch_level, ch_start, ch_len = [], [], [], [], []
    carriers = []
    dmask = []                                 # (clone, tape, level) -> packed bits
    dmask_key = []
    for ci, (cl, cells) in enumerate(sorted(by_clone.items())):
        n = len(cells)
        clone_names.append(cl); clone_sizes.append(n)
        # D masks: determined(cell, tape, L) = tape recovered AND (L <= read_len OR stopped by None)
        for ti in range(len(tapes)):
            rl = np.array([cells[c][ti][1] for c in range(n)], dtype=np.int16)
            tm = np.array([cells[c][ti][2] for c in range(n)], dtype=np.int8)
            recovered = tm != ABSENT
            for L in range(1, N_SITES + 1):
                det = recovered & ((rl >= L) | (tm == UNEDITED) | (tm == COMPLETE))
                if det.any():
                    dmask.append(np.packbits(det))
                    dmask_key.append((ci, ti, L))
        # characters
        for ti in range(len(tapes)):
            trie = defaultdict(list)
            for c in range(n):
                pref = cells[c][ti][0]
                for L in range(1, len(pref) + 1):
                    trie[pref[:L]].append(c)
            for pref, mem in trie.items():
                if len(mem) < MIN_CARRIERS:
                    continue
                ch_clone.append(ci); ch_tape.append(ti); ch_level.append(len(pref))
                ch_start.append(len(carriers)); ch_len.append(len(mem))
                carriers.extend(mem)

    out = RES / f"characters_{tbl}.npz"
    np.savez_compressed(out,
        tapes=np.array(tapes), clone_names=np.array(clone_names),
        clone_sizes=np.array(clone_sizes, dtype=np.int32),
        ch_clone=np.array(ch_clone, dtype=np.int32),
        ch_tape=np.array(ch_tape, dtype=np.int16),
        ch_level=np.array(ch_level, dtype=np.int8),
        ch_start=np.array(ch_start, dtype=np.int64),
        ch_len=np.array(ch_len, dtype=np.int32),
        carriers=np.array(carriers, dtype=np.int32),
        dmask=np.concatenate(dmask) if dmask else np.array([], dtype=np.uint8),
        dmask_lens=np.array([len(x) for x in dmask], dtype=np.int32),
        dmask_key=np.array(dmask_key, dtype=np.int32))
    lv = np.array(ch_len)
    print(f"  characters (>= {MIN_CARRIERS} carriers): {len(ch_len):,}"
          f" | >=3 carriers: {(lv>=3).sum():,}"
          f" | clones {len(clone_names):,} (sizes: median {np.median(clone_sizes):.0f},"
          f" max {max(clone_sizes):,})")
    print(f"  -> {out.name} ({out.stat().st_size/1e6:.1f} MB)", flush=True)
    return dict(table=tbl, cells_kept=int(sum(clone_sizes)), clones=len(clone_names),
                characters_ge2=int(len(ch_len)), characters_ge3=int((lv>=3).sum()),
                dropped_no_barcode=n_nobc, dropped_low_tape=n_lowtape)


if __name__ == "__main__":
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    clone = load_clones()
    print(f"clone assignments: {len(clone):,} cells")
    summ = [run(t, clone, xi_all) for t in (sys.argv[1:] or list(TAPE_FILTER))]
    (RES / "character_summary.json").write_text(json.dumps(summ, indent=2))
    print("\nwrote results/character_summary.json")
