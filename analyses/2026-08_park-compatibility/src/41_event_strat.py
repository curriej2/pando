#!/usr/bin/env python3
r"""
Event catalogue, stratified by clone size -- because detection power is not
uniform and the pooled percentage is diluted by clones that could never be tested.

gamma_{C,z} is fitted from the m cells of clone C, so for a 5-cell clone it
absorbs essentially all of that clone's structure and NO subclade can be called.
Reporting "3.64% of missing entries" over all clones therefore mixes clones where
an event was detectable with clones where it was impossible.  Here both numerator
and denominator are restricted to the same clone-size band.

For each band: clones, cells, missing entries in those clones, events found, the
cell x tape slots inside events, and the missing entries inside them as a share of
the band's own missing entries.
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
ARMS = ["Mouse1", "Mouse2", "Mouse3", "Initial", "Subclone"]
THR = {"Initial": 100, "Subclone": 100, "Mouse1": 20, "Mouse2": 20, "Mouse3": 20}
BANDS = [(5, 10), (10, 20), (20, 50), (50, 200), (200, 10**9)]

out = {}
for arm in (sys.argv[1:] or ARMS):
    z = np.load(RES / f"dropout_matrix_{arm}.npz", allow_pickle=False)
    Y, clone = z["recovered"], z["clone"]
    s = (Y.sum(1) >= THR[arm]) & (clone >= 0)
    Y, clone = Y[s], clone[s]
    _, g = np.unique(clone, return_inverse=True)
    Gc = g.max() + 1
    sizes = np.bincount(g)
    miss_per_clone = np.array([(~Y[g == i]).sum() for i in range(Gc)])

    rows = [l.split("\t") for l in
            gzip.open(RES / f"events_{arm}.tsv.gz", "rt").read().splitlines()[1:]]
    ev_clone = np.array([int(r[0]) for r in rows], dtype=int)
    ev_slots = np.array([int(r[5]) for r in rows], dtype=int)
    ev_miss = np.array([int(r[6]) for r in rows], dtype=int)
    ev_exp = np.array([float(r[7]) for r in rows])
    ev_lam = np.array([float(r[9]) for r in rows])

    rec = []
    print(f"\n=== {arm}")
    print("  clone size |  clones    cells    missing | events  slots  missing_in_ev"
          "   share of band   excess   Lambda")
    for lo, hi in BANDS:
        cl = np.flatnonzero((sizes >= lo) & (sizes < hi))
        if cl.size == 0:
            continue
        m_band = int(miss_per_clone[cl].sum())
        sel = np.isin(ev_clone, cl)
        share = ev_miss[sel].sum() / max(m_band, 1)
        rec.append({"lo": lo, "hi": None if hi > 10**8 else hi,
                    "clones": int(cl.size), "cells": int(sizes[cl].sum()),
                    "missing_entries": m_band, "events": int(sel.sum()),
                    "slots_in_events": int(ev_slots[sel].sum()),
                    "missing_in_events": int(ev_miss[sel].sum()),
                    "share_of_band_missing": float(share),
                    "excess": float(ev_miss[sel].sum() - ev_exp[sel].sum()),
                    "lambda": float(ev_lam[sel].sum())})
        print(f"  {lo:>5}-{'inf' if hi>10**8 else hi:<5} | {cl.size:>6,} {sizes[cl].sum():>8,} "
              f"{m_band:>10,} | {sel.sum():>6,} {ev_slots[sel].sum():>6,} "
              f"{ev_miss[sel].sum():>13,} {100*share:>13.2f}% "
              f"{ev_miss[sel].sum()-ev_exp[sel].sum():>8.0f} {ev_lam[sel].sum():>9,.0f}")
    out[arm] = rec
(RES / "event_strat.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/event_strat.json")
