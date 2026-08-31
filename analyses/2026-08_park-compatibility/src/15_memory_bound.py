#!/usr/bin/env python3
r"""Analytic peak-memory bound for the compatibility run (script 14).

Peak memory is set by nnz(M M^T) for the SINGLE LARGEST clone -- clones are
processed one at a time. A cell shared by k characters contributes C(k,2)
co-occurring pairs, so

    nnz(M M^T)  <=  sum_over_cells  k*(k-1)/2

at ~16 bytes per non-zero (int64 indices + data in scipy's COO/CSR).

This is an UPPER bound: a pair sharing several cells is counted once per shared
cell. It is the right quantity for sizing --mem because it cannot be exceeded.

Vectorised: the naive version (a Python list of per-character slices) took >280 s
on Subclone and was killed. Here the flat carrier ranges are expanded with the
standard repeat/arange trick, so it is a few seconds per table.
"""
import sys
from pathlib import Path
import numpy as np

RES = Path(__file__).resolve().parents[1] / "results"
TABLES = sys.argv[1:] or ["Mouse3", "Mouse2", "Mouse1", "Initial", "Subclone"]
MIN_CARRIERS = 3


print(f"{'table':>9} {'clone':>7} {'chars':>10} {'cells':>8} {'bound nnz':>16} "
      f"{'bound mem':>10} {'suggest --mem':>14}")
for t in TABLES:
    z = np.load(RES / f"characters_{t}.npz", allow_pickle=True)
    k = np.flatnonzero(z["ch_len"] >= MIN_CARRIERS)
    cl = z["ch_clone"][k]; st = z["ch_start"][k].astype(np.int64); ln = z["ch_len"][k].astype(np.int64)
    car = z["carriers"]; sizes = z["clone_sizes"]
    worst = None
    for ci in np.unique(cl):
        m = cl == ci
        n = int(sizes[ci])
        offs = np.repeat(st[m], ln[m]) + (np.arange(int(ln[m].sum()), dtype=np.int64)
                                          - np.repeat(np.cumsum(ln[m]) - ln[m], ln[m]))
        cnt = np.bincount(car[offs], minlength=n).astype(np.int64)
        b = int((cnt * (cnt - 1) // 2).sum())
        if worst is None or b > worst[1]:
            worst = (int(ci), b, int(m.sum()), n)
    ci, b, C, n = worst
    gb = b * 16 / 1e9
    sug = max(16, int(np.ceil(gb * 3 / 16) * 16))     # 3x headroom, round to 16G
    print(f"{t:>9} {ci:>7} {C:>10,} {n:>8,} {b:>16,} {gb:>9.1f}G {str(sug)+'G':>14}")
