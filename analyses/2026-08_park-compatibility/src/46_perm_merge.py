#!/usr/bin/env python3
r"""
================================================================================
 Pool the permutation parts into one null  (B = 1,000)
================================================================================

40_event_catalogue.py --permpart i/N and 43_soft_events.py --permpart i/N each
write ONE file per array task:

    results/permparts/permcounts_{script}_{arm}{tag}_p{i}.npz
        grid        (G,)     the fixed threshold grid -- identical in every part
        counts      (b, G)   #{candidates >= t} for each permutation in the part
        perm_index  (b,)     which of the B permutations those are
        meta        (5-6,)   B, N, i, seed, lam0 [, maxd]

This pools them into

    results/permnull_{script}_{arm}{tag}.npz     grid, counts (B, G), perm_index
    results/permnull_{script}_{arm}{tag}.json    committable summary

and CHECKS THE SET IS COMPLETE AND DISJOINT.  That check is the whole reason this
is a separate script: a silently missing part would shrink the null and inflate
every q-value's denominator, and a silently duplicated one would do the reverse.
Because permutation b is seeded from (SEED, b) rather than from a stream advanced
in order, duplicates are detectable by index -- and would have been detectable by
value too.

Nothing is averaged here.  The per-permutation count vectors are kept whole,
because the two quantities downstream need different things:

  * q-values need the MEAN null count at each threshold;
  * the global p-value needs the whole distribution --
    p = (1 + #{b : C_b >= C_obs}) / (B + 1).

Usage: 46_perm_merge.py <40|43> <arm> [--tag _d6] [--expect 1000]
================================================================================
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"

script = sys.argv[1]
arm = sys.argv[2]
TAG = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else ""
EXPECT = int(sys.argv[sys.argv.index("--expect") + 1]) if "--expect" in sys.argv else 1000

stem = f"permcounts_{script}_{arm}{TAG}_p"
parts = sorted((RES / "permparts").glob(stem + "*.npz"))
assert parts, f"no parts matching results/permparts/{stem}*.npz"

grid = None
counts, index, metas = [], [], []
for f in parts:
    z = np.load(f, allow_pickle=False)
    if grid is None:
        grid = z["grid"]
    assert np.array_equal(z["grid"], grid), f"{f.name}: grid differs from {parts[0].name}"
    assert z["counts"].shape[0] == z["perm_index"].size, f"{f.name}: ragged"
    assert z["counts"].shape[1] == grid.size, f"{f.name}: wrong grid width"
    counts.append(z["counts"])
    index.append(z["perm_index"])
    metas.append(z["meta"])

index = np.concatenate(index)
counts = np.concatenate(counts, 0)
seeds = {float(m[3]) for m in metas}
Bdecl = {float(m[0]) for m in metas}
assert len(seeds) == 1, f"parts disagree on the seed: {seeds}"
assert len(Bdecl) == 1, f"parts disagree on B: {Bdecl}"
B = int(Bdecl.pop())

uniq, cnt = np.unique(index, return_counts=True)
dup = uniq[cnt > 1]
assert dup.size == 0, f"permutation indices appear twice: {dup[:10].tolist()}"
missing = np.setdiff1d(np.arange(B), uniq)
assert missing.size == 0, (
    f"{missing.size} of B={B} permutations missing, e.g. {missing[:10].tolist()} "
    f"-- resubmit the array tasks that own them")
assert B == EXPECT, f"parts declare B={B}, expected {EXPECT}"

o = np.argsort(index)                      # canonical order, so the file is stable
counts, index = counts[o], index[o]

out = RES / f"permnull_{script}_{arm}{TAG}.npz"
np.savez_compressed(out, grid=grid, counts=counts.astype(np.int64),
                    perm_index=index.astype(np.int64),
                    meta=np.array([B, len(parts), float(list(seeds)[0])], dtype=float))

tot = counts[:, 0].astype(float)           # candidates above the scan floor
summ = {"script": script, "arm": arm, "tag": TAG, "B": B, "n_parts": len(parts),
        "seed": float(list(seeds)[0]), "lambda_scan_floor": float(grid[0]),
        "grid_points": int(grid.size),
        "null_total_candidates": {"mean": float(tot.mean()), "sd": float(tot.std(ddof=1)),
                                  "min": float(tot.min()), "max": float(tot.max()),
                                  "median": float(np.median(tot))},
        "null_mean_curve": {"lambda": grid.tolist(),
                            "mean": counts.mean(0).tolist(),
                            "max": counts.max(0).tolist()}}
(RES / f"permnull_{script}_{arm}{TAG}.json").write_text(json.dumps(summ, indent=1))
print(f"{script} {arm}{TAG}: pooled {len(parts)} parts -> B={B} permutations, "
      f"{grid.size} thresholds")
print(f"  null candidates above the {grid[0]:.0f}-nat scan floor: "
      f"mean {tot.mean():,.1f}  sd {tot.std(ddof=1):,.1f}  max {tot.max():,.0f}")
print(f"  wrote {out.relative_to(ROOT)} and {out.with_suffix('.json').name}")
