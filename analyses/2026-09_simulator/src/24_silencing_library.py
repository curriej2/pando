#!/usr/bin/env python3
r"""
24_silencing_library.py -- the rung-2 silencing library: layer D on top of the rung-1 library.

Approved by Justin 2026-10-01 (PROPOSAL (2026-10-01, rev. 1) §4b, ../CLAUDE.md, plus Hbar = 0.01
added at approval); generated only after 23_validate_silencing.py PASSED (seed 3102026).
One array task per tree cell (n, rho, theta), same task numbering as 19_dropout_library.py.
NOTHING is concluded from the library; the only summary is a manifest.

LAYER D.  Rung 1's layers A (20 library trees per cell) x B (5 panel/decoration replicates b) are
reused unchanged; silencing depends on the TREE only, so one set of loss draws E (Exp(1) per
branch x tape, 22's sampler) is drawn per (tree, b) -- independent across replicates, and shared by
every Lam_T, panel and mask of rung 1 at that (tree, b), and by every silencing setting (common
random numbers: nested across Hbar at a fixed shape).
GRID: Hbar in {0, .01, .02, .05, .10, .20, .33} x zeta_0 in
  constant;  early = 2 on [0, .4), 1/3 on [.4, 1)  (80% of the hazard in the first 40%);
  late = 1/3 on [0, .6), 2 on [.6, 1)              (80% in the last 40%)   -> 21 settings.

SEEDS.  E: (SEED, 10, task, tree, b).

STORED, one file per cell, results/silencing_library/sl_n{n}_rho{rho}_th{theta}.npz:
  Hp    (5, 20, 21, ceil(n k / 8)) uint8   tip H per setting, np.packbits of H.reshape(n*k)
  rec_b, rec_tree, rec_e, rec_z, rec_E     sparse loss record: every (branch, tape) whose draw E
                                           fires at Hbar_max under ANY shape -- with the tree,
                                           this rebuilds H at any setting with Hbar <= Hbar_max
  settings: list of (shape name, Hbar), in Hp's axis-2 order
Compose with a rung-1 file at the same (b, tree):  Y = ~X & ~H & (D >= 1)  (see compose()).

usage:  24_silencing_library.py <task 0..29> [--quick]
        24_silencing_library.py --collect          (manifests -> results/silencing_library_summary.json)
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, _HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


lib1 = _load("dropout_library", "19_dropout_library.py")
sl = _load("silencing", "22_silencing.py")
ed = lib1.ed

OUT = _HERE.parent / "results" / "silencing_library"
SUMMARY = _HERE.parent / "results" / "silencing_library_summary.json"
SEED = 1012026
K, TREES, REPS, CELLS = lib1.K, lib1.TREES, lib1.REPS, lib1.CELLS
HBARS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.33)
SHAPES = {"constant": ((0.0, 1.0), (1.0,)),
          "early": ((0.0, 0.4, 1.0), (2.0, 1.0 / 3.0)),
          "late": ((0.0, 0.6, 1.0), (1.0 / 3.0, 2.0))}
SETTINGS = [(s, h) for s in SHAPES for h in HBARS]


def params(shape, hbar):
    kn, lv = SHAPES[shape]
    return sl.SilencingParams(hbar, kn, lv)


def H_from_record(parent, t, n, E_sparse, sp):
    """Rebuild tip H at setting sp from the sparse record (E = +inf where nothing was stored)."""
    return sl.silence(parent, t, E_sparse, sp)["H"][:n]


def compose(X, D, H):
    """Rung-2 observed view: recovered AND not silenced AND edited at least once."""
    return sl.observe(X, D, H)


def run(task, quick):
    n, rho, th = CELLS[task]
    reps, ntree = (1, 2) if quick else (REPS, TREES)
    trees = ed.load_library_trees(n, rho, th)[:ntree]
    assert len(trees) == ntree
    OUT.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()
    nb = (n * K + 7) // 8
    Hp = np.zeros((reps, ntree, len(SETTINGS), nb), np.uint8)
    recs = {k_: [] for k_ in ("b", "tree", "e", "z", "E")}
    viol = {k_: [0, 0] for k_ in ("solve_exact", "dollo_monotone", "nested_across_Hbar",
                                  "size_weighted_V_equals_H", "record_rebuilds_H",
                                  "replay_identical", "Hbar0_empty")}
    for s, h in SETTINGS:
        sp = params(s, h)
        viol["solve_exact"][0] += int(abs(1 - np.exp(-sp.Gam_T) - h) > 1e-12)
        viol["solve_exact"][1] += 1
    rows = []
    hmax = max(HBARS)
    for b in range(reps):
        for ti, (parent, branch) in enumerate(trees):
            t = ed.node_times(parent, branch, n)
            sE = np.random.SeedSequence([SEED, 10, task, ti, b])
            E = sl.loss_draws(parent.size, K, np.random.default_rng(sE))
            sz = sl.clade_sizes(parent, n)
            # sparse record: fires at Hbar_max under any shape
            gmax = np.zeros(parent.size)
            for s in SHAPES:
                gmax = np.maximum(gmax, sl.silence(parent, t, E[:, :1], params(s, hmax))["Gam_e"])
            keep = E <= gmax[:, None]
            ee, zz = np.nonzero(keep)
            for key, v in (("b", np.full(ee.size, b)), ("tree", np.full(ee.size, ti)), ("e", ee),
                           ("z", zz), ("E", E[ee, zz])):
                recs[key].append(v)
            E_sparse = np.where(keep, E, np.inf)
            prev = {}
            for si, (s, h) in enumerate(SETTINGS):
                sp = params(s, h)
                out = sl.silence(parent, t, E, sp)
                Hn, V = out["H"], out["V"]
                Ht = Hn[:n]
                Hp[b, ti, si] = np.packbits(Ht.reshape(-1))
                nr = np.arange(parent.size - 1)
                viol["dollo_monotone"][0] += int(np.sum(Hn[parent[nr]] & ~Hn[nr]))
                viol["dollo_monotone"][1] += nr.size * K
                viol["size_weighted_V_equals_H"][0] += int(np.sum((sz[:, None] * V).sum(0)
                                                                  != Ht.sum(0)))
                viol["size_weighted_V_equals_H"][1] += K
                Hr = H_from_record(parent, t, n, E_sparse, sp)
                viol["record_rebuilds_H"][0] += int(np.sum(Hr != Ht))
                viol["record_rebuilds_H"][1] += Ht.size
                if h == 0:
                    viol["Hbar0_empty"][0] += int(Hn.sum()); viol["Hbar0_empty"][1] += Hn.size
                if s in prev:
                    viol["nested_across_Hbar"][0] += int(np.sum(prev[s] & ~Ht))
                    viol["nested_across_Hbar"][1] += Ht.size
                prev[s] = Ht
                if ti == 0 and b == 0:
                    E2 = sl.loss_draws(parent.size, K, np.random.default_rng(sE))
                    o2 = sl.silence(parent, t, E2, sp)
                    viol["replay_identical"][0] += int(np.sum(o2["H"] != Hn))
                    viol["replay_identical"][1] += Hn.size
                sil = int(Ht.sum())
                her = int(np.sum((sz[:, None] * V)[sz >= 2])) if sil else 0
                rows.append({"b": b, "tree": ti, "shape": s, "Hbar": h,
                             "Hbar_realised": float(Ht.mean()),
                             "heritable_share_realised": her / sil if sil else None,
                             "heritable_share_expected": sl.heritable_share(parent, t, n, sp),
                             "visible_events_per_tape": float(V.sum() / K),
                             "stem_event_frac": float(V[parent.size - 1].mean())})
    rec = {k_: np.concatenate(v) if v else np.zeros(0) for k_, v in recs.items()}
    fname = OUT / f"sl_n{n}_rho{rho:g}_th{th:g}{'_quick' if quick else ''}.npz"
    np.savez_compressed(
        fname, Hp=Hp, rec_b=rec["b"].astype(np.int8), rec_tree=rec["tree"].astype(np.int8),
        rec_e=rec["e"].astype(np.int32), rec_z=rec["z"].astype(np.int8),
        rec_E=rec["E"].astype(np.float64),
        settings=np.array([f"{s}:{h:g}" for s, h in SETTINGS]),
        meta=json.dumps({"cell": [n, rho, th], "task": task, "seed": SEED, "reps": reps,
                         "trees": ntree, "k": K, "Hbar_max": hmax,
                         "settings": [{"shape": s, "Hbar": h, "knots": SHAPES[s][0],
                                       "levels": SHAPES[s][1], "Gam_T": params(s, h).Gam_T}
                                      for s, h in SETTINGS],
                         "layout": "Hp (b, tree, setting, packbits of H.reshape(n*k)); "
                                   "trees = 14.load_library_trees(n, rho, theta)[:20], as in 19",
                         "compose": "Y = ~X & ~H & (D >= 1) with rung-1 file at same (b, tree)"}))
    manifest = {"cell": {"n": n, "rho": rho, "theta": th, "task": task}, "seed": SEED,
                "rows": rows, "assertions": viol, "record_entries": int(rec["e"].size),
                "file_MB": round(fname.stat().st_size / 1e6, 3),
                "seconds": round(time.perf_counter() - t_start, 1)}
    tag = "_quick" if quick else ""
    (OUT / f"manifest_task{task:02d}{tag}.json").write_text(json.dumps(manifest, indent=1,
                                                                       default=float))
    bad = sum(v[0] for v in viol.values())
    print(f"task {task} {CELLS[task]}: {manifest['seconds']} s, {manifest['file_MB']} MB, "
          f"record {manifest['record_entries']} entries, assertion violations {bad}")
    assert bad == 0, viol


def collect():
    """Aggregate the per-task manifests (sanity summary only; nothing concluded)."""
    mans = [json.loads(p.read_text()) for p in sorted(OUT.glob("manifest_task??.json"))]
    assert len(mans) == len(CELLS), f"{len(mans)} of {len(CELLS)} manifests"
    viol = sum(sum(v[0] for v in m["assertions"].values()) for m in mans)
    table = []
    for m in mans:
        c = m["cell"]
        for s, h in SETTINGS:
            rr = [r for r in m["rows"] if r["shape"] == s and r["Hbar"] == h]
            hr = np.array([r["Hbar_realised"] for r in rr])
            hs = [r["heritable_share_realised"] for r in rr if r["heritable_share_realised"]
                  is not None]
            he = [r["heritable_share_expected"] for r in rr if r["heritable_share_expected"]
                  == r["heritable_share_expected"]]
            table.append({"n": c["n"], "rho": c["rho"], "theta": c["theta"], "shape": s,
                          "Hbar": h, "Hbar_realised_mean": float(hr.mean()),
                          "Hbar_realised_sd": float(hr.std(ddof=1)),
                          "heritable_share_realised_mean": float(np.mean(hs)) if hs else None,
                          "heritable_share_expected_mean": float(np.mean(he)) if he else None,
                          "visible_events_per_tape_per_tree":
                              float(np.mean([r["visible_events_per_tape"] for r in rr])),
                          "stem_event_frac": float(np.mean([r["stem_event_frac"] for r in rr]))})
    out = {"tasks": len(mans), "files": len(mans), "assertion_violations": int(viol),
           "size_MB": round(sum(m["file_MB"] for m in mans), 1),
           "max_task_seconds": max(m["seconds"] for m in mans),
           "trees_x_reps_per_cell": len(mans[0]["rows"]) // len(SETTINGS),
           "settings": [f"{s}:{h:g}" for s, h in SETTINGS], "by_cell_setting": table,
           "note": "manifest sanity check only; nothing concluded. Realised Hbar is the share of "
                   "tip x tape entries silenced (true state, before mask/censoring/filters); "
                   "heritable share = silenced entries in visible events of >= 2 cells / all "
                   "silenced entries."}
    SUMMARY.write_text(json.dumps(out, indent=1, default=float))
    print(f"{len(mans)} manifests, {viol} assertion violations, {out['size_MB']} MB")


if __name__ == "__main__":
    if sys.argv[1] == "--collect":
        collect()
    else:
        run(int(sys.argv[1]), "--quick" in sys.argv)
