#!/usr/bin/env python3
r"""
19_dropout_library.py -- the rung-1 dropout library: library trees x tape panels x masks.

Approved 2026-09-30 ("PROPOSAL (2026-09-29)", amendment 6, ../CLAUDE.md); generated only after
18_validate_dropout.py PASSED.  One array task per tree cell (n, rho, theta).  NOTHING is concluded
from the library; the only summary is a manifest of realised vs target technical missingness.

THREE LAYERS, so costly work is reused.
  A  trees     20 stored library trees per cell (no new trees -- rung 1 is observation-only)
  B  panel +   per (Lam_T, panel type, replicate b): a tape panel (u, r) and one decoration of each
     editing   tree by 14_editing.py.  A panel is SHARED by every tree of the library at the same
               b, as all clones of one experiment share their integrations.
  C  masks     30 technical masks per decorated tree, from per-tree common random numbers
               (a, V, W), so masks are nested as pbar rises and comparable across settings.

GRID.  Lam_T in {2.9, 5.5} edits/tape; panel types none (phi 0) / park (phi .15, r_cl .6) /
strong (phi .3, r_cl .4), all with s_r 0.3, Delta_beta 2; 5 replicates; masks
pbar {0, .1, .25, .44, .6} x s_beta {0, 1, 2} x s_alpha {0, .7}.  k = 30, N = 6, lineage channel
only (Park mice pooled xi^L), constant rate, Lam_pre = 0.

SEEDS.  Panel draws: (SEED, 1, b) -- identical across Lam_T and panel type.  Decoration:
(SEED, 2, cell, tree, b).  Mask draws: (SEED, 3, cell, tree, b).  Neither depends on Lam_T or
panel type, so the library is replayable from the stored settings alone.

STORED, one file per (cell, Lam_T, panel type), axes (b, tree, ...):
  D    (5, 20, n, k) int8        true tip depth
  sym  (5, 20, n, k, N) int16    true tip symbols (-1 empty)
  X    (5, 20, 30, ceil(n k / 8)) uint8   technical masks, np.packbits of X.reshape(n*k)
  a, W (5, 20, n) float32        cell capture draw (alpha = s_alpha a) and ClonalBC uniform
  panel: U, eps, v, u, r (5, k); beta (5, 30, k); mu_beta (5, 30)
Observed view: Y = ~X & (D >= 1) (censoring on) or ~X (off).  Filters are applied at export
(17.export_arm), never stored.

usage:  19_dropout_library.py <task 0..29> [--quick]
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
import time

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("dropout", _HERE / "17_dropout.py")
dl = importlib.util.module_from_spec(_spec)
sys.modules["dropout"] = dl
_spec.loader.exec_module(dl)
ed = dl.ed

OUT = _HERE.parent / "results" / "dropout_library"
SEED = 30092026
K, N, TREES, REPS = 30, 6, 20, 5
CELLS = [(n, rho, th) for n in (25, 64, 210, 626, 1976) for rho in (0.002, 0.1, 0.25)
         for th in (0.3, 0.7)]
LAMS = (2.9, 5.5)
PANELS = {"none": dl.PanelParams(0.0, 1.0, 0.3, 2.0),
          "park": dl.PanelParams(0.15, 0.6, 0.3, 2.0),
          "strong": dl.PanelParams(0.3, 0.4, 0.3, 2.0)}
MASKS = [dl.MaskParams(p, sb, sa) for p in (0.0, 0.1, 0.25, 0.44, 0.6)
         for sb in (0.0, 1.0, 2.0) for sa in (0.0, 0.7)]


def main():
    task = int(sys.argv[1])
    quick = "--quick" in sys.argv
    n, rho, th = CELLS[task]
    reps, ntree = (1, 2) if quick else (REPS, TREES)
    trees = ed.load_library_trees(n, rho, th)
    assert len(trees) >= ntree, f"cell {CELLS[task]} has {len(trees)} trees"
    trees = trees[:ntree]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    panels = [dl.panel_draws(K, np.random.default_rng(np.random.SeedSequence([SEED, 1, b])))
              for b in range(reps)]
    mdraws = [[dl.mask_draws(n, K, np.random.default_rng(
        np.random.SeedSequence([SEED, 3, task, ti, b]))) for ti in range(ntree)]
        for b in range(reps)]
    manifest = {"cell": {"n": n, "rho": rho, "theta": th, "task": task}, "seed": SEED,
                "files": [], "assertions": {}}
    viol = {"Y0_where_D0": [0, 0], "replay_identical": [0, 0], "depth_le_N": [0, 0]}
    nb = (n * K + 7) // 8
    for lam in LAMS:
        for pname, pp in PANELS.items():
            D = np.zeros((reps, ntree, n, K), np.int8)
            S = np.full((reps, ntree, n, K, N), -1, np.int16)
            Xp = np.zeros((reps, ntree, len(MASKS), nb), np.uint8)
            U = np.zeros((reps, K)); EPS = np.zeros((reps, K)); V = np.zeros((reps, K))
            UU = np.zeros((reps, K), np.int8); R = np.zeros((reps, K))
            BETA = np.zeros((reps, len(MASKS), K)); MU = np.zeros((reps, len(MASKS)))
            for b in range(reps):
                pd_ = panels[b]
                u, r = dl.panel_tapes(pd_, pp)
                U[b], EPS[b], V[b], UU[b], R[b] = pd_["U"], pd_["eps"], pd_["v"], u, r
                for mi, mp in enumerate(MASKS):
                    BETA[b, mi], MU[b, mi] = dl.panel_beta(pd_, u, pp, mp)
                rec = ed.Recorder(k=K, Lam_T=lam, r=r, shares=((0.0,),))
                for ti, (parent, branch) in enumerate(trees):
                    sdec = np.random.SeedSequence([SEED, 2, task, ti, b])
                    d_, s_, _, _ = dl.decorate_tips(parent, branch, n, rec,
                                                    np.random.default_rng(sdec))
                    D[b, ti], S[b, ti] = d_, s_
                    viol["depth_le_N"][0] += int(np.sum(d_ > N)); viol["depth_le_N"][1] += d_.size
                    if ti == 0:
                        d2, s2, _, _ = dl.decorate_tips(parent, branch, n, rec,
                                                        np.random.default_rng(sdec))
                        viol["replay_identical"][0] += int(np.sum(d2 != d_) + np.sum(s2 != s_))
                        viol["replay_identical"][1] += d_.size * (1 + N)
                    md = mdraws[b][ti]
                    for mi, mp in enumerate(MASKS):
                        X, _ = dl.technical_mask(md, BETA[b, mi], mp.s_alpha)
                        Xp[b, ti, mi] = np.packbits(X.reshape(-1))
                        Y = dl.observe(X, d_)
                        viol["Y0_where_D0"][0] += int(np.sum(Y & (d_ == 0)))
                        viol["Y0_where_D0"][1] += int((d_ == 0).sum())
                # manifest: realised vs target, pooled over this replicate's trees
                for mi, mp in enumerate(MASKS):
                    Xall = np.stack([np.unpackbits(Xp[b, ti, mi])[:n * K].reshape(n, K)
                                     for ti in range(ntree)]).astype(bool)
                    Yall = (~Xall) & (D[b] >= 1)
                    manifest["files"].append({
                        "Lam_T": lam, "panel": pname, "rep": b, "mask": dl.asdict(mp),
                        "closed_tapes": int(u.sum()), "mean_r": float(r.mean()),
                        "realised_technical_missing": float(Xall.mean()),
                        "observed_missing": float((~Yall).mean()),
                        "unedited_share": float((D[b] == 0).mean())})
            fname = OUT / f"dl_n{n}_rho{rho:g}_th{th:g}_L{lam:g}_{pname}{'_quick' if quick else ''}.npz"
            np.savez_compressed(
                fname, D=D, sym=S, X=Xp,
                a=np.stack([[m["a"] for m in mb] for mb in mdraws]).astype(np.float32),
                W=np.stack([[m["W"] for m in mb] for mb in mdraws]).astype(np.float32),
                U=U, eps=EPS, v=V, u=UU, r=R, beta=BETA, mu_beta=MU,
                meta=dl.settings_json(cell=[n, rho, th], task=task, Lam_T=lam, panel=pp,
                                      panel_name=pname, masks=[dl.asdict(m) for m in MASKS],
                                      seed=SEED, reps=reps, trees=ntree, k=K, N=N,
                                      M=ed.Recorder(k=K, shares=((0.0,),)).M,
                                      layout="(b, tree, ...); X packbits of X.reshape(n*k)"))
            print(f"{fname.name}: {time.perf_counter() - t0:.0f} s", flush=True)
    manifest["assertions"] = viol
    manifest["seconds"] = round(time.perf_counter() - t0, 1)
    tag = "_quick" if quick else ""
    (OUT / f"manifest_task{task:02d}{tag}.json").write_text(json.dumps(manifest, indent=1,
                                                                       default=float))
    bad = sum(v[0] for v in viol.values())
    print(f"task {task} {CELLS[task]}: {manifest['seconds']} s, assertion violations {bad}")
    assert bad == 0, viol


if __name__ == "__main__":
    main()
