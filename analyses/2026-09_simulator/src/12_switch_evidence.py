#!/usr/bin/env python3
r"""
12_switch_evidence.py -- how far can a lineage-specific switch be detected?

Approved by Justin 2026-09-24 ("PROPOSAL (2026-09-24b)" in ../CLAUDE.md), with his two
additions: a clone size n = 8, and capture fractions up to 0.8.  Supersedes the
information-budget pair (temporal rank + per-branch evidence) proposed earlier the same day.

THE EVENT.  A switch: at a point (branch x0, time s0) a cell turns a pathway on, every
descendant inherits the state, and all of them turn it off at s1 = min(s0 + D, 1), with
D = Lam_e / Lam_T the duration as a fraction of the experiment.  Switch locations are
uniform along the tree's total length (a constant per-lineage switch rate), STEM
INCLUDED: recording starts at the clone's founding, s = 0.  Time s = t/T in [0, 1].

THE EVIDENCE.  Expected log-likelihood ratio in favour of the switch when it happened
(nats), per tape per unit (1-d):
    I1_branch = sum_b W_b KL(pbar_b || p0)          edits placed only to their branch
    I1_known  = KL(p1 || p0) sum_b Won_b            edit times known (upper bound)
  W_b    = E(e_b) - E(a_b), expected edits per tape written on branch b = [a_b, e_b]
  E(s)   = E[min(Poisson(Lam_T s), N)] = sum_{j=1..N} P(j, Lam_T s)  (regularised lower gamma)
  Won_b  = the same over the part of b that is on
  pbar_b = p0 + (Won_b / W_b)(p1 - p0)  -- the signal share of an edit known only to sit on b
  KL     = Bernoulli (each edit is the signal symbol or not): ENGRAM moves xi, not lambda,
           so the total count on a branch carries no information about the state.
Full evidence is I = k (1-d) I1, LINEAR in k, so I1 is stored and any k evaluated afterwards.
  ⚠ The branch-only statistic uses the unconditional U-weighting of edit times on a branch;
  the exact branch-level likelihood would also condition on the tape depth at the branch
  start, so I1_branch is a slight LOWER bound on branch-level information.
  ⚠ By convexity of KL in its first argument, I1_branch <= I1_known always -- asserted.

CALLABLE.  I >= 3 + ln H nats, H = (2n-1) * |Lam_e grid| candidate switches (a Bonferroni
charge for searching; conservative, since positions are correlated).  The pseudo-LR has
E_H0[LR] = 1 exactly, per edit, so P_H0(LLR >= c) <= e^-c holds for it.  "Callable" means
expected evidence clears the bar -- roughly 50% power, not certainty.

TREES.  rho in {0.002, 0.1, 0.25} from the library (theta = 0.5); rho in {0.5, 0.8} simulated
fresh here with 04's run_cell (same code, same check B), into results/tree_library_hi_rho/
(gitignored as *.npz).  ⚠ The library stores branch LENGTHS with the root's own branch = 0,
so absolute times are rebuilt from the tips (all at s = 1); the stem is 0 -> t_root.

VERIFICATION FIRST (written to the JSON before any result): E(s) against Monte Carlo, and the
closed-form I1 against the mean of the realised LLR over simulated tapes on a real tree,
both modes, reported in standard errors.

usage:  12_switch_evidence.py [--smoke]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import sys
import time

import numpy as np
from scipy.special import gammainc

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE.parent / "results"
LIB = RES / "tree_library"
HI = RES / "tree_library_hi_rho"

LAM_T = 5.5          # edits per tape at harvest (mice 5.39-5.69, log session 13)
N_SITES = 6
D_DROP = 0.44        # per-tape dropout, Park
THETA = 0.5          # turnover, middle of the 0.46-0.60 growth-arithmetic range
NS = [8, 25, 64, 210, 626, 1976]
RHOS_LIB = [0.002, 0.1, 0.25]
RHOS_NEW = [0.5, 0.8]
LAM_E = np.array([0.25, 0.5, 1.0, 2.0])      # switch duration, edits per tape
DS = LAM_E / LAM_T                            # ... as a fraction of the experiment
FOLDS = np.array([2.0, 5.0, 15.0, 24.0])
P1S = np.array([0.05, 0.1, 0.2, 0.3, 0.4])    # signal share when on
KS = np.array([8, 11, 20, 30, 50, 100, 166])  # tapes per cell evaluated
K_HEAD = 30
S_PER_TREE = 1000
NREP = 20
SEED_NEW = 24_092_026
MB_EDGES = np.array([1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 10 ** 9])
S_EDGES = np.linspace(0.0, 1.0, 11)
MSTAR_MIN_COUNT = 50  # a clade-size bin needs this many switch locations to define m*

P0 = P1S[None, :] / FOLDS[:, None]            # (F, P)
P1 = np.broadcast_to(P1S[None, :], P0.shape)  # (F, P)


def E_edits(s):
    """Expected edits per tape by time s: E[min(Poisson(Lam_T s), N)]."""
    x = LAM_T * np.asarray(s, dtype=np.float64)
    out = np.zeros_like(x)
    for j in range(1, N_SITES + 1):
        out += gammainc(j, x)
    return out


def kl_bern(p, q):
    return p * np.log(p / q) + (1 - p) * np.log((1 - p) / (1 - q))


KL1 = kl_bern(P1, P0)                          # (F, P) nats per edit, timing known


# --------------------------------------------------------------------------
# trees
# --------------------------------------------------------------------------

def prep_tree(parent, branch, n):
    parent = parent.astype(np.int64)
    branch = branch.astype(np.float64)
    nn = parent.size
    root = nn - 1
    assert parent[root] == -1 and nn == 2 * n - 1
    t = np.full(nn, np.nan)
    t[:n] = 1.0
    worst = 0.0
    for x in range(nn - 1):                    # children have smaller ids than parents
        tp = t[x] - branch[x]
        p = parent[x]
        if np.isnan(t[p]):
            t[p] = tp
        else:
            worst = max(worst, abs(t[p] - tp))
    assert worst < 1e-4, f"tree not ultrametric to 1e-4 ({worst:g})"
    assert 0.0 < t[root] < 1.0, f"root time {t[root]}"
    a = np.empty(nn)
    a[:root] = t[parent[:root]]
    a[root] = 0.0
    e = t.copy()
    assert np.all(e - a >= -1e-6)
    kids = [[] for _ in range(nn)]
    for x in range(nn - 1):
        kids[parent[x]].append(x)
    order, stack = [], [root]
    while stack:                               # preorder: every subtree is contiguous
        x = stack.pop()
        order.append(x)
        stack.extend(kids[x])
    order = np.asarray(order)
    pos = np.empty(nn, np.int64)
    pos[order] = np.arange(nn)
    m = np.zeros(nn, np.int64)
    m[:n] = 1
    for x in range(nn - 1):
        m[parent[x]] += m[x]
    return dict(n=n, parent=parent, a=a, e=e, L=np.maximum(e - a, 0.0), order=order,
                pos=pos, m=m, size=2 * m - 1, W=E_edits(e) - E_edits(a), root=root)


def subtree(tr, x0):
    p = tr["pos"][x0]
    return tr["order"][p:p + tr["size"][x0]]


def load_trees(n, rho):
    d = HI if rho in RHOS_NEW else LIB
    files = list(d.glob(f"trees_n{n}_rho{rho:g}_th{THETA:g}_r*.npz"))
    files.sort(key=lambda f: int(re.search(r"_r(\d+)\.npz$", f.name).group(1)))
    par, br = [], []
    for f in files:
        z = np.load(f)
        meta = json.loads(str(z["meta"]))
        assert not meta["faults"], f"{f.name}: {meta['faults']}"
        par.append(z["parent"]); br.append(z["branch"])
    assert par, f"no trees for n={n} rho={rho}"
    par, br = np.concatenate(par), np.concatenate(br)
    assert par.shape[0] >= NREP, f"n={n} rho={rho}: only {par.shape[0]} trees"
    return [prep_tree(par[i], br[i], n) for i in range(NREP)]


def make_hi_rho_trees(ns):
    spec = importlib.util.spec_from_file_location("treelib", HERE / "04_tree_library.py")
    lib = importlib.util.module_from_spec(spec)
    sys.modules["treelib"] = lib
    spec.loader.exec_module(lib)
    HI.mkdir(exist_ok=True)
    made = []
    for rho in RHOS_NEW:
        for n in ns:
            f = HI / f"trees_n{n}_rho{rho:g}_th{THETA:g}_r0.npz"
            if f.exists():
                continue
            seed0 = SEED_NEW + 7919 * n + int(round(1000 * rho))
            _, faults, sec = lib.run_cell(n, rho, THETA, 0, NREP, seed0, HI)
            assert not faults, f"n={n} rho={rho}: {faults}"
            made.append({"n": n, "rho": rho, "seed0": seed0, "seconds": round(sec, 1)})
            print(f"  simulated n={n} rho={rho}: {sec:.1f} s", flush=True)
    return made


# --------------------------------------------------------------------------
# evidence for one switch
# --------------------------------------------------------------------------

def switch_evidence(tr, x0, s0):
    """(I1_known, I1_branch), each (D, F, P), nats per tape per unit (1-d)."""
    sub = subtree(tr, x0)
    lo = tr["a"][sub].copy()
    lo[0] = s0
    s1 = np.minimum(s0 + DS, 1.0)                               # (D,)
    keep = lo < s1.max()
    lo, e, W = lo[keep], tr["e"][sub][keep], tr["W"][sub][keep]
    hi = np.minimum(e[None, :], s1[:, None])                    # (D, L)
    on = hi > lo[None, :]
    cols = on.any(0)
    lo, hi, on, W = lo[cols], hi[:, cols], on[:, cols], W[cols]
    Won = np.where(on, E_edits(hi) - E_edits(np.broadcast_to(lo, hi.shape)), 0.0)
    Won = np.clip(Won, 0.0, None)
    known = Won.sum(1)[:, None, None] * KL1[None]              # (D, F, P)
    f = np.divide(Won, W[None, :], out=np.zeros_like(Won), where=W[None, :] > 0)
    f = np.minimum(f, 1.0)
    pbar = P0[None, None] + f[:, :, None, None] * (P1 - P0)[None, None]   # (D, L, F, P)
    branch = (W[None, :, None, None] * kl_bern(pbar, P0[None, None])).sum(1)
    return known, branch


def sample_switches(tr, S, rng):
    cum = np.cumsum(tr["L"])
    tot = cum[-1]
    u = (np.arange(S) + rng.random()) / S * tot                 # systematic, uniform in length
    x0 = np.minimum(np.searchsorted(cum, u, side="right"), tr["L"].size - 1)
    s0 = np.clip(tr["e"][x0] - (cum[x0] - u), tr["a"][x0], tr["e"][x0])
    return x0, s0, tot


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------

def check_E(rng, R=1_000_000):
    rows = []
    for s in (0.05, 0.3, 0.6, 1.0):
        x = np.minimum(rng.poisson(LAM_T * s, R), N_SITES)
        mc, se = x.mean(), x.std() / np.sqrt(R)
        cf = float(E_edits(s))
        rows.append({"s": s, "closed_form": cf, "mc": float(mc), "z": float((mc - cf) / se)})
    return rows


def mc_llr(tr, x0, s0, D, F, p1, R, rng):
    """Realised LLR over R simulated tapes (k = 1, d = 0) vs the closed form."""
    p0 = p1 / F
    s1 = min(s0 + D, 1.0)
    nn = tr["a"].size
    in_sub = np.zeros(nn, bool)
    in_sub[subtree(tr, x0)] = True
    lo = tr["a"].copy()
    lo[x0] = s0
    hi = np.minimum(tr["e"], s1)
    on = in_sub & (hi > lo)
    Won = np.where(on, E_edits(hi) - E_edits(lo), 0.0).clip(0)
    f = np.divide(Won, tr["W"], out=np.zeros(nn), where=tr["W"] > 0).clip(0, 1)
    pbar = p0 + f * (p1 - p0)
    cf_known = float(kl_bern(p1, p0) * Won.sum())
    cf_branch = float((tr["W"] * kl_bern(pbar, p0)).sum())

    c_end = np.zeros((nn, R), np.int16)
    llr_k, llr_b = np.zeros(R), np.zeros(R)
    for x in tr["order"]:
        c = np.zeros(R, np.int16) if x == tr["root"] else c_end[tr["parent"][x]].copy()
        tcur = np.full(R, tr["a"][x])
        active = c < N_SITES
        start_on = s0 if x == x0 else tr["a"][x]
        for _ in range(N_SITES):
            tcur = tcur + rng.exponential(1.0 / LAM_T, R)
            ev = active & (tcur < tr["e"][x])
            if not ev.any():
                break
            is_on = in_sub[x] & (tcur >= start_on) & (tcur < s1)
            p = np.where(is_on, p1, p0)
            X = rng.random(R) < p
            llr_k += ev * np.where(X, np.log(p / p0), np.log((1 - p) / (1 - p0)))
            pb = pbar[x]
            llr_b += ev * np.where(X, np.log(pb / p0), np.log((1 - pb) / (1 - p0)))
            c = c + ev
            active = ev & (c < N_SITES)
        c_end[x] = c
    out = {}
    for name, v, cf in (("known", llr_k, cf_known), ("branch", llr_b, cf_branch)):
        mc, se = float(v.mean()), float(v.std() / np.sqrt(R))
        out[name] = {"closed_form": cf, "mc": mc, "se": se,
                     "z": (mc - cf) / se if se > 0 else 0.0}
    return out


def verify(rng, smoke):
    tr = load_trees(25, 0.25)[0]
    root = tr["root"]
    internal = [x for x in range(tr["n"], root) if 3 <= tr["m"][x] <= 12]
    terminal = int(np.argmax(tr["L"][:tr["n"]]))
    sites = [("stem", root, 0.5 * tr["e"][root]),
             ("internal", internal[0], 0.5 * (tr["a"][internal[0]] + tr["e"][internal[0]])),
             ("terminal", terminal, 0.5 * (tr["a"][terminal] + tr["e"][terminal]))]
    R = 4000 if smoke else 40000
    rows = []
    for label, x0, s0 in sites:
        for lam_e in (0.5, 2.0):
            for F, p1 in ((24.0, 0.3), (2.0, 0.3), (5.0, 0.05)):
                r = mc_llr(tr, x0, s0, lam_e / LAM_T, F, p1, R, rng)
                rows.append({"where": label, "m": int(tr["m"][x0]), "lam_e": lam_e, "F": F,
                             "p1": p1, **{f"{k}_{q}": v for k, d in r.items()
                                          for q, v in d.items()}})
    zmax = max(max(abs(r["known_z"]), abs(r["branch_z"])) for r in rows)
    return {"E_edits": check_E(rng), "llr": rows, "R_tapes": R, "max_abs_z_llr": zmax}


# --------------------------------------------------------------------------
# one (n, rho) cell
# --------------------------------------------------------------------------

def _clean(a, nd=4):
    a = np.asarray(a, dtype=np.float64)
    return [None if not np.isfinite(v) else round(float(v), nd) for v in a.ravel()] \
        if a.ndim == 1 else [_clean(x, nd) for x in a]


def run_cell(n, rho, S, rng):
    trees = load_trees(n, rho)
    thr = 3.0 + np.log((2 * n - 1) * LAM_E.size)
    I1, M, S0, TID = [], [], [], []
    tot_len, term_share, stem = [], [], []
    for i, tr in enumerate(trees):
        x0s, s0s, tot = sample_switches(tr, S, rng)
        tot_len.append(tot)
        term_share.append(tr["L"][:n].sum() / tot)
        stem.append(tr["e"][tr["root"]])
        for x0, s0 in zip(x0s, s0s):
            k_, b_ = switch_evidence(tr, int(x0), float(s0))
            assert np.all(b_ <= k_ * (1 + 1e-9) + 1e-12), "branch-only exceeded known timing"
            I1.append(np.stack([k_, b_]))                         # (2, D, F, P)
        M.append(tr["m"][x0s]); S0.append(s0s); TID.append(np.full(S, i))
    I1 = np.asarray(I1)                                            # (Nsw, 2, D, F, P)
    M, S0, TID = np.concatenate(M), np.concatenate(S0), np.concatenate(TID)

    with np.errstate(divide="ignore"):
        kneed = thr / ((1 - D_DROP) * I1)                          # tapes per cell needed
        kneed3 = 3.0 / ((1 - D_DROP) * I1)
    call = kneed[..., None] <= KS                                  # (Nsw, 2, D, F, P, K)
    frac_k = call.mean(0)
    head = kneed <= K_HEAD                                         # (Nsw, 2, D, F, P)
    per_tree = np.stack([head[TID == i].mean(0) for i in range(len(trees))])
    mb = np.searchsorted(MB_EDGES, M, side="right") - 1
    sb = np.clip(np.searchsorted(S_EDGES, S0, side="right") - 1, 0, S_EDGES.size - 2)
    nmb, nsb = MB_EDGES.size - 1, S_EDGES.size - 1
    mcount = np.bincount(mb, minlength=nmb)
    scount = np.bincount(sb, minlength=nsb)
    by_m = np.full((nmb,) + head.shape[1:], np.nan)
    by_s = np.full((nsb,) + head.shape[1:], np.nan)
    for j in range(nmb):
        if mcount[j]:
            by_m[j] = head[mb == j].mean(0)
    for j in range(nsb):
        if scount[j]:
            by_s[j] = head[sb == j].mean(0)
    ok = (mcount >= MSTAR_MIN_COUNT)[:, None, None, None, None] & (by_m >= 0.5)
    first = np.where(ok.any(0), ok.argmax(0), -1)
    m_star = np.where(first >= 0, MB_EDGES[np.maximum(first, 0)], np.inf)
    return {
        "n": n, "rho": rho, "threshold_nats": round(float(thr), 4),
        "n_switches": int(M.size), "mean_tree_length_T": round(float(np.mean(tot_len)), 4),
        "terminal_length_share": round(float(np.mean(term_share)), 4),
        "mean_stem_T": round(float(np.mean(stem)), 4),
        "frac_callable_by_k": _clean(np.moveaxis(frac_k, -1, 0)),       # [K][mode][D][F][P]
        "frac_callable_k30_tree_mean": _clean(per_tree.mean(0)),        # [mode][D][F][P]
        "frac_callable_k30_tree_sd": _clean(per_tree.std(0, ddof=1)),
        "frac_callable_k30_at_3nats": _clean((kneed3 <= K_HEAD).mean(0)),
        "k50_tapes": _clean(np.median(kneed, axis=0), 1),               # [mode][D][F][P]
        "m_star_k30": _clean(m_star, 0),                                # [mode][D][F][P]
        "clade_bins_left_edge": MB_EDGES[:-1].tolist(),
        "clade_bin_share": _clean(mcount / M.size),
        "frac_callable_k30_by_clade": _clean(by_m),                     # [bin][mode][D][F][P]
        "time_bins_left_edge": S_EDGES[:-1].round(2).tolist(),
        "time_bin_share": _clean(scount / M.size),
        "frac_callable_k30_by_time": _clean(by_s),                      # [bin][mode][D][F][P]
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    ns = [8, 25] if args.smoke else NS
    rhos = [0.25, 0.5] if args.smoke else RHOS_LIB + RHOS_NEW
    S = 100 if args.smoke else S_PER_TREE
    rng = np.random.default_rng(20260924)
    t0 = time.perf_counter()

    out = {"settings": {
        "Lam_T_edits_per_tape": LAM_T, "N_sites": N_SITES, "dropout_d": D_DROP, "theta": THETA,
        "n": ns, "rho": rhos, "lam_e_edits_per_tape": LAM_E.tolist(), "fold_F": FOLDS.tolist(),
        "p1": P1S.tolist(), "k_grid": KS.tolist(), "k_head": K_HEAD,
        "switches_per_tree": S, "trees_per_cell": NREP,
        "threshold": "3 + ln((2n-1) * |lam_e grid|) nats",
        "axes": {"mode": ["known_timing", "branch_only"], "D": "lam_e", "F": "fold_F",
                 "P": "p1", "K": "k_grid"}}}

    print("verification ...", flush=True)
    out["verification"] = verify(rng, args.smoke)
    v = out["verification"]
    print(f"  E(s) z: {[round(r['z'], 2) for r in v['E_edits']]}", flush=True)
    print(f"  LLR closed form vs MC: max |z| = {v['max_abs_z_llr']:.2f} over "
          f"{len(v['llr'])} x 2 checks", flush=True)

    print("hi-rho trees ...", flush=True)
    out["hi_rho_trees_simulated"] = make_hi_rho_trees(ns)

    out["cells"] = []
    for rho in rhos:
        for n in ns:
            tc = time.perf_counter()
            c = run_cell(n, rho, S, rng)
            out["cells"].append(c)
            # headline, for the log: branch-only, lam_e 0.5, F 24, p1 0.3, k 30
            d, f, p = 1, 3, 3
            print(f"  n={n:5d} rho={rho:<6g} thr={c['threshold_nats']:.2f}  "
                  f"callable(k30, 24x, p1 .3, lam_e .5): known "
                  f"{c['frac_callable_k30_tree_mean'][0][d][f][p]:.3f}  branch "
                  f"{c['frac_callable_k30_tree_mean'][1][d][f][p]:.3f}  "
                  f"[{time.perf_counter() - tc:.0f} s]", flush=True)

    out["seconds"] = round(time.perf_counter() - t0, 1)
    name = "switch_evidence_smoke.json" if args.smoke else "switch_evidence.json"
    (RES / name).write_text(json.dumps(out))
    print(f"wrote results/{name} in {out['seconds']} s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
