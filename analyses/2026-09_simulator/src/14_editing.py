#!/usr/bin/env python3
r"""
14_editing.py -- the editing layer: decorate a stored tree with sequential Typewriter edits.

Approved by Justin 2026-09-27 ("PROPOSAL (2026-09-24c)" in ../CLAUDE.md) with Lambda_pre = 0:
tapes are EMPTY at the clone founder (Park to be revisited). Theory: notes/sciphy_notes.md
§S4.9. Validation: 15_validate_editing.py. Nothing scientific is concluded here.

THE MODEL (SciPhy's recorder, §1a).  A tape has N slots filled left to right.  While it has an
open slot it is edited as a Poisson process at rate lambda_z(t); each edit writes one symbol,
drawn from xi(t), into the next open slot; at N it is full.  At a division both daughters copy
the parent's tapes.  Given the tree, tapes are independent (row A5).

NOTATION (§S4.9).  t in [0, 1] is calendar time as a fraction of the experiment (T = 1 in the
library).  A branch is e (NEVER b -- b is the birth rate); branch e is the edge INTO node e and
spans [t0_e, t1_e], the root's edge being the stem [0, t_root].  z = tape, i = symbol, j = slot.
  lambda_z(t) = r_z * Lam_T * lambda_0(t)     edits / tape / experiment
  Lam_z(t)    = r_z * Lam_T * W(t),  W(t) = int_0^t lambda_0      (Lambda-time, edits / tape)
lambda_0 is piecewise constant: level omega_g on [kappa_{g-1}, kappa_g), normalised so W(1) = 1.
Composition: A signal channels (share p_a(t), own composition xi^(a)) plus the lineage channel
(share 1 - sum_a p_a(t), composition xi^L).  Each symbol is in exactly one channel.  Symbol ids
are laid out signal channels first, lineage last.  Shares are piecewise constant on their OWN
knots, independent of the rate knots.

THE SAMPLER -- route (ii), §S4.3.  For one branch, vectorised over tapes: a tape arriving at
depth d has c = N - d open slots.  Draw waits E ~ Exp(1) in Lambda-time, S = cumsum(E); the
edits are the S_nu <= mu_{e,z} = Lam_z(t1) - Lam_z(t0) among the first c; each maps back to
calendar time by the exact two-step inverse of §S4.9.2; its symbol is drawn channel-then-symbol
at that time.  ⚑ For vectorisation a fixed (k, N) block of waits is drawn and only the first c
of each row are used, so the law is route (ii)'s exactly; only RNG consumption differs.
Route (i) -- count ~ Poisson(mu), uniforms on [0, mu], sort, keep the first c -- is the
independent second implementation, selected by route="i".

STORAGE -- per BRANCH, the truth written once (§S4.9, decision "edit times stored").
  bcount[e, z]      edits written on branch e, tape z
  bslot0[e, z]      the slot (0-based) the generator believed the branch's first edit went into,
                    i.e. its handed-off parent depth -- stored so the validator can check it
                    against the depth it recomputes independently from bcount along the path
  bsym[e, z, :]     symbol ids of those edits, -1 padded      (int16)
  btime[e, z, :]    their calendar times, NaN padded           (float32)
A cell's tape is the branch edits along its root-to-tip path, joined in order (assemble()).
"""
from __future__ import annotations

import json
import pathlib
from dataclasses import dataclass, field

import numpy as np

_HERE = pathlib.Path(__file__).resolve().parent
RES = _HERE.parent / "results"
LIB = RES / "tree_library"
XI_FILE = _HERE.parents[1] / "2026-08_park-compatibility" / "results" / "xi_vectors.json"
MICE = ("Mouse1", "Mouse2", "Mouse3")


def park_xi_L() -> np.ndarray:
    """Park's lineage composition: the three mouse arms' trimmed xi pooled by edit count.

    Only the frequency VALUES are returned (sorted, descending) -- never the symbol sequences.
    Mice because Lam_T = 5.5 is the mouse editing depth (5.39-5.69, log session 13).
    """
    d = json.loads(XI_FILE.read_text())
    tot: dict[str, float] = {}
    wsum = 0.0
    for arm in MICE:
        w = float(d[arm]["n_edits"])
        wsum += w
        for sym, v in d[arm]["xi"].items():
            tot[sym] = tot.get(sym, 0.0) + w * float(v)
    xi = np.sort(np.array(list(tot.values())) / wsum)[::-1]
    return xi / xi.sum()


@dataclass
class Recorder:
    k: int = 30
    N: int = 6
    Lam_T: float = 5.5
    r: np.ndarray | None = None                   # (k,) tape speeds, mean 1; None -> all 1
    rate_knots: tuple = (0.0, 1.0)                 # kappa_0 .. kappa_G
    rate_levels: tuple = (1.0,)                    # omega_1 .. omega_G (normalised here)
    comp_knots: tuple = (0.0, 1.0)                 # knots of the share schedule
    shares: tuple = ((0.05,),)                     # [epoch][channel] -> p_a
    signal_xi: tuple = ((1.0,),)                   # per channel: its within-channel composition
    xi_L: np.ndarray | None = None                 # lineage composition; None -> Park mice pooled
    _cache: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.r = np.ones(self.k) if self.r is None else np.asarray(self.r, float)
        assert self.r.shape == (self.k,) and np.all(self.r > 0)
        kap = np.asarray(self.rate_knots, float)
        om = np.asarray(self.rate_levels, float)
        assert kap[0] == 0.0 and kap[-1] == 1.0 and np.all(np.diff(kap) > 0)
        assert om.size == kap.size - 1 and np.all(om >= 0) and om.sum() > 0
        om = om / np.sum(om * np.diff(kap))       # int_0^1 lambda_0 = 1
        self.kappa, self.omega = kap, om
        self.W = np.concatenate([[0.0], np.cumsum(om * np.diff(kap))])
        self.W[-1] = 1.0
        # last epoch with positive rate: where tau = Lam_z(1) exactly maps to
        self._last_pos = int(np.nonzero(om > 0)[0][-1])

        ck = np.asarray(self.comp_knots, float)
        sh = np.asarray(self.shares, float)
        assert ck[0] == 0.0 and ck[-1] == 1.0 and np.all(np.diff(ck) > 0)
        assert sh.ndim == 2 and sh.shape[0] == ck.size - 1
        assert np.all(sh >= 0) and np.all(sh.sum(1) <= 1.0)
        self.ck, self.sh = ck, sh
        self.A = sh.shape[1]
        sx = [np.asarray(x, float) / np.sum(x) for x in self.signal_xi]
        assert len(sx) == self.A
        self.xi_L = park_xi_L() if self.xi_L is None else np.asarray(self.xi_L, float)
        self.xi_L = self.xi_L / self.xi_L.sum()
        self.chan_xi = sx + [self.xi_L]            # channel A = lineage
        sizes = [x.size for x in self.chan_xi]
        self.offsets = np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(int)
        self.M = int(np.sum(sizes))
        self.chan_cdf = [np.cumsum(x) for x in self.chan_xi]
        self.q_chan = np.array([np.sum(x ** 2) for x in self.chan_xi])   # q_a ..., q_L
        self.channel_of = np.concatenate([np.full(s, a) for a, s in enumerate(sizes)])

    # ---- the clock ---------------------------------------------------------------------
    def W_of(self, t):
        """W(t) = int_0^t lambda_0: straight line inside each rate epoch."""
        t = np.asarray(t, float)
        g = np.clip(np.searchsorted(self.kappa, t, side="right"), 1, self.omega.size)
        return self.W[g - 1] + self.omega[g - 1] * (t - self.kappa[g - 1])

    def Lam(self, t):
        """Lam_z(t) for every tape: shape t.shape + (k,) ... broadcast as (k,) for scalar t."""
        return self.Lam_T * self.r * self.W_of(t)[..., None] if np.ndim(t) else \
            self.Lam_T * self.r * float(self.W_of(t))

    def W_inv(self, w):
        """Exact inverse of W (§S4.9.2).  Zero-rate epochs occupy an empty interval of W and
        are never selected; the boundary tie resolves into the NEXT active epoch (side='right')."""
        w = np.asarray(w, float)
        g = np.clip(np.searchsorted(self.W, w, side="right"), 1, self.omega.size)
        g = np.where(self.omega[g - 1] > 0, g, self._last_pos + 1)   # only at w >= W(1)
        return self.kappa[g - 1] + (w - self.W[g - 1]) / self.omega[g - 1]

    def Lam_inv(self, tau, z=None):
        """Calendar time at which tape z (all tapes if None, last axis) reached Lambda-time tau."""
        r = self.r if z is None else self.r[z]
        return self.W_inv(np.asarray(tau, float) / (self.Lam_T * r))

    # ---- composition ---------------------------------------------------------------------
    def share_at(self, t):
        """Signal shares p_a(t), shape t.shape + (A,)."""
        h = np.clip(np.searchsorted(self.ck, np.asarray(t, float), side="right"), 1,
                    self.sh.shape[0]) - 1
        return self.sh[h]

    def inner(self, t1, t2):
        """<xi(t1), xi(t2)> = sum_a p_a(t1)p_a(t2)q_a + (1-sum p(t1))(1-sum p(t2)) q_L (§S4.9.4)."""
        p1, p2 = self.share_at(t1), self.share_at(t2)
        sig = np.sum(p1 * p2 * self.q_chan[:-1], axis=-1)
        return sig + (1 - p1.sum(-1)) * (1 - p2.sum(-1)) * self.q_chan[-1]

    def draw_symbols(self, t, rng):
        """Channel first (probabilities p_a(t), then 1 - sum p), symbol second."""
        t = np.asarray(t, float)
        if t.size == 0:
            return np.zeros(0, np.int16)
        cum = np.cumsum(self.share_at(t), axis=-1)                 # (n, A)
        u = rng.random(t.size)
        ch = np.sum(u[:, None] >= cum, axis=1)                     # 0..A, A = lineage
        out = np.empty(t.size, np.int64)
        v = rng.random(t.size)
        for a in range(self.A + 1):
            m = ch == a
            if m.any():
                idx = np.searchsorted(self.chan_cdf[a], v[m], side="right")
                out[m] = self.offsets[a] + np.minimum(idx, self.chan_xi[a].size - 1)
        return out.astype(np.int16)

    def settings(self) -> dict:
        return {"k": self.k, "N": self.N, "Lam_T": self.Lam_T, "r": self.r.round(6).tolist(),
                "rate_knots": self.kappa.tolist(), "rate_levels": self.omega.tolist(),
                "comp_knots": self.ck.tolist(), "shares": self.sh.tolist(),
                "channel_sizes": [int(x.size) for x in self.chan_xi],
                "q_channel": self.q_chan.tolist(), "M": self.M, "Lam_pre": 0.0}


# --------------------------------------------------------------------------
# trees
# --------------------------------------------------------------------------

def node_times(parent, branch, n):
    """Absolute node times with tips at t = 1 (the library stores the root's branch as 0)."""
    parent = np.asarray(parent, np.int64)
    branch = np.asarray(branch, np.float64)
    nn = parent.size
    t = np.full(nn, np.nan)
    t[:n] = 1.0
    worst = 0.0
    for x in range(nn - 1):                          # children have smaller ids than parents
        tp = t[x] - branch[x]
        p = parent[x]
        if np.isnan(t[p]):
            t[p] = tp
        else:
            worst = max(worst, abs(t[p] - tp))
    assert worst < 1e-4, f"tree not ultrametric to 1e-4 ({worst:g})"
    assert 0.0 < t[nn - 1] < 1.0 and parent[nn - 1] == -1
    return t


def load_library_trees(n, rho, theta):
    """All stored trees of one library cell, as a list of (parent, branch)."""
    files = sorted(LIB.glob(f"trees_n{n}_rho{rho:g}_th{theta:g}_r*.npz"))
    assert files, f"no library cell n={n} rho={rho} theta={theta}"
    out = []
    for f in files:
        z = np.load(f)
        meta = json.loads(str(z["meta"]))
        assert not meta["faults"], f"{f.name}: {meta['faults']}"
        out += [(z["parent"][i].astype(np.int64), z["branch"][i].astype(np.float64))
                for i in range(z["parent"].shape[0])]
    return out


# --------------------------------------------------------------------------
# the sampler
# --------------------------------------------------------------------------

def _waits_route_ii(mu, c, N, rng):
    """Lambda-time offsets of the edits on one branch, per tape: (k, N) with +inf = no edit."""
    S = np.cumsum(rng.exponential(1.0, size=(mu.size, N)), axis=1)
    use = (np.arange(N)[None, :] < c[:, None]) & (S <= mu[:, None])
    return np.where(use, S, np.inf)


def _waits_route_i(mu, c, N, rng):
    """Route (i): count ~ Poisson(mu), that many uniforms on [0, mu], sort, keep the first c."""
    k = mu.size
    n = rng.poisson(mu)
    m = max(int(n.max()), N)
    U = rng.random((k, m)) * mu[:, None]
    U = np.where(np.arange(m)[None, :] < n[:, None], U, np.inf)
    S = np.sort(U, axis=1)[:, :N]
    keep = np.arange(N)[None, :] < np.minimum(n, c)[:, None]
    return np.where(keep, S, np.inf)


def decorate(parent, branch, n, rec: Recorder, rng, route="ii"):
    """Write edits onto every branch of one tree.  Parents before children (descending id)."""
    t = node_times(parent, branch, n)
    nn, k, N = parent.size, rec.k, rec.N
    root = nn - 1
    t0 = np.where(np.arange(nn) == root, 0.0, t[np.maximum(parent, 0)])
    depth = np.zeros((nn, k), np.int16)                 # the generator's own hand-off state
    bcount = np.zeros((nn, k), np.int8)
    bslot0 = np.zeros((nn, k), np.int8)
    bsym = np.full((nn, k, N), -1, np.int16)
    btime = np.full((nn, k, N), np.nan, np.float32)
    waits = _waits_route_ii if route == "ii" else _waits_route_i
    for e in range(nn - 1, -1, -1):
        d = np.zeros(k, np.int16) if e == root else depth[parent[e]]
        c = N - d
        L0 = rec.Lam(t0[e])                             # (k,)
        mu = rec.Lam(t[e]) - L0
        S = waits(mu, c, N, rng)                        # (k, N), prefix-shaped
        used = np.isfinite(S)
        ne = used.sum(1)
        if ne.any():
            zz, jj = np.nonzero(used)
            te = rec.Lam_inv(L0[zz] + S[zz, jj], zz)
            te = np.clip(te, t0[e], t[e])
            bsym[e, zz, jj] = rec.draw_symbols(te, rng)
            btime[e, zz, jj] = te
        bcount[e] = ne
        bslot0[e] = d
        depth[e] = d + ne
    return {"n": n, "parent": np.asarray(parent, np.int64), "t": t, "t0": t0,
            "bcount": bcount, "bslot0": bslot0, "bsym": bsym, "btime": btime,
            "gen_depth": depth, "route": route}


def assemble(out):
    """Tape state at every node from the per-branch record ALONE (bcount/bsym/btime), never the
    generator's hand-off state -- so a hand-off bug shows as a mismatch, not a silent copy.
    Returns (pathdepth (nn,k), tape_sym (nn,k,N), tape_time (nn,k,N), tape_origin (nn,k,N)),
    tape_origin = the branch that wrote each slot (-1 empty).  Two cells whose slot has the same
    origin share that edit -- which is what makes cell-level comparisons redundant."""
    parent, bcount = out["parent"], out["bcount"].astype(np.int64)
    nn, k, N = out["bsym"].shape
    root = nn - 1
    pdepth = np.zeros((nn, k), np.int64)
    tsym = np.full((nn, k, N), -1, np.int16)
    ttime = np.full((nn, k, N), np.nan, np.float32)
    torig = np.full((nn, k, N), -1, np.int32)
    zk = np.arange(k)
    for e in range(nn - 1, -1, -1):
        if e != root:
            p = parent[e]
            tsym[e], ttime[e], torig[e], base = tsym[p], ttime[p], torig[p], pdepth[p]
        else:
            base = np.zeros(k, np.int64)
        for j in range(int(bcount[e].max(initial=0))):
            m = bcount[e] > j
            slot = base[m] + j
            ok = slot < N                                # an overflow is a defect; D counts it
            zm, sm = zk[m][ok], slot[ok]
            tsym[e, zm, sm] = out["bsym"][e, zm, j]
            ttime[e, zm, sm] = out["btime"][e, zm, j]
            torig[e, zm, sm] = e
        pdepth[e] = base + bcount[e]
    return pdepth, tsym, ttime, torig


def save_decorated(path, out, rec: Recorder, seed_info: dict):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, parent=out["parent"], t=out["t"], bcount=out["bcount"],
                        bslot0=out["bslot0"], bsym=out["bsym"], btime=out["btime"],
                        meta=json.dumps({"recorder": rec.settings(), "route": out["route"],
                                         "n": out["n"], **seed_info}))


def load_decorated(path):
    z = np.load(path)
    out = {key: z[key] for key in ("parent", "t", "bcount", "bslot0", "bsym", "btime")}
    out["meta"] = json.loads(str(z["meta"]))
    out["n"] = out["meta"]["n"]
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser(description="decorate one library tree and summarise it")
    ap.add_argument("--n", type=int, default=210)
    ap.add_argument("--rho", type=float, default=0.25)
    ap.add_argument("--theta", type=float, default=0.5)
    ap.add_argument("--rep", type=int, default=0)
    ap.add_argument("--k", type=int, default=30)
    ap.add_argument("--seed", type=int, default=24092026)
    a = ap.parse_args()
    rec = Recorder(k=a.k)
    parent, branch = load_library_trees(a.n, a.rho, a.theta)[a.rep]
    out = decorate(parent, branch, a.n, rec, np.random.default_rng(a.seed))
    pd, _, _, _ = assemble(out)
    tip = pd[:a.n]
    print(f"n={a.n} k={a.k} M={rec.M} q_L={rec.q_chan[-1]:.5f}  "
          f"mean tip depth {tip.mean():.3f} (closed form "
          f"{np.sum(np.minimum(np.arange(60), rec.N) * _pois(np.arange(60), rec.Lam_T)):.3f})")


def _pois(x, mu):
    from scipy.stats import poisson
    return poisson.pmf(x, mu)


if __name__ == "__main__":
    main()
