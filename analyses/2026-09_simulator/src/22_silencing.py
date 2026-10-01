#!/usr/bin/env python3
r"""
22_silencing.py -- dropout RUNG 2: heritable silencing as complete Dollo loss.

Approved by Justin 2026-10-01 ("PROPOSAL (2026-10-01, rev. 1)" in ../CLAUDE.md).
Validation: 23_validate_silencing.py.  Nothing scientific is concluded here.

THE MODEL.  Silencing strikes each tape along each lineage as a Poisson process in calendar time,
    zeta(t) = Gam_T zeta_0(t),   int_0^1 zeta_0 = 1,   Gam(t) = Gam_T W_zeta(t),
the same hazard on every tape (no tape heterogeneity) and no silencing before the clone founder
(one tree = one clone; tapes intact at t = 0).  Once silenced a tape is lost in every descendant
and never returns -- complete Dollo loss.  zeta_0 is piecewise constant on its OWN knots.

THE KNOB.  Hbar = the expected fraction of tip x tape entries silenced at harvest (t = 1), the TRUE
state, before the technical mask, depth-0 censoring and the cell filters.  Every tip's path carries
exposure Gam_T, so Hbar = 1 - exp(-Gam_T) and Gam_T = -ln(1 - Hbar), whatever the tree and
whatever zeta_0.  Solved on the EXPECTATION; the realised fraction is reported, never forced.

THE SAMPLER (route (ii) with one slot).  Per (branch e, tape z) one E ~ Exp(1), drawn from a stream
of its own:  L_ez = 1[E_ez <= Gam_e],  Gam_e = Gam_T [W_zeta(t1_e) - W_zeta(t0_e)],  and the first
event's time t* = W_zeta^{-1}(W_zeta(t0_e) + E_ez / Gam_T).  For fixed draws L is monotone in Hbar,
so silenced sets are NESTED across Hbar at a fixed zeta_0 (not across shapes).

INHERITANCE.  H_e = H_pa(e) OR L_e (node state at the lower end of branch e); nothing above the
stem.  A VISIBLE event tops a lineage's first loss: V_e = L_e AND NOT H_pa(e).

OBSERVATION.  Y = (1 - X)(1 - H) 1[D >= 1] -- observe().  Editing (D) is NOT altered by silencing:
under complete loss with no trans limb, post-silencing edits are never observed.  17's export_arm
takes Y unchanged, so rung-2 data are written in exactly rung 1's cache formats.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class SilencingParams:
    Hbar: float = 0.0                 # target expected silenced fraction of tip x tape entries
    knots: tuple = (0.0, 1.0)         # zeta_0's knots, 0 .. 1
    levels: tuple = (1.0,)            # zeta_0's levels (normalised so int zeta_0 = 1)

    @property
    def Gam_T(self) -> float:
        assert 0.0 <= self.Hbar < 1.0
        return float(-np.log1p(-self.Hbar))


class Shape:
    """zeta_0 piecewise constant; W_zeta its cumulative, with the exact two-step inverse of
    notes §S4.9.2 (a zero-hazard window occupies an empty interval of W and is never selected;
    a boundary tie resolves into the NEXT active epoch)."""

    def __init__(self, knots, levels):
        kap = np.asarray(knots, float)
        om = np.asarray(levels, float)
        assert kap[0] == 0.0 and kap[-1] == 1.0 and np.all(np.diff(kap) > 0)
        assert om.size == kap.size - 1 and np.all(om >= 0) and om.sum() > 0
        om = om / np.sum(om * np.diff(kap))
        self.kappa, self.omega = kap, om
        self.W = np.concatenate([[0.0], np.cumsum(om * np.diff(kap))])
        self.W[-1] = 1.0
        self._last_pos = int(np.nonzero(om > 0)[0][-1])

    def W_of(self, t):
        t = np.asarray(t, float)
        g = np.clip(np.searchsorted(self.kappa, t, side="right"), 1, self.omega.size)
        return self.W[g - 1] + self.omega[g - 1] * (t - self.kappa[g - 1])

    def W_inv(self, w):
        w = np.asarray(w, float)
        g = np.clip(np.searchsorted(self.W, w, side="right"), 1, self.omega.size)
        g = np.where(self.omega[g - 1] > 0, g, self._last_pos + 1)
        return self.kappa[g - 1] + (w - self.W[g - 1]) / self.omega[g - 1]

    def zero_windows(self):
        """Open intervals of t where zeta_0 = 0."""
        return [(self.kappa[g], self.kappa[g + 1]) for g in range(self.omega.size)
                if self.omega[g] == 0]


def loss_draws(nn, k, rng) -> np.ndarray:
    """The loss layer's common random numbers: E (nn, k) ~ Exp(1), one per (branch, tape)."""
    return rng.exponential(1.0, size=(nn, k))


def branch_ends(parent, t):
    """(t0, t1) of every branch e (the edge INTO node e); the root's edge is the stem [0, t_root]."""
    parent = np.asarray(parent, np.int64)
    root = parent.size - 1
    t0 = np.where(np.arange(parent.size) == root, 0.0, t[np.maximum(parent, 0)])
    return t0, np.asarray(t, float)


def silence(parent, t, E, sp: SilencingParams):
    """The loss layer on one tree.  parent: library convention (children have smaller ids, root
    = last, parent[root] = -1); t: node times with tips at 1 (14's node_times); E: loss_draws().
    Returns L (fired on branch), tstar (first event time, NaN where L = 0), H (node state at the
    lower end of each branch), V (visible event tops), Gam_e (branch hazards)."""
    parent = np.asarray(parent, np.int64)
    nn, k = E.shape
    assert parent.size == nn
    shape = Shape(sp.knots, sp.levels)
    t0, t1 = branch_ends(parent, t)
    W0, W1 = shape.W_of(t0), shape.W_of(t1)
    G = sp.Gam_T
    Gam_e = G * (W1 - W0)                                       # (nn,)
    L = (E <= Gam_e[:, None]) & (Gam_e[:, None] > 0)
    tstar = np.full((nn, k), np.nan)
    if L.any():
        ee, zz = np.nonzero(L)
        ts = shape.W_inv(W0[ee] + E[ee, zz] / G)
        tstar[ee, zz] = ts
    H = np.zeros((nn, k), bool)
    V = np.zeros((nn, k), bool)
    root = nn - 1
    for e in range(nn - 1, -1, -1):                              # parents before children
        above = np.zeros(k, bool) if e == root else H[parent[e]]
        H[e] = above | L[e]
        V[e] = L[e] & ~above
    return {"L": L, "tstar": tstar, "H": H, "V": V, "Gam_e": Gam_e, "t0": t0, "t1": t1,
            "Gam_T": G}


def clade_sizes(parent, n):
    """n_e = number of tips below node e."""
    parent = np.asarray(parent, np.int64)
    size = np.zeros(parent.size, np.int64)
    size[:n] = 1
    for x in range(parent.size - 1):                             # children before parents
        size[parent[x]] += size[x]
    return size


def p_visible(parent, t, sp: SilencingParams):
    """P(V_e = 1) = (1 - Hbar)^{W_zeta(t0_e)} [1 - (1 - Hbar)^{W_zeta(t1_e) - W_zeta(t0_e)}] (G3)."""
    shape = Shape(sp.knots, sp.levels)
    t0, t1 = branch_ends(parent, t)
    W0, W1 = shape.W_of(t0), shape.W_of(t1)
    G = sp.Gam_T
    return np.exp(-G * W0) * (-np.expm1(-G * (W1 - W0)))


def heritable_share(parent, t, n, sp: SilencingParams):
    """Expected fraction of silenced entries sitting in a visible event of >= 2 cells."""
    if sp.Hbar == 0:
        return float("nan")
    pv = p_visible(parent, t, sp)
    sz = clade_sizes(parent, n)
    return float(np.sum((sz * pv)[sz >= 2]) / (n * sp.Hbar))


def observe(X, D, H):
    """Y = recovered AND not silenced AND edited at least once (three independent factors)."""
    return (~X) & (~H) & (D >= 1)


def settings_json(sp: SilencingParams) -> str:
    return json.dumps({**asdict(sp), "Gam_T": sp.Gam_T})
