#!/usr/bin/env python3
r"""
17_dropout.py -- the dropout layer, RUNG 1: technical dropout with no lineage structure.

Approved by Justin 2026-09-30 ("PROPOSAL (2026-09-29)" in ../CLAUDE.md, as amended).
Validation: 18_validate_dropout.py.  Library: 19_dropout_library.py.  Nothing scientific is
concluded here.

THE MODEL.  Four parts, applied to tapes decorated by 14_editing.py.
  1. A TAPE PANEL, drawn once per simulated experiment and shared by every cell of every tree in
     it (all Park clones share the same 166 integrations).  Per tape z:
       u_z    = 1[U_z < phi]                       closed integration site (the latent class)
       beta_z = mu_beta + Delta_beta u_z + s_beta eps_z      logit of P(technically missing)
       r_z    = rbar_0 r_cl^{u_z} exp(s_r v_z - s_r^2/2),   rbar_0 = 1/(1 - phi + phi r_cl)
     U, eps, v are the panel's common random numbers: every (phi, r_cl, s_beta, pbar) setting is a
     deterministic function of the same draws, so the closed set is nested as phi rises and the
     worst-recovered tapes stay the worst as s_beta changes.  rbar_0 holds E[r] = 1 (the AVERAGE
     tape's speed), so Lam_T keeps its meaning.  mu_beta is SOLVED so that the panel's realised
     tapes give mean technical missingness pbar (averaged over alpha ~ N(0, s_alpha^2)).
  2. The TECHNICAL MASK.  X_cz = 1[V_cz < sigma(alpha_c + beta_z)], alpha_c = s_alpha a_c,
     a_c ~ N(0,1) iid over cells and blind to the tree -- the definition of rung 1.  1 = missing
     (the ladder's convention).  V, a are per-tree common random numbers.
  3. DEPTH-0 CENSORING.  Y_cz = (1 - X_cz) 1[D_cz >= 1]: a recovered but unedited tape reads as
     missing, exactly as 23_dropout_matrix.py codes a tape whose first site is None as ABSENT.
  4. CELL FILTERS, at export: QC keeps R_c = sum_z Y_cz >= R_min (applied downstream by the
     scripts themselves, via --thr); ClonalBC retention keeps a cell with probability psi(R_c)
     (clone = -1 otherwise, which is how 23/44/77 exclude a cell).

EXPORT.  export_arm() writes dropout_matrix_{arm}.npz (23's format) and prefix_codes6_{arm}.npz
(44's format, rows = cells with R_c >= R_min and clone >= 0, in dropout-matrix order), so
../../2026-08_park-compatibility/src/73, 74, 77 read simulated data through identical code.
Symbols are simulator ids, never Park sequences.  There is no JUNK and no internal gap in
simulated tapes, so 44's `seen` flag is exactly Y.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys
from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import brentq

_HERE = pathlib.Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("editing", _HERE / "14_editing.py")
ed = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("editing", ed)
_spec.loader.exec_module(ed)

RES = _HERE.parent / "results"
ABSENT, UNEDITED, JUNK, COMPLETE = 0, 1, 2, 3          # 23_dropout_matrix.py's term codes
_GH_X, _GH_W = np.polynomial.hermite_e.hermegauss(80)   # E_{N(0,1)}[f] = sum w f / sqrt(2 pi)
_GH_W = _GH_W / np.sqrt(2 * np.pi)


def sigmoid(x):
    x = np.asarray(x, float)
    with np.errstate(over="ignore"):
        return np.where(np.isneginf(x), 0.0, 1.0 / (1.0 + np.exp(-x)))


# --------------------------------------------------------------------------
# 1. the tape panel
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class PanelParams:
    phi: float = 0.0          # fraction of closed tapes
    r_cl: float = 1.0         # closed-tape speed relative to open (1 = no slowdown)
    s_r: float = 0.0          # ordinary tape-to-tape speed scatter (log scale)
    delta_beta: float = 0.0   # closed-tape read penalty, logits


@dataclass(frozen=True)
class MaskParams:
    pbar: float = 0.0         # target mean technical missingness (0 = no technical mask)
    s_beta: float = 0.0       # tape-to-tape recovery scatter within a class, logits
    s_alpha: float = 0.0      # cell-to-cell capture scatter, logits


def panel_draws(k, rng) -> dict:
    """The panel's common random numbers: U (class), eps (recovery), v (speed)."""
    return {"U": rng.random(k), "eps": rng.standard_normal(k), "v": rng.standard_normal(k)}


def panel_tapes(draws, pp: PanelParams):
    """(u, r): the class of every tape and its speed, E[r] = 1 by construction."""
    u = (draws["U"] < pp.phi).astype(np.int8)
    rbar0 = 1.0 / (1.0 - pp.phi + pp.phi * pp.r_cl)
    r = rbar0 * pp.r_cl ** u * np.exp(pp.s_r * draws["v"] - pp.s_r ** 2 / 2)
    return u, r


def beta_offset(draws, u, pp: PanelParams, mp: MaskParams):
    """beta_z - mu_beta = Delta_beta u_z + s_beta eps_z."""
    return pp.delta_beta * u + mp.s_beta * draws["eps"]


def mean_missing(mu_beta, base, s_alpha):
    """Mean over the panel's tapes of E_alpha sigma(mu_beta + base_z + alpha)."""
    x = mu_beta + base[:, None] + s_alpha * _GH_X[None, :]
    return float(np.mean(sigmoid(x) @ _GH_W))


def solve_mu_beta(base, mp: MaskParams):
    """mu_beta giving mean technical missingness pbar on this panel; -inf when pbar = 0."""
    if mp.pbar <= 0.0:
        return -np.inf
    assert mp.pbar < 1.0
    return brentq(lambda m: mean_missing(m, base, mp.s_alpha) - mp.pbar, -80.0, 80.0,
                  xtol=1e-12, rtol=1e-12)


def panel_beta(draws, u, pp, mp):
    base = beta_offset(draws, u, pp, mp)
    mu = solve_mu_beta(base, mp)
    return mu + base, mu


# --------------------------------------------------------------------------
# 2-3. the mask and the observed view
# --------------------------------------------------------------------------

def mask_draws(n, k, rng) -> dict:
    """Per-tree common random numbers: a (cell capture), V (entries), W (ClonalBC retention)."""
    return {"a": rng.standard_normal(n), "V": rng.random((n, k)), "W": rng.random(n)}


def technical_mask(md, beta, s_alpha):
    """X (n,k) bool, 1 = technically missing; pi = its probability sigma(alpha_c + beta_z)."""
    pi = sigmoid(s_alpha * md["a"][:, None] + np.asarray(beta)[None, :])
    return md["V"] < pi, pi


def observe(X, D):
    """Y = recovered AND edited at least once: 23's ABSENT rule (depth 0 reads as missing)."""
    return (~X) & (D >= 1)


def psi_of(R, psi):
    """ClonalBC retention probability at R_c: psi = ((lo_1, p_1), (lo_2, p_2), ...), lo_1 = 0,
    a step function; None = keep everyone."""
    if psi is None:
        return np.ones(np.shape(R))
    lo = np.array([s[0] for s in psi], float)
    p = np.array([s[1] for s in psi], float)
    assert lo[0] == 0 and np.all(np.diff(lo) > 0)
    return p[np.searchsorted(lo, R, side="right") - 1]


def retained(Y, md, psi):
    """ClonalBC retention per cell (QC is applied downstream via --thr)."""
    return md["W"] < psi_of(Y.sum(1), psi)


# --------------------------------------------------------------------------
# editing
# --------------------------------------------------------------------------

def decorate_tips(parent, branch, n, rec, rng):
    """Tip depth D (n,k) int8 and tip symbols (n,k,N) int16 (-1 empty), plus the full record."""
    out = ed.decorate(parent, branch, n, rec, rng)
    pd, tsym, _, _ = ed.assemble(out)
    return pd[:n].astype(np.int8), tsym[:n].astype(np.int16), out, pd


# --------------------------------------------------------------------------
# 4. export in the Park caches' formats
# --------------------------------------------------------------------------

def prefix_codes(S, A):
    """44_prefix_codes6.py's code construction, verbatim in logic: code at depth d identifies the
    symbol prefix S[..., :d+1]; -1 where any of those slots is empty."""
    n, K, MAX_D = S.shape
    codes = np.full((n, K, MAX_D), -1, dtype=np.int32)
    prev = np.zeros((n, K), dtype=np.int64)
    ok = np.ones((n, K), dtype=bool)
    for d in range(MAX_D):
        ok &= S[:, :, d] >= 0
        prev = prev * A + (S[:, :, d].astype(np.int64) + 1)
        flat = np.where(ok, prev, -1).ravel()
        _, inv = np.unique(flat, return_inverse=True)
        codes[:, :, d] = np.where(ok, inv.reshape(n, K).astype(np.int32), -1)
    return codes


def export_arm(arm, parts, R_min, M, outdir, N=6):
    """parts: list of dicts with D (n,k), sym (n,k,N), Y (n,k) bool, kept (n,) bool, name str --
    one per clone (tree).  Writes dropout_matrix_{arm}.npz and prefix_codes6_{arm}.npz."""
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    D = np.concatenate([p["D"] for p in parts]).astype(np.int8)
    Y = np.concatenate([p["Y"] for p in parts]).astype(bool)
    S = np.concatenate([p["sym"] for p in parts]).astype(np.int16)
    kept = np.concatenate([p["kept"] for p in parts]).astype(bool)
    cid = np.concatenate([np.full(p["D"].shape[0], i, np.int32) for i, p in enumerate(parts)])
    n, K = D.shape
    depth = np.where(Y, D, 0).astype(np.int8)
    term = np.where(~Y, ABSENT, np.where(D >= N, COMPLETE, UNEDITED)).astype(np.int8)
    clone = np.where(kept, cid, -1).astype(np.int32)
    tapes = np.array([f"tape{z:03d}" for z in range(K)])
    np.savez_compressed(
        outdir / f"dropout_matrix_{arm}.npz", tapes=tapes, depth=depth, term=term, recovered=Y,
        sample=np.zeros(n, np.int16), sample_names=np.array(["sim"]),
        clone=clone, clone_names=np.array([p["name"] for p in parts]))
    rows = (Y.sum(1) >= R_min) & kept                    # 44: ClonalBC present and ntapes >= THR
    slot = np.arange(N)[None, None, :]
    Sv = np.where(Y[..., None] & (slot < D[..., None]), S, -1)[rows]
    codes = prefix_codes(Sv, M + 1)
    np.savez_compressed(outdir / f"prefix_codes6_{arm}.npz", codes=codes)
    return {"cells": int(n), "tapes": int(K), "clones": len(parts), "rows_prefix": int(rows.sum()),
            "R_min": int(R_min)}


def settings_json(**kw):
    def conv(x):
        if isinstance(x, (PanelParams, MaskParams)):
            return asdict(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        return x
    return json.dumps({k: conv(v) for k, v in kw.items()}, default=float)
