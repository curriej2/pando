#!/usr/bin/env python3
r"""
02_engine.py -- the likelihood engine, in JAX: exact per-tape log-likelihood of one clone, with
gradients. Same mathematics as the reference twin 01_engine_ref.py (and the write-up,
../writeup/likelihood_engine.tex), written for speed and automatic differentiation.

Approved by Justin 2026-10-06 ("PROPOSAL (2026-10-06), PART 2", ../CLAUDE.md). Runs in the separate
environment /data1/choij10/justin/envs/pando-jax. Double precision throughout.

DESIGN (notes/sciphy_notes.md H.8).
  * Depth slots: every (node, tape) carries N+2 numbers (depths 0..N, then S = silenced).
  * Topology as DATA: parent/children index arrays and per-(node, tape) anchors are arrays whose
    shapes depend only on (n, k, N), so a new topology never triggers recompilation.
  * Internal nodes are processed in id order (children have smaller ids than parents) by lax.scan,
    all tapes at once; each step rescales its vector and accumulates the log scale.
  * Branch kernel = product over composition epochs of the closed-form one-segment kernel
    (zero-length epochs give the identity); expected edits mu and silencing exposure G come from the
    piecewise-linear cumulatives W(t), W_zeta(t) by exact interpolation.

PUBLIC FUNCTIONS
  prepare(parent, n, D, sym, Y, N)      -> Data   (numpy; topology- and reads-dependent)
  make_loglik(N, rate_knots, comp_knots, sil_knots) -> loglik(t, prm, data) -> (k,) per-tape nats
  times_from_x(x, data), log_jacobian(x, data), log_tree_prior(t, b, delta, rho, data)
  make_objective(...)                    -> J(x, prm, data, sample=False)
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import jax
import jax.numpy as jnp
from jax.scipy.special import gammainc, gammaln, xlogy

jax.config.update("jax_enable_x64", True)


class Data(NamedTuple):
    """Everything topology- and read-dependent, as arrays (a JAX pytree)."""
    parent: jnp.ndarray       # (nn,)  parent id, -1 for the root
    children: jnp.ndarray     # (nn, 2) child ids for internal nodes (tips: -1)
    internal: jnp.ndarray     # (n-1,) internal node ids in increasing order (a post-order)
    lcp_sym: jnp.ndarray      # (nn, k, N) lcp symbols below each node (-1 padded)
    lcp_len: jnp.ndarray      # (nn, k) lcp length; -1 = floating (nothing read below)
    tip_obs: jnp.ndarray      # (n, k) bool, tape read
    tip_len: jnp.ndarray      # (n, k) read length (0 if missing)


def prepare(parent, n, D, sym, Y, N):
    """Build the Data arrays from a library tree and the simulator's reads (D depth, sym symbols,
    Y observed). Longest common prefixes are combined child -> parent, vectorised over tapes."""
    parent = np.asarray(parent, np.int64)
    nn = parent.size
    k = Y.shape[1]
    seq = -np.ones((nn, k, N), np.int64)
    ln = -np.ones((nn, k), np.int64)
    slot = np.arange(N)[None, :]
    seq[:n] = np.where(Y[..., None] & (slot[None] < D[..., None]), sym, -1)
    ln[:n] = np.where(Y, D, -1)
    for v in range(nn - 1):                                     # children before parents
        pv = parent[v]
        has = ln[v] >= 0
        if not has.any():
            continue
        fresh = has & (ln[pv] < 0)
        seq[pv, fresh], ln[pv, fresh] = seq[v, fresh], ln[v, fresh]
        both = has & ~fresh
        if both.any():
            m = np.minimum(ln[pv, both], ln[v, both])
            eq = (seq[pv, both] == seq[v, both]) & (slot < m[:, None])
            newlen = np.cumprod(eq, axis=1).sum(axis=1)
            ln[pv, both] = newlen
            s = seq[pv, both]
            s[slot >= newlen[:, None]] = -1
            seq[pv, both] = s
    children = -np.ones((nn, 2), np.int64)
    cnt = np.zeros(nn, np.int64)
    for v in range(nn - 1):
        children[parent[v], cnt[parent[v]]] = v
        cnt[parent[v]] += 1
    assert np.all(cnt[n:] == 2), "every internal node must have exactly two children"
    return Data(parent=jnp.asarray(parent), children=jnp.asarray(children),
                internal=jnp.arange(n, nn), lcp_sym=jnp.asarray(seq), lcp_len=jnp.asarray(ln),
                tip_obs=jnp.asarray(Y), tip_len=jnp.asarray(np.where(Y, D, 0)))


def _cumulative(knots, levels):
    """Knot values of W(t) = int_0^t shape, the shape normalised to integrate to 1."""
    dk = jnp.diff(knots)
    om = levels / jnp.sum(levels * dk)
    return jnp.concatenate([jnp.zeros(1), jnp.cumsum(om * dk)])


def _segment_kernel(mu, G, cl, ln, N):
    """One-segment kernel on depth slots for one tape. mu: expected edits; G: silencing exposure;
    cl: (N+1,) cumulative log symbol probabilities along the child's lcp (cl[0] = 0);
    ln: child's lcp length (-1 = floating). Returns (N+2, N+2); slot N+1 = S."""
    d0 = jnp.arange(N + 1)[:, None]
    d1 = jnp.arange(N + 1)[None, :]
    j = d1 - d0
    jj = jnp.maximum(j, 0).astype(mu.dtype)                     # float: xlogy differentiates both args
    pos = mu > 0
    mu_s = jnp.where(pos, mu, 1.0)                              # safe value: no 0 * inf gradients
    pmf = jnp.where(pos, jnp.exp(-mu_s + xlogy(jj, mu_s) - gammaln(jj + 1.0)), (jj == 0) * 1.0)
    a = jnp.maximum(N - d0, 1).astype(mu.dtype)                 # edits needed to fill
    tail = jnp.where(d0 < N, jnp.where(pos, gammainc(a, mu_s), 0.0), 1.0)   # P(>= N - d0 edits)
    C = jnp.where(j < 0, 0.0, jnp.where(d1 < N, pmf, tail))
    anch = ln >= 0
    sym = jnp.where((d1 <= ln) & (j >= 0), jnp.exp(cl[jnp.minimum(d1, N)] - cl[d0]), 0.0)
    Xi = jnp.where(anch, sym, 1.0)
    surv = jnp.exp(-G)
    K = jnp.zeros((N + 2, N + 2), mu.dtype)
    K = K.at[:N + 1, :N + 1].set(surv * C * Xi)
    K = K.at[:N + 1, N + 1].set(jnp.where(anch, 0.0, 1.0 - surv))
    return K.at[N + 1, N + 1].set(1.0)


def make_loglik(N, rate_knots, comp_knots, sil_knots):
    """Build the jitted per-tape log-likelihood for fixed knot positions (fixed model structure).
    prm keys: Lam_T, r (k,), rate_levels (G,), logxi (E, M), Gam_T, sil_levels (Gs,),
              alpha (n,), beta (k,), m (nn,)."""
    rate_knots = jnp.asarray(rate_knots, float)
    comp_knots = jnp.asarray(comp_knots, float)
    sil_knots = jnp.asarray(sil_knots, float)
    E = comp_knots.size - 1
    seg = jax.vmap(_segment_kernel, in_axes=(0, None, 0, 0, None))   # over tapes

    def branch_kernels(t0, t1, child, prm, data, Wk, Wz, cl_all):
        """(k, N+2, N+2) kernel of branch `child` over [t0, t1], all tapes."""
        K = None
        for e in range(E):
            a = jnp.clip(t0, comp_knots[e], comp_knots[e + 1])
            b = jnp.clip(t1, comp_knots[e], comp_knots[e + 1])
            dW = jnp.interp(b, rate_knots, Wk) - jnp.interp(a, rate_knots, Wk)
            mu = prm["r"] * prm["Lam_T"] * prm["m"][child] * dW            # (k,)
            G = prm["Gam_T"] * (jnp.interp(b, sil_knots, Wz) - jnp.interp(a, sil_knots, Wz))
            Ke = seg(mu, G, cl_all[e][child], data.lcp_len[child], N)
            K = Ke if K is None else jnp.einsum("zij,zjk->zik", K, Ke)
        return K

    @jax.jit
    def loglik(t, prm, data: Data):
        nn, k = data.lcp_len.shape
        n = data.tip_obs.shape[0]
        S = N + 1
        Wk = _cumulative(rate_knots, prm["rate_levels"])
        Wz = _cumulative(sil_knots, prm["sil_levels"])
        # cumulative log symbol probabilities along every node's lcp, per composition epoch
        sym = jnp.maximum(data.lcp_sym, 0)
        valid = data.lcp_sym >= 0
        cl_all = []
        for e in range(E):
            lx = jnp.where(valid, prm["logxi"][e][sym], 0.0)               # (nn, k, N)
            cl_all.append(jnp.concatenate([jnp.zeros((nn, k, 1)), jnp.cumsum(lx, axis=2)], axis=2))
        # tips
        pi = jax.nn.sigmoid(prm["alpha"][:, None] + prm["beta"][None, :])  # (n, k)
        dslot = jnp.arange(N + 2)[None, None, :]
        miss = jnp.where(dslot == 0, 1.0, jnp.where(dslot == S, 1.0, pi[..., None]))
        obs = jnp.where(dslot == data.tip_len[..., None], (1.0 - pi)[..., None], 0.0)
        tipvec = jnp.where(data.tip_obs[..., None], obs, miss)
        L0 = jnp.zeros((nn, k, N + 2)).at[:n].set(tipvec)
        ls0 = jnp.zeros((nn, k))

        def step(carry, v):
            L, ls = carry
            vec = jnp.ones((k, N + 2))
            lsv = jnp.zeros(k)
            for slot in range(2):
                c = data.children[v, slot]
                K = branch_kernels(t[v], t[c], c, prm, data, Wk, Wz, cl_all)
                vec = vec * jnp.einsum("zij,zj->zi", K, L[c])
                lsv = lsv + ls[c]
            lnv = data.lcp_len[v]
            keep = jnp.where((lnv >= 0)[:, None],
                             (dslot[0] <= lnv[:, None]) & (dslot[0] < S), True)
            vec = jnp.where(keep, vec, 0.0)
            s = jnp.max(vec, axis=1)
            s = jnp.where(s > 0, s, 1.0)
            return (L.at[v].set(vec / s[:, None]), ls.at[v].set(lsv + jnp.log(s))), None

        (L, ls), _ = jax.lax.scan(step, (L0, ls0), data.internal)
        root = nn - 1
        Kst = branch_kernels(0.0, t[root], root, prm, data, Wk, Wz, cl_all)
        tot = jnp.einsum("zj,zj->z", Kst[:, 0, :], L[root])
        return ls[root] + jnp.log(tot)

    return loglik


# -------------------------------------------------------------------------------------------------
# coordinates, Jacobian, tree prior, objective
# -------------------------------------------------------------------------------------------------

def times_from_x(x, data: Data):
    """Ratio transform. x: (n-1,) one per internal node in id order. Returns t (nn,), tips at 1."""
    nn = data.parent.shape[0]
    n = nn - x.shape[0]

    def step(h, v):
        pv = data.parent[v]
        hp = jnp.where(pv < 0, 1.0, h[jnp.maximum(pv, 0)])
        return h.at[v].set(jax.nn.sigmoid(x[v - n]) * hp), None

    h, _ = jax.lax.scan(step, jnp.zeros(nn), data.internal[::-1])   # parents before children
    return 1.0 - h


def log_jacobian(x, data: Data):
    t = times_from_x(x, data)
    h = 1.0 - t
    nn = data.parent.shape[0]
    n = nn - x.shape[0]
    pv = data.parent[n:]
    hp = jnp.where(pv < 0, 1.0, h[jnp.maximum(pv, 0)])
    return jnp.sum(jnp.log(hp) + jax.nn.log_sigmoid(x) + jax.nn.log_sigmoid(-x))


def log_tree_prior(t, b, delta, rho, data: Data):
    """Conditioned on n cells and the founder at height 1, the n-1 branching heights are iid with
    density p_1(h) / int_0^1 p_1 on (0, 1). Up to a constant depending only on n."""
    n = data.tip_obs.shape[0]
    h = 1.0 - t[n:]
    r = b - delta
    K = b * (1.0 - rho) - delta
    e = jnp.exp(-r * h)
    logp1 = jnp.log(rho) + 2.0 * jnp.log(r) - r * h - 2.0 * jnp.log(b * rho + K * e)
    norm = rho * (-jnp.expm1(-r)) / (b * rho + K * jnp.exp(-r))   # int_0^1 p_1, no 1/K
    return jnp.sum(logp1) - h.size * jnp.log(norm)


def make_objective(loglik):
    """J(x, prm, data, tree=(b, delta, rho), sample) = sum_z ell_z + log tree prior
    [+ log Jacobian when sampling]. Parameter priors are added by the caller."""
    def J(x, prm, data, tree, sample=False):
        t = times_from_x(x, data)
        val = jnp.sum(loglik(t, prm, data)) + log_tree_prior(t, *tree, data)
        return val + (log_jacobian(x, data) if sample else 0.0)
    return J


# -------------------------------------------------------------------------------------------------
# helpers to build engine inputs from the simulator's objects
# -------------------------------------------------------------------------------------------------

def comp_epochs_from_recorder(rec):
    """(comp_knots, xi (E, M)) from a 14_editing.Recorder: full symbol vector per share epoch."""
    E = rec.sh.shape[0]
    xi = np.zeros((E, rec.M))
    for e in range(E):
        p = rec.sh[e]
        for a in range(rec.A):
            xi[e, rec.offsets[a]:rec.offsets[a] + rec.chan_xi[a].size] = p[a] * rec.chan_xi[a]
        xi[e, rec.offsets[rec.A]:] = (1.0 - p.sum()) * rec.chan_xi[rec.A]
    return np.asarray(rec.ck, float), xi


def safe_logxi(xi):
    return np.log(np.maximum(xi, 1e-300))
