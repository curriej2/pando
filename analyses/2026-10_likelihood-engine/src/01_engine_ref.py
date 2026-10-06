#!/usr/bin/env python3
r"""
01_engine_ref.py -- the NumPy REFERENCE TWIN of the likelihood engine.

Approved by Justin 2026-10-06 ("PROPOSAL (2026-10-06), PART 2", ../CLAUDE.md). Written to be READ:
plain loops, one tape at a time, every formula of the write-up (../writeup/likelihood_engine.tex)
spelled out. Slow and obviously correct; the JAX engine (02) is checked against it.

Also here, for check V1: an independent DENSE calculator that works on the FULL state space (every
string of length <= N, plus the silenced state) with no slot or lcp logic at all, and a literal
enumeration of every hidden history for the smallest trees.

CONVENTIONS (the tree library's): nodes 0..n-1 are the cells (tips); internal nodes have larger ids
than their children; the root is the last node; parent[root] = -1. Node times t: tips at 1, the
founder (origin) at 0, so the root's branch is the STEM [0, t_root]. Slots 0..N are depths; slot N+1
is the silenced state S.

MODEL (rung 2 of the simulator). Per tape z, along every lineage: editing at rate
r_z Lam_T lambda_0(t) m_branch, symbol i with probability xi_i(t); silencing at rate
Gam_T zeta_0(t). Readout: silenced -> missing; else lost with pi_cz = sigma(alpha_c + beta_z);
captured but blank -> missing; else read exactly.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import gammainc, gammaln
from scipy.stats import poisson


def sigmoid(y):
    return 1.0 / (1.0 + np.exp(-np.asarray(y, float)))


# ---------------------------------------------------------------------------------------------
# parameters: piecewise-constant shapes, composition epochs
# ---------------------------------------------------------------------------------------------

@dataclass
class PWShape:
    """A piecewise-constant shape on [0,1] with integral 1; W(t) = its cumulative (piecewise linear)."""
    knots: np.ndarray            # 0 = k_0 < ... < k_G = 1
    levels: np.ndarray           # G levels (normalised here)

    def __post_init__(self):
        self.knots = np.asarray(self.knots, float)
        lv = np.asarray(self.levels, float)
        self.levels = lv / np.sum(lv * np.diff(self.knots))
        self.Wk = np.concatenate([[0.0], np.cumsum(self.levels * np.diff(self.knots))])
        self.Wk[-1] = 1.0

    def W(self, t):
        return np.interp(t, self.knots, self.Wk)


@dataclass
class Params:
    """Everything the likelihood needs besides the tree and the reads."""
    N: int
    Lam_T: float                      # expected edits / tape / experiment for an average tape
    r: np.ndarray                     # (k,) tape speeds
    rate: PWShape                     # lambda_0
    comp_knots: np.ndarray            # composition epochs 0 = c_0 < ... < c_E = 1
    xi: np.ndarray                    # (E, M) symbol probabilities per composition epoch
    Gam_T: float                      # total silencing exposure per lineage (= -ln(1 - Hbar))
    sil: PWShape                      # zeta_0
    alpha: np.ndarray                 # (n,) cell offsets (logits)
    beta: np.ndarray                  # (k,) tape offsets (logits)
    b: float = 1.0                    # tree prior: division rate
    delta: float = 0.0                # death rate
    rho: float = 1.0                  # sampling fraction
    m: np.ndarray | None = None       # (n_nodes,) per-branch rate multipliers (hook; default 1)


# ---------------------------------------------------------------------------------------------
# branch amounts: expected edits mu and silencing exposure G on each composition epoch
# ---------------------------------------------------------------------------------------------

def branch_epochs(t0, t1, z, p: Params, mult=1.0):
    """List of (mu_e, G_e, xi_e) for the composition epochs overlapped by [t0, t1] (zero-length
    epochs are dropped: they would contribute the identity)."""
    out = []
    ck = np.asarray(p.comp_knots, float)
    for e in range(len(ck) - 1):
        a, b = max(t0, ck[e]), min(t1, ck[e + 1])
        if b <= a:
            continue
        mu = p.r[z] * p.Lam_T * mult * (p.rate.W(b) - p.rate.W(a))
        G = p.Gam_T * (p.sil.W(b) - p.sil.W(a))
        out.append((mu, G, p.xi[e]))
    return out


# ---------------------------------------------------------------------------------------------
# the one-segment kernel on depth slots (write-up eq. kernel) and its chaining over segments
# ---------------------------------------------------------------------------------------------

def count_prob(d0, d1, mu, N):
    """C_N(d0 -> d1; mu): probability a tape at depth d0 is at depth d1 after mu expected edits."""
    if d1 < d0:
        return 0.0
    if d1 < N:
        return float(poisson.pmf(d1 - d0, mu))
    if d0 == N:
        return 1.0
    return float(poisson.sf(N - d0 - 1, mu))         # P(at least N - d0 edits)


def segment_kernel(mu, G, xi, child_lcp, N):
    """(N+2) x (N+2) kernel for one segment. child_lcp = tuple of symbols (anchored child) or None
    (floating child). Slot N+1 = S."""
    K = np.zeros((N + 2, N + 2))
    surv = math.exp(-G)
    top = N if child_lcp is None else len(child_lcp)
    pmf = poisson.pmf(np.arange(N + 1), mu)          # P(exactly j edits), j = 0..N
    tail = poisson.sf(np.arange(N + 1) - 1, mu)       # tail[c] = P(at least c edits)
    for d0 in range(N + 1):
        for d1 in range(d0, top + 1):
            if child_lcp is None:
                sym_p = 1.0
            else:
                sym_p = float(np.prod([xi[s] for s in child_lcp[d0:d1]])) if d1 > d0 else 1.0
            if d1 < N:
                cnt = pmf[d1 - d0]
            else:
                cnt = 1.0 if d0 == N else tail[N - d0]
            K[d0, d1] = surv * cnt * sym_p
        if child_lcp is None:
            K[d0, N + 1] = 1.0 - surv
    K[N + 1, N + 1] = 1.0
    return K


def branch_kernel(t0, t1, z, child_lcp, p: Params, mult=1.0):
    """Product of the segment kernels in time order (identity if the branch has zero length)."""
    K = np.eye(p.N + 2)
    for mu, G, xi in branch_epochs(t0, t1, z, p, mult):
        K = K @ segment_kernel(mu, G, xi, child_lcp, p.N)
    return K


# ---------------------------------------------------------------------------------------------
# reads, anchors, lcps
# ---------------------------------------------------------------------------------------------

def reads_from_sim(D, sym, Y):
    """Simulator arrays -> reads[c][z] = observed tuple (depth >= 1) or None (missing)."""
    n, k = Y.shape
    return [[tuple(int(s) for s in sym[c, z, :D[c, z]]) if Y[c, z] else None for z in range(k)]
            for c in range(n)]


def _lcp(a, b):
    m = 0
    while m < min(len(a), len(b)) and a[m] == b[m]:
        m += 1
    return a[:m]


def anchors(parent, n, reads_z):
    """lcp[v] = longest common prefix of the reads below v (None = floating), for one tape."""
    nn = len(parent)
    lcp = [None] * nn
    for c in range(n):
        lcp[c] = reads_z[c]
    for v in range(nn - 1):                          # children before parents
        pv = parent[v]
        if lcp[v] is None:
            continue
        lcp[pv] = lcp[v] if lcp[pv] is None else _lcp(lcp[pv], lcp[v])
    return lcp


def children_of(parent):
    ch = [[] for _ in range(len(parent))]
    for v, pv in enumerate(parent):
        if pv >= 0:
            ch[pv].append(v)
    return ch


# ---------------------------------------------------------------------------------------------
# per-tape log-likelihood by pruning on depth slots
# ---------------------------------------------------------------------------------------------

def loglik_tape(parent, t, n, reads_z, z, p: Params):
    """ell_z = log P(reads on tape z | tree, params), pruning on depth slots with rescaling."""
    N, nn = p.N, len(parent)
    S = N + 1
    m = np.ones(nn) if p.m is None else p.m
    lcp = anchors(parent, n, reads_z)
    ch = children_of(parent)
    L = [None] * nn
    logscale = 0.0
    for c in range(n):                               # tips: emission on slots
        pi = sigmoid(p.alpha[c] + p.beta[z])
        vec = np.zeros(N + 2)
        if reads_z[c] is None:
            vec[:N + 1] = pi
            vec[0] = 1.0                             # pi + (1 - pi) at depth 0
            vec[S] = 1.0
        else:
            vec[len(reads_z[c])] = 1.0 - pi
        L[c] = vec
    for v in range(n, nn):                           # internal nodes, children first
        vec = np.ones(N + 2)
        for w in ch[v]:
            K = branch_kernel(t[v], t[w], z, lcp[w], p, m[w])
            vec = vec * (K @ L[w])
        if lcp[v] is not None:                       # anchored: only prefixes of the lcp; never S
            vec[len(lcp[v]) + 1:N + 1] = 0.0
            vec[S] = 0.0
        s = vec.max()
        if s > 0:
            vec = vec / s
            logscale += math.log(s)
        L[v] = vec
    root = nn - 1
    Kstem = branch_kernel(0.0, t[root], z, lcp[root], p, m[root])
    tot = float(Kstem[0] @ L[root])
    return logscale + math.log(tot) if tot > 0 else -math.inf


def loglik(parent, t, n, reads, p: Params):
    """Per-tape array (ell_1..ell_k), nats."""
    k = len(reads[0])
    return np.array([loglik_tape(parent, t, n, [reads[c][z] for c in range(n)], z, p)
                     for z in range(k)])


# ---------------------------------------------------------------------------------------------
# coordinates and the tree prior
# ---------------------------------------------------------------------------------------------

def times_from_x(x, parent, n):
    """Ratio transform: h_root = sigma(x_root), h_v = sigma(x_v) h_parent; t = 1 - h; tips at t = 1.
    x has one entry per internal node, in node-id order (x[0] belongs to node n)."""
    nn = len(parent)
    h = np.zeros(nn)
    for v in range(nn - 1, n - 1, -1):               # parents before children
        hp = 1.0 if parent[v] < 0 else h[parent[v]]
        h[v] = sigmoid(x[v - n]) * hp
    return 1.0 - h


def x_from_times(t, parent, n):
    nn = len(parent)
    h = 1.0 - np.asarray(t, float)
    x = np.zeros(nn - n)
    for v in range(n, nn):
        hp = 1.0 if parent[v] < 0 else h[parent[v]]
        f = h[v] / hp
        x[v - n] = math.log(f / (1.0 - f))
    return x


def log_jacobian(x, parent, n):
    t = times_from_x(x, parent, n)
    h = 1.0 - t
    tot = 0.0
    for v in range(n, len(parent)):
        hp = 1.0 if parent[v] < 0 else h[parent[v]]
        s = sigmoid(x[v - n])
        tot += math.log(hp) + math.log(s) + math.log(1.0 - s)
    return tot


def bds_hfun(h, b, delta, rho):
    """1 - p_0(h): probability a lineage alive at height h leaves >= 1 sampled descendant."""
    r = b - delta
    K = b * (1 - rho) - delta
    return rho * r / (b * rho + K * np.exp(-r * np.asarray(h, float)))


def bds_p1(h, b, delta, rho):
    """Probability a lineage alive at height h leaves exactly one sampled descendant."""
    r = b - delta
    K = b * (1 - rho) - delta
    e = np.exp(-r * np.asarray(h, float))
    return rho * r * r * e / (b * rho + K * e) ** 2


def log_tree_prior(t, parent, n, b, delta, rho):
    """Birth-death-sampling prior on the internal node heights, conditioned on n cells and on the
    founder at height 1: the n-1 branching heights are iid with density p_1(h) / int_0^1 p_1, on
    (0, 1). (int_0^1 p_1 = [H(1) - H(0)] / K with H = 1 - p_0, K = b(1-rho) - delta, H(0) = rho.)
    Up to a constant that depends only on n (the count of ranked topologies)."""
    h = 1.0 - np.asarray(t, float)[n:]
    r = b - delta
    K = b * (1 - rho) - delta
    # int_0^1 p_1 = [H(1) - rho] / K, and H(1) - rho = rho K (1 - e^{-r}) / (b rho + K e^{-r}),
    # so the K cancels: no division by K (which is 0 when delta = b(1 - rho)).
    norm = rho * (-math.expm1(-r)) / (b * rho + K * math.exp(-r))
    return float(np.sum(np.log(bds_p1(h, b, delta, rho))) - len(h) * math.log(norm))


# ---------------------------------------------------------------------------------------------
# V1 reference: the FULL state space, no slots, no lcps
# ---------------------------------------------------------------------------------------------

def full_states(N, M):
    st = [()] + [s for j in range(1, N + 1) for s in itertools.product(range(M), repeat=j)]
    return st + ['S']


def full_kernel(states, mu_G_xi, N):
    """Dense kernel over every state, from the model's definition, chained over segments."""
    nS = len(states)
    idx = {s: i for i, s in enumerate(states)}
    Ktot = np.eye(nS)
    for mu, G, xi in mu_G_xi:
        K = np.zeros((nS, nS))
        surv = math.exp(-G)
        for a in states:
            if a == 'S':
                K[idx[a], idx['S']] = 1.0
                continue
            K[idx[a], idx['S']] = 1.0 - surv
            for bb in states:
                if bb == 'S' or len(bb) < len(a) or bb[:len(a)] != a:
                    continue
                j = len(bb) - len(a)
                if len(bb) < N:
                    cnt = math.exp(-mu) * mu ** j / math.factorial(j) if mu > 0 or j == 0 else 0.0
                else:
                    cnt = 1.0 - sum(math.exp(-mu) * mu ** i / math.factorial(i)
                                    for i in range(N - len(a)))
                    if len(a) == N:
                        cnt = 1.0
                K[idx[a], idx[bb]] = surv * cnt * float(np.prod([xi[s] for s in bb[len(a):]]))
        Ktot = Ktot @ K
    return Ktot


def full_emission(states, read, pi):
    e = np.zeros(len(states))
    for i, s in enumerate(states):
        if read is None:
            e[i] = 1.0 if s == 'S' else (pi + (1 - pi) * (len(s) == 0))
        else:
            e[i] = 0.0 if s == 'S' else (1 - pi) * (s == read)
    return e


def loglik_tape_full(parent, t, n, reads_z, z, p: Params, M):
    """ell_z on the FULL state space by dense pruning (independent of slots/lcps)."""
    states = full_states(p.N, M)
    nn = len(parent)
    m = np.ones(nn) if p.m is None else p.m
    ch = children_of(parent)
    L = [None] * nn
    for c in range(n):
        L[c] = full_emission(states, reads_z[c], sigmoid(p.alpha[c] + p.beta[z]))
    for v in range(n, nn):
        vec = np.ones(len(states))
        for w in ch[v]:
            vec = vec * (full_kernel(states, branch_epochs(t[v], t[w], z, p, m[w]), p.N) @ L[w])
        L[v] = vec
    root = nn - 1
    K = full_kernel(states, branch_epochs(0.0, t[root], z, p, m[root]), p.N)
    return math.log(float(K[states.index(())] @ L[root]))


def loglik_tape_enumerate(parent, t, n, reads_z, z, p: Params, M):
    """The most literal check: list every assignment of states to the internal nodes and add up the
    probability of each complete story. Feasible only for 2-3 internal nodes."""
    states = full_states(p.N, M)
    nn = len(parent)
    idx = {s: i for i, s in enumerate(states)}
    m = np.ones(nn) if p.m is None else p.m
    Kb = {w: full_kernel(states, branch_epochs(t[parent[w]] if parent[w] >= 0 else 0.0, t[w], z, p,
                                               m[w]), p.N) for w in range(nn)}
    tipterm = {c: Kb[c] @ full_emission(states, reads_z[c], sigmoid(p.alpha[c] + p.beta[z]))
               for c in range(n)}
    internal = list(range(n, nn))
    tot = 0.0
    for assign in itertools.product(range(len(states)), repeat=len(internal)):
        x = dict(zip(internal, assign))
        pr = Kb[nn - 1][idx[()], x[nn - 1]]
        for v in internal:
            if v == nn - 1 or pr == 0.0:
                continue
            pr *= Kb[v][x[parent[v]], x[v]]
        if pr == 0.0:
            continue
        for c in range(n):
            pr *= tipterm[c][x[parent[c]]]
        tot += pr
    return math.log(tot)
