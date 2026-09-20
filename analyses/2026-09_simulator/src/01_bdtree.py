#!/usr/bin/env python3
r"""
01_bdtree.py -- forward birth-death-sampling tree simulator.  TREE LAYER ONLY.

Simulates, literally, the process defined in notes/sciphy_notes.md S3.3:

    start with one lineage at t = 0; every extant lineage independently divides
    at rate b and dies at rate delta; run to T; sample each survivor
    independently with probability rho; prune to the RECONSTRUCTED tree on the
    sampled tips (delete non-ancestors, then suppress degree-2 nodes).

Editing and dropout are NOT here -- separate layers, separate scripts.

UNITS.  b, delta are per-lineage rates, 1/time.  T is a time.  rho is a
dimensionless capture fraction in (0, 1].  n is a tip count in cells.

THE GILLESPIE LOOP.  With n lineages alive the total event propensity is
A = n(b + delta); the wait to the next event anywhere is Exp(A); and given that
an event happened it is a birth with probability b/(b+delta), independently of
when.  Those three facts ARE the simulator.  Exact -- no time discretisation,
no step-size parameter, no wasted "nothing happened" iterations.

TWO IMPLEMENTATIONS, deliberately.
  simulate_tree_naive  one event at a time, keep every node, prune at the end.
                       Slow and obviously correct.  The oracle that
                       02_validate_tree.py checks the fast path against.
  simulate_tree        production path.  Two passes over the SAME count stream:
                       pass 1 gets N(T) with no genealogy at all (vectorised in
                       numpy), so an attempt that will miss the target tip count
                       is rejected before a single node is allocated; pass 2
                       replays the identical trajectory from the same seed and
                       builds the genealogy ONLY on acceptance.

WHY THE TWO-PASS SPLIT IS EXACT.  At each event the lineage it lands on is
chosen uniformly and independently of everything else, so the law factorises as
(count trajectory) x (uniform lineage assignments | count trajectory).  The
acceptance test reads only N(T) and a Binomial(N(T), rho) draw -- never the
genealogy.  Conditioning on the count and building the genealogy afterwards
therefore samples the same distribution, and the ~2.7n rejected attempts cost
only the cheap vectorised pass.  Given Binomial(N, rho) = k, the captured set is
a uniform k-subset of the survivors, which is how pass 2 draws it.

PARAMETERISATION FOR COST WORK (rate_params).  alpha and beta depend on
(b, delta, T) only through u = e^{(b-delta)T} and the turnover theta = delta/b,
so the count process -- hence the whole rejection cost -- has exactly two free
numbers at fixed rho.  T is a pure scale choice and is fixed at 1.  (For
ADJUDICATION against Park that is no longer true: there T is pinned by
Lambda = lambda*T at the measured Lambda-hat.  It is true for the cost study.)
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------
# Closed forms from S3.3.2, used only to CHECK the simulator (02), never by it.
# --------------------------------------------------------------------------


def bd_alpha_beta(b: float, delta: float, T: float) -> tuple[float, float]:
    r"""alpha = P(extinct by T); beta from P(N=n) = (1-alpha)(1-beta)beta^{n-1}."""
    r = b - delta
    if abs(r) * T < 1e-12:                      # critical case, by continuity
        return delta * T / (1 + delta * T), b * T / (1 + b * T)
    u = np.exp(r * T)
    den = b * u - delta
    return delta * (u - 1) / den, b * (u - 1) / den


def bd_pmf(n, b: float, delta: float, T: float):
    r"""P(N(T) = n) for n >= 1.  Geometric tail -- exponentially bounded (S3.3.2)."""
    a, bet = bd_alpha_beta(b, delta, T)
    return (1 - a) * (1 - bet) * bet ** (np.asarray(n) - 1)


def bd_mean(b: float, delta: float, T: float) -> float:
    r"""E[N(T)] = e^{(b-delta)T}.  Sees only the difference (S3.3.1)."""
    return float(np.exp((b - delta) * T))


def rate_params(n_target: int, rho: float, turnover: float, T: float = 1.0):
    r"""(b, delta, T) tuned so E[N(T)] = n_target/rho at the requested turnover.

    Without this the process is not centred on the target and acceptance is
    tiny for a reason that has nothing to do with the question being asked.
    """
    if not 0.0 <= turnover < 1.0:
        raise ValueError("turnover = delta/b must lie in [0, 1)")
    r = np.log(n_target / rho) / T
    b = r / (1.0 - turnover)
    return float(b), float(turnover * b), float(T)


# --------------------------------------------------------------------------
# NOT USED BY THE SIMULATOR.  The BDS branching-time machinery (S3.3.5), kept
# here so the deferred forward-vs-density equivalence check needs no new code.
# --------------------------------------------------------------------------


def bds_h(t, b: float, delta: float, rho: float):
    r"""h(t) = 1 - p_0(t) = P(a lineage t before the present leaves a sampled descendant)."""
    r = b - delta
    K = b * (1 - rho) - delta
    return rho * r / (b * rho + K * np.exp(-r * np.asarray(t)))


def bds_G(x, b: float, delta: float, rho: float, t_or: float):
    r"""CDF of one branching time, measured before the present.  h(0) = rho."""
    h0 = rho
    return (bds_h(x, b, delta, rho) - h0) / (bds_h(t_or, b, delta, rho) - h0)


# --------------------------------------------------------------------------
# Pass 1 -- the count trajectory, vectorised.  No genealogy.
# --------------------------------------------------------------------------


@dataclass
class Trajectory:
    n_T: int                      # lineages alive at T; 0 = extinct before T
    n_events: int                 # births + deaths that occurred by T
    steps: np.ndarray | None = None   # +1 birth / -1 death, in order
    times: np.ndarray | None = None   # event times, same length


def simulate_counts(b, delta, T, rng, record=False,
                    chunk0=256, chunk_max=1 << 20, max_events=1 << 33) -> Trajectory:
    r"""Gillespie loop over the POPULATION SIZE only, in geometrically growing chunks.

    Draws are made in the identical order whether or not `record` is set, so a
    replay from the same seed reproduces the trajectory exactly -- which is what
    licenses the two-pass design.
    """
    total = b + delta
    p_birth = b / total
    n, t, done, chunk = 1, 0.0, 0, chunk0
    rs: list[np.ndarray] = []
    rt: list[np.ndarray] = []

    while True:
        if done > max_events:
            raise RuntimeError(f"max_events exceeded: {done} events without reaching T")

        steps = np.where(rng.random(chunk) < p_birth, np.int8(1), np.int8(-1))
        counts = n + np.cumsum(steps, dtype=np.int64)       # size AFTER each event

        z = np.flatnonzero(counts == 0)                     # extinction inside this chunk
        if z.size:
            k = int(z[0]) + 1
            steps, counts = steps[:k], counts[:k]

        pre = counts - steps                                # size BEFORE each event, >= 1
        times = t + np.cumsum(rng.exponential(1.0, steps.size) / (pre * total))

        over = np.flatnonzero(times > T)
        if over.size:                                       # T reached inside this chunk
            k = int(over[0])
            n_T = int(counts[k - 1]) if k else n
            if record:
                rs.append(steps[:k]); rt.append(times[:k])
            return _pack(n_T, done + k, record, rs, rt)

        if record:
            rs.append(steps); rt.append(times)
        done += steps.size

        if counts[-1] == 0:                                 # died out before T
            return _pack(0, done, record, rs, rt)

        n, t = int(counts[-1]), float(times[-1])
        chunk = min(chunk * 2, chunk_max)


def _pack(n_T, n_events, record, rs, rt) -> Trajectory:
    if not record:
        return Trajectory(n_T, n_events)
    s = np.concatenate(rs) if rs else np.empty(0, np.int8)
    u = np.concatenate(rt) if rt else np.empty(0, np.float64)
    return Trajectory(n_T, n_events, s, u)


# --------------------------------------------------------------------------
# Pass 2 -- the genealogy, replayed from the recorded trajectory.
# --------------------------------------------------------------------------


def build_genealogy(steps, times, rng_lin):
    r"""Assign each recorded event to a uniformly chosen live lineage.

    Node ids are allocated in creation order, so a child ALWAYS has a larger id
    than its parent.  `reconstruct` relies on that to order its accumulation.
    Returns (parent, tbirth, live) where tbirth[x] is when lineage x came into
    existence -- so the DIVISION time of a node is tbirth of either child.
    """
    n_births = int(np.count_nonzero(steps == 1))
    n_nodes = 1 + 2 * n_births
    parent = np.full(n_nodes, -1, np.int64)
    tbirth = np.zeros(n_nodes, np.float64)

    max_live = int((1 + np.cumsum(steps, dtype=np.int64)).max()) if steps.size else 1
    live = np.empty(max(max_live, 1), np.int64)
    live[0] = 0
    nl, nxt = 1, 1

    u = rng_lin.random(steps.size)               # one array draw, not one call per event
    for k in range(steps.size):
        i = int(u[k] * nl)
        x = live[i]
        if steps[k] == 1:
            parent[nxt] = x; tbirth[nxt] = times[k]
            parent[nxt + 1] = x; tbirth[nxt + 1] = times[k]
            live[i] = nxt
            live[nl] = nxt + 1
            nl += 1; nxt += 2
        else:
            nl -= 1
            live[i] = live[nl]                   # swap-with-last removal, O(1)
    return parent, tbirth, live[:nl].copy()


# --------------------------------------------------------------------------
# Step 5 -- the reconstructed tree.  The bug-prone part; 02 hammers it.
# --------------------------------------------------------------------------


@dataclass
class ReconTree:
    n_tips: int
    parent: np.ndarray          # length 2*n_tips-1; leaves 0..n_tips-1, root last, root parent -1
    time: np.ndarray            # node times measured FORWARD from the origin t=0
    T: float
    n_ancestors: int = 0        # nodes on the sampled tips' ancestor paths, BEFORE suppression.
                                # (n_ancestors - 1)/(n_tips - 1) is the mean unary-chain length --
                                # i.e. how much work degree-2 suppression actually did.

    @property
    def branching_times(self) -> np.ndarray:
        r"""Internal node ages measured BEFORE THE PRESENT -- the quantity G is a CDF for."""
        return self.T - self.time[self.n_tips:]

    def newick(self, labels=None) -> str:
        lab = labels or [f"t{i}" for i in range(self.n_tips)]
        kids: dict[int, list[int]] = {}
        for x in range(self.parent.size - 1):
            kids.setdefault(int(self.parent[x]), []).append(x)
        root = self.parent.size - 1

        def emit(x):
            if x < self.n_tips:
                return f"{lab[x]}:{self.T - self.time[int(self.parent[x])]:.8g}"
            inner = ",".join(emit(c) for c in kids[x])
            if x == root:
                return f"({inner})"
            return f"({inner}):{self.time[x] - self.time[int(self.parent[x])]:.8g}"

        return emit(root) + ";"


def reconstruct(parent, tbirth, live, sampled, T) -> ReconTree | None:
    r"""Delete every node with no sampled descendant, then suppress degree-2 nodes.

    A birth at which only ONE child left sampled descendants is not a visible
    branching -- it is spliced through.  That suppression is the whole reason
    the reconstructed tree is not itself a birth-death process.
    """
    sampled = np.asarray(sampled, dtype=np.int64)
    n_tips = int(sampled.size)
    if n_tips < 2:
        return None

    # ancestors of the sampled tips: walk up, stop the moment we revisit
    seen = np.zeros(parent.size, dtype=bool)
    order: list[int] = []
    for tip in sampled:
        x = int(tip)
        while x != -1 and not seen[x]:
            seen[x] = True
            order.append(x)
            x = int(parent[x])
    anc = np.sort(np.asarray(order, dtype=np.int64))       # ids increase down the tree

    # how many sampled tips sit below each ancestor (children first => reversed)
    ndesc = np.zeros(parent.size, dtype=np.int64)
    ndesc[sampled] = 1
    for x in anc[::-1]:
        p = int(parent[x])
        if p != -1:
            ndesc[p] += ndesc[x]

    kids: dict[int, list[int]] = {}
    for x in anc:
        p = int(parent[x])
        if p != -1:
            kids.setdefault(p, []).append(int(x))

    tip_slot = {int(t): i for i, t in enumerate(sampled)}

    def descend(x: int) -> tuple[int, bool]:
        """Follow the unary chain from x to the next branch point or sampled tip."""
        while True:
            if x in tip_slot:
                return x, True
            c = kids.get(x, [])
            if len(c) == 2:
                return x, False
            if len(c) == 1:
                x = c[0]
                continue
            raise RuntimeError(f"ancestor {x} has no sampled descendant -- pruning is wrong")

    n_nodes = 2 * n_tips - 1
    out_parent = np.full(n_nodes, -1, np.int64)
    out_time = np.empty(n_nodes, np.float64)
    out_time[:n_tips] = T

    root, is_tip = descend(int(anc[0]))
    if is_tip:
        return None
    slot = n_nodes - 1
    stack = [(root, -1)]
    while stack:
        x, pslot = stack.pop()
        s = slot; slot -= 1
        out_parent[s] = pslot
        out_time[s] = tbirth[kids[x][0]]                   # x divided when its children were born
        for c in kids[x]:
            y, y_is_tip = descend(c)
            if y_is_tip:
                out_parent[tip_slot[y]] = s
            else:
                stack.append((y, s))

    return ReconTree(n_tips, out_parent, out_time, T, int(anc.size))


# --------------------------------------------------------------------------
# The two entry points.
# --------------------------------------------------------------------------


def simulate_tree_naive(b, delta, T, rho, rng):
    r"""Oracle.  One event at a time, nothing pruned until the end, no conditioning.

    Returns (tree_or_None, n_T).  Use only at small n -- it is O(all events) in
    both time and memory by design, so that there is nothing clever to be wrong.
    """
    total = b + delta
    p_birth = b / total
    parent, tbirth, live, t = [-1], [0.0], [0], 0.0
    while live:
        n = len(live)
        t += float(rng.exponential(1.0 / (n * total)))
        if t > T:
            break
        i = int(rng.integers(n))
        x = live[i]
        if rng.random() < p_birth:
            c1 = len(parent); parent.append(x); tbirth.append(t)
            c2 = len(parent); parent.append(x); tbirth.append(t)
            live[i] = c1; live.append(c2)
        else:
            live[i] = live[-1]; live.pop()
    if not live:
        return None, 0
    n_T = len(live)
    keep = [x for x in live if rng.random() < rho]
    if len(keep) < 2:
        return None, n_T
    return reconstruct(np.asarray(parent), np.asarray(tbirth), live, keep, T), n_T


@dataclass
class Stats:
    attempts: int = 0
    extinct: int = 0
    events_scan: int = 0        # pass-1 lineage-events: the term that dominates cost
    events_build: int = 0       # pass-2 lineage-events: paid once per acceptance
    n_T: int = 0                # survivors at T in the accepted replicate
    peak_live: int = 0
    seconds: float = 0.0
    accepted: bool = False
    rejects: list = field(default_factory=list)   # sampled tip counts of rejected attempts


def simulate_tree(b, delta, T, rho, n_target, seed,
                  max_attempts=10 ** 8, tol=0, keep_rejects=False):
    r"""Rejection-sample a reconstructed tree with |tips - n_target| <= tol.

    tol = 0 is exact conditioning on the observed clone size, which is what
    adjudication against Park needs and what the cost measurement must report.
    tol > 0 is the tolerance-window discount -- a legitimate softening, but it
    changes what the conditional means and must be declared wherever it is used.
    """
    st = Stats()
    t0 = time.perf_counter()
    for a in range(max_attempts):
        st.attempts += 1
        c_ss, s_ss, l_ss = np.random.SeedSequence([seed, a]).spawn(3)

        tr = simulate_counts(b, delta, T, np.random.default_rng(c_ss))
        st.events_scan += tr.n_events
        if tr.n_T == 0:
            st.extinct += 1
            continue

        rng_s = np.random.default_rng(s_ss)
        k = int(rng_s.binomial(tr.n_T, rho))
        if abs(k - n_target) > tol or k < 2:
            if keep_rejects:
                st.rejects.append(k)
            continue

        # accepted: replay the identical count stream, this time keeping it
        tr2 = simulate_counts(b, delta, T, np.random.default_rng(c_ss), record=True)
        assert tr2.n_T == tr.n_T and tr2.n_events == tr.n_events, "replay diverged"
        parent, tbirth, livearr = build_genealogy(tr2.steps, tr2.times,
                                                  np.random.default_rng(l_ss))
        idx = rng_s.choice(livearr.size, size=k, replace=False, shuffle=False)
        tree = reconstruct(parent, tbirth, livearr, livearr[idx], T)

        st.events_build = tr2.n_events
        st.n_T = tr.n_T
        st.peak_live = int(livearr.size)
        st.accepted = tree is not None
        st.seconds = time.perf_counter() - t0
        return tree, st

    st.seconds = time.perf_counter() - t0
    return None, st
