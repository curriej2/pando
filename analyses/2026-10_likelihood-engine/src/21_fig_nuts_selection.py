"""Teaching figure for Betancourt §5.2 / App. A (notes §G.4.31-32): why NUTS picks a state from the whole trajectory
with probability proportional to exp(-H), and why that is robust. 1-D standard normal posterior, unit mass.

Trajectory construction used throughout (Betancourt A.2, "static multinomial" HMC, the core of NUTS): given the
current state z0, choose uniformly where z0 sits in a trajectory of L+1 states, integrate backwards and forwards
with leapfrog accordingly, then pick ONE state from the list.

a  One trajectory (the worked 8-state example): the chance each state is picked, exp(-H)-weighted vs uniform.
b  The crowd test: 400,000 chains already distributed exactly as the target (q and p drawn from it); each makes
   ONE transition with a deliberately coarse step (eps = 1.2, 4 steps). Pick by exp(-H) -> crowd stays correct;
   pick uniformly -> crowd is distorted; take the endpoint with no correction -> distorted.
c  Robust to trajectory length: lag-1 correlation of successive samples vs number of steps (eps = 0.3, one full
   swing = 21 steps). Endpoint + accept/reject swings between +1 and -1 (the periodicity traps); picking from
   the whole trajectory stays low.
(A panel d comparing "fraction of transitions that do not move" was built and dropped: re-picking the current
state has baseline chance ~1/(number of states), so for short trajectories the whole-trajectory pick stays put
MORE often than endpoint + accept/reject (11-64% vs 0.4-39% at eps 0.2-1.8). Its advantages are exactness
(panel b) and robustness to trajectory length (panel c), not fewer stays. Notes §G.4.32.)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, MUTED, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#a3a29c", "#e4e3df", "#fcfcfb", "#efeeea"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})
rng = np.random.default_rng(21102026)


def subtitle(a, txt):
    a.text(0, 1.03, txt, transform=a.transAxes, fontsize=8.4, color=INK2, va="bottom")


Hf = lambda q, p: 0.5 * q * q + 0.5 * p * p


def lf(q, p, eps):                                     # vectorised leapfrog step, force = -q
    p = p - 0.5 * eps * q; q = q + eps * p; p = p - 0.5 * eps * q
    return q, p


def build(q0, p0, eps, L, k):
    """States of a trajectory with L+1 states in which z0 sits at index k (k steps backwards, L-k forwards).
    Vectorised over chains: q0, p0, k are arrays. Returns Q, P of shape (L+1, n)."""
    n = q0.shape[0]
    Q = np.empty((L + 1, n)); P = np.empty((L + 1, n))
    # backwards part: integrate with -eps, fill indices k-1 ... 0
    qb, pb = q0.copy(), p0.copy()
    for j in range(1, L + 1):
        qb, pb = lf(qb, pb, -eps)
        idx = k - j
        m = idx >= 0
        Q[np.where(m, idx, 0), np.arange(n)] = np.where(m, qb, Q[np.where(m, idx, 0), np.arange(n)])
        P[np.where(m, idx, 0), np.arange(n)] = np.where(m, pb, P[np.where(m, idx, 0), np.arange(n)])
    Q[k, np.arange(n)] = q0; P[k, np.arange(n)] = p0
    qf, pf = q0.copy(), p0.copy()
    for j in range(1, L + 1):
        qf, pf = lf(qf, pf, eps)
        idx = k + j
        m = idx <= L
        Q[np.where(m, idx, L), np.arange(n)] = np.where(m, qf, Q[np.where(m, idx, L), np.arange(n)])
        P[np.where(m, idx, L), np.arange(n)] = np.where(m, pf, P[np.where(m, idx, L), np.arange(n)])
    return Q, P


def pick(Q, P, rule):
    n = Q.shape[1]
    if rule == "uniform":
        j = rng.integers(0, Q.shape[0], n)
    else:
        W = np.exp(-(Hf(Q, P) - Hf(Q, P).min(axis=0)))
        W = W / W.sum(axis=0)
        u = rng.random(n); j = (np.cumsum(W, axis=0) < u).sum(axis=0)
    return Q[j, np.arange(n)]


fig = plt.figure(figsize=(16, 11))
gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.22)

# ---------------------------------------------------------------- a: one trajectory's weights
q0, p0, eps = 1.2, -0.8, 0.6
traj = []
q, p = q0, p0
back = []
for i in range(3):
    q, p = lf(q, p, -eps); back.append((q, p))
traj = back[::-1] + [(q0, p0)]
q, p = q0, p0
for i in range(4):
    q, p = lf(q, p, eps); traj.append((q, p))
Hs = np.array([Hf(*z) for z in traj]); w = np.exp(-Hs) / np.exp(-Hs).sum()
ax = fig.add_subplot(gs[0, 0])
lab = [f"{k-3:+d}" for k in range(len(traj))]
bars = ax.bar(range(len(traj)), w, color=[C[0] if k != 3 else INK for k in range(len(traj))], width=0.7)
ax.axhline(1 / len(traj), color=C[1], lw=1.6, ls="--", label=f"uniform pick (1/{len(traj)} = {1/len(traj):.3f})")
for k, (wk, h) in enumerate(zip(w, Hs)):
    ax.text(k, wk + 0.0015, f"ΔH {h - Hs[3]:+.3f}", ha="center", fontsize=7.6, color=INK2)
ax.set_xticks(range(len(traj))); ax.set_xticklabels(lab)
ax.set_ylim(0.11, 0.137)
ax.set_xlabel("state along the trajectory  (0 = the current state, black)")
ax.set_ylabel("chance of being picked")
ax.set_title("a  One trajectory: pick probability ∝ e^(−H)", loc="left", pad=28)
subtitle(ax, "energy errors ΔH relative to the start are small, so the pick is nearly uniform;\n"
         "a state where leapfrog's energy went UP (less probable under the target) is picked slightly less")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)

# ---------------------------------------------------------------- b: the crowd test
n = 400_000; epsb, Lb = 1.2, 4
qc = rng.standard_normal(n); pc = rng.standard_normal(n)        # crowd already exactly at the target
k = rng.integers(0, Lb + 1, n)
Q, P = build(qc, pc, epsb, Lb, k)
q_w = pick(Q, P, "weighted"); q_u = pick(Q, P, "uniform")
qe, pe = qc.copy(), pc.copy()
for i in range(Lb):
    qe, pe = lf(qe, pe, epsb)
axb = fig.add_subplot(gs[0, 1])
xg = np.linspace(-4.5, 4.5, 400)
bins = np.linspace(-4.5, 4.5, 61)
axb.plot(xg, stats.norm.pdf(xg), color=INK, lw=2.4, label="target posterior (where the crowd started)")
for arr, col, labl in [(q_w, C[0], f"pick ∝ e^(−H):  SD {q_w.std():.3f}"),
                       (q_u, C[1], f"pick uniformly:  SD {q_u.std():.3f}"),
                       (qe, C[2], f"endpoint, no correction:  SD {qe.std():.3f}")]:
    h, e = np.histogram(arr, bins=bins, density=True)
    axb.step(e[:-1], h, where="post", color=col, lw=1.8, label=labl)
axb.set_xlabel("q after ONE transition"); axb.set_ylabel("density")
axb.set_title("b  The crowd test: which rule keeps a correct crowd correct?", loc="left", pad=28)
subtitle(axb, "400,000 chains start EXACTLY at the target (true SD 1); each makes one transition with a coarse\n"
         "step (ε = 1.2, 4 steps) so leapfrog's energy errors are large. Only the e^(−H) pick leaves the crowd as it was")
axb.legend(frameon=False, fontsize=8, loc="upper right")
axb.grid(axis="y", color=GRID, lw=0.6); axb.set_axisbelow(True)

# ---------------------------------------------------------------- c: robustness to trajectory length
epsc = 0.3; Ls = np.arange(1, 43); n_ch = 4000; n_it = 60
r_end, r_mult = [], []
for L in Ls:
    q = rng.standard_normal(n_ch); prev = []; nxt = []
    qa = q.copy()                                     # endpoint + accept/reject
    for it in range(n_it):
        p = rng.standard_normal(n_ch); H0 = Hf(qa, p); qn, pn = qa.copy(), p.copy()
        for i in range(L):
            qn, pn = lf(qn, pn, epsc)
        acc = rng.random(n_ch) < np.exp(np.minimum(0, H0 - Hf(qn, pn)))
        qnew = np.where(acc, qn, qa)
        if it >= 10:
            prev.append(qa); nxt.append(qnew)
        qa = qnew
    r_end.append(np.corrcoef(np.concatenate(prev), np.concatenate(nxt))[0, 1])
    qm = q.copy(); prev = []; nxt = []
    for it in range(n_it):
        p = rng.standard_normal(n_ch); kk = rng.integers(0, L + 1, n_ch)
        Qm, Pm = build(qm, p, epsc, L, kk)
        qnew = pick(Qm, Pm, "weighted")
        if it >= 10:
            prev.append(qm); nxt.append(qnew)
        qm = qnew
    r_mult.append(np.corrcoef(np.concatenate(prev), np.concatenate(nxt))[0, 1])
axc = fig.add_subplot(gs[1, :])
axc.axhline(0, color=INK2, lw=0.8)
axc.plot(Ls, r_end, "-o", color=C[1], ms=3.5, label="endpoint + accept/reject")
axc.plot(Ls, r_mult, "-o", color=C[0], ms=3.5, label="pick from the whole trajectory ∝ e^(−H)")
for frac, txt in [(0.5, "half a swing:\nlands on the mirror image"), (1.0, "full swing:\nback to the start")]:
    x = frac * 2 * np.pi / epsc
    axc.axvline(x, color=MUTED, ls=":", lw=1)
    axc.text(x + 0.4, 0.62 if frac == 1.0 else -0.62, txt, fontsize=8, color=INK2)
axc.set_ylim(-1.05, 1.05)
axc.set_xlabel("number of leapfrog steps per trajectory  (ε = 0.3; one full swing ≈ 21 steps)")
axc.set_ylabel("correlation of successive samples")
axc.set_title("c  Robust to trajectory length", loc="left", pad=28)
subtitle(axc, "0 = successive samples independent (ideal); ±1 = no new information. The endpoint rule is best at exactly a quarter swing "
         "(≈5 steps) and terrible at half and full swings;\nthe whole-trajectory pick needs a longer trajectory (≈¾ swing) to decorrelate, "
         "then STAYS decorrelated at any length — useful because the right length is never known in advance")
axc.legend(frameon=False, fontsize=8, loc="lower left")
axc.grid(color=GRID, lw=0.6); axc.set_axisbelow(True)

fig.suptitle("Why NUTS picks one state from the whole trajectory with probability ∝ e^(−H) — and why that is robust",
             x=0.06, y=1.0, ha="left", fontsize=13, fontweight="bold")
fig.text(0.06, 0.965, "1-D standard normal posterior.  e^(−H) is the target probability of a (position, momentum) state; "
         "picking in proportion to it keeps an already-correct crowd of chains correct.", fontsize=9, color=INK2)
fig.savefig(OUT / "fig_hmc_nuts_selection.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("a: pick probs", np.round(w, 3))
print(f"b: SD after one transition: weighted {q_w.std():.4f}, uniform {q_u.std():.4f}, endpoint-uncorrected {qe.std():.4f}; "
      f"weighted p-check SD start {qc.std():.4f}")
print("c: corr endpoint at L=10,11,21 ->", [round(r_end[i - 1], 3) for i in (10, 11, 21)],
      " multinomial ->", [round(r_mult[i - 1], 3) for i in (10, 11, 21)], " max |r| mult", round(np.max(np.abs(r_mult)), 3))
print("wrote", OUT / "fig_hmc_nuts_selection.png")
