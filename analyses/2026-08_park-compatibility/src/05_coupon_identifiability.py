#!/usr/bin/env python3
"""Is m recoverable from s? -- the coupon-inversion identifiability check.

m  = number of independent write events at a site within a prefix clade (unobservable)
s  = number of DISTINCT symbols those events produced (directly countable)

Under SciPhy's jump chain the m events are iid draws from xi (Sec 1a.3), so

    E[s | m] = sum_i [ 1 - (1 - xi_i)^m ]

is monotone increasing in m and therefore invertible: measure s, read off m-hat,
and (m-hat - s) estimates the number of homoplastic recurrences.

The question this script answers is over WHAT RANGE OF m that inversion actually
carries information. Once m is large enough that nearly every symbol has been
drawn at least once, s saturates at M_obs and stops responding to m -- so the
estimator fails exactly in the high-homoplasy regime we care about.

To make that honest we also need the sampling spread of s:

    Var[s | m] = sum_i p_i(1-p_i)
               + sum_{i != j} [ (1-xi_i-xi_j)^m - (1-xi_i)^m (1-xi_j)^m ]
    where p_i = 1 - (1-xi_i)^m

(the cross term is the occupancy covariance, which is negative -- seeing one
symbol makes another marginally less likely at fixed m). A true m is then only
distinguishable from the m' values whose expected s lies within +-1 sd, so
inverting the band [E-sd, E+sd] gives the recoverable interval [m_lo, m_hi].

Outputs results/coupon_identifiability.json and figures/coupon_identifiability.png
"""
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
MAIN = "Mouse1"
MIN_COUNT = 1000     # symbols below this are the artifact tail: 97 of 197 symbols
                     # carry <0.05% of edits between them, and they are almost
                     # certainly sequencing error rather than real pegRNAs. They
                     # matter enormously here -- being rare, they keep s growing
                     # and make the estimator LOOK well-conditioned at large m.
BIRTHDAY_Q = None                  # filled from the data

# Park clade sizes worth marking (notes 1735: ~74 cells/clone; subclones 313-10,506)
MARKS = [(74, "clone\n~74 cells"), (313, "subclone\nmin 313"), (10506, "subclone\nmax 10,506")]

xi_all = json.loads((RES / "xi_vectors.json").read_text())


def moments(xi, m, want_sd=True):
    """E[s|m] and sd[s|m] for m iid draws from xi. m may be an array.

    want_sd=False skips the O(M^2) occupancy-covariance term. Needed for the
    untrimmed vector, which now carries ~10^4 values (every value seen >=2x), where
    the pair array would be terabytes. Only the mean is used for that curve."""
    xi = np.asarray(xi)
    m = np.atleast_1d(np.asarray(m, dtype=float))
    # a[k,i] = (1 - xi_i)^m_k
    a = np.power(1.0 - xi[None, :], m[:, None])
    p = 1.0 - a
    exp_s = p.sum(axis=1)
    if not want_sd:
        return exp_s, np.zeros_like(exp_s)

    pair = 1.0 - xi[None, :, None] - xi[None, None, :]          # 1 - xi_i - xi_j
    B = np.power(np.clip(pair, 0.0, None), m[:, None, None])
    cov = B - a[:, :, None] * a[:, None, :]
    diag = np.power(np.clip(1.0 - 2 * xi, 0.0, None), m[:, None]) - a ** 2
    cross = cov.sum(axis=(1, 2)) - diag.sum(axis=1)             # i != j only
    var = (p * (1 - p)).sum(axis=1) + cross
    return exp_s, np.sqrt(np.clip(var, 0.0, None))


def invert(m_grid, exp_s, target):
    """Smallest m whose E[s|m] reaches `target`; np.inf if never (saturated)."""
    out = np.full_like(target, np.inf, dtype=float)
    ok = target <= exp_s[-1]
    out[ok] = np.interp(target[ok], exp_s, m_grid)
    out[target <= exp_s[0]] = m_grid[0]
    return out


# dense log grid for inversion
m_grid = np.unique(np.round(np.logspace(0, 6, 4000)).astype(int)).astype(float)

def trim(d):
    """script 04 now emits both: xi = kept alphabet, xi_untrimmed = every value seen >=2x"""
    xi_real = np.array(list(d["xi"].values())); xi_real = xi_real / xi_real.sum()
    xi_full = np.array(list(d["xi_untrimmed"].values())); xi_full = xi_full / xi_full.sum()
    return xi_real, xi_full

curves = {}
for tbl, d in xi_all.items():
    xi_real, xi_full = trim(d)
    e, sd = moments(xi_real, m_grid)
    e_full, _ = moments(xi_full, m_grid, want_sd=False)
    curves[tbl] = dict(xi=xi_real, exp_s=e, sd=sd, exp_s_full=e_full,
                       M_real=len(xi_real), M_obs=d["M_obs"],
                       q=float((xi_real**2).sum()),
                       share=float(xi_full[np.round(xi_full*d["n_edits"])>=MIN_COUNT].sum()))

c = curves[MAIN]
q = c["q"]
BIRTHDAY_Q = np.sqrt(2 / q)

m_lo = invert(m_grid, c["exp_s"], c["exp_s"] - c["sd"])
m_hi = invert(m_grid, c["exp_s"], c["exp_s"] + c["sd"])
finite = np.isfinite(m_hi)

def first_failure(ok):
    """Largest m below which `ok` holds continuously (not merely the last True)."""
    bad = np.flatnonzero(~ok)
    if bad.size == 0:
        return m_grid[-1]
    return m_grid[bad[0] - 1] if bad[0] > 0 else np.nan

m_break = first_failure(finite)                  # above this, s has saturated
ratio = np.where(finite, m_hi / np.maximum(m_lo, 1e-9), np.inf)
m_fac2 = first_failure(ratio <= 2.0)             # recoverable to better than 2x

summary = dict(
    table=MAIN, q=q, M_obs=c["M_real"], birthday_threshold=float(BIRTHDAY_Q),
    m_at_half_alphabet=float(np.interp(c["M_real"] / 2, c["exp_s"], m_grid)),
    m_at_90pct_alphabet=float(np.interp(0.9 * c["M_real"], c["exp_s"], m_grid)),
    m_recoverable_within_factor2=float(m_fac2),
    m_upper_bound_lost_above=float(m_break),
    E_s_at=({str(k): float(np.interp(k, m_grid, c["exp_s"])) for k in [11, 74, 313, 10506]}),
)
print(json.dumps(summary, indent=2))
RES.mkdir(exist_ok=True)
(RES / "coupon_identifiability.json").write_text(json.dumps(summary, indent=2))

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
GRID, AXIS, BLUE, ORANGE = "#e1e0d9", "#c3c2b7", "#2a78d6", "#eb6834"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": INK2, "ytick.color": INK2, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
})

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.2, 4.5))

# --- Panel A: the saturation curve
axA.plot(m_grid, c["exp_s_full"], color=AXIS, lw=1.6, zorder=1)
axA.text(2.4e4, 118, "including the rare\nartifact tail\n(197 symbols)",
         color=INK2, fontsize=8, va="center")
axA.text(2.4e4, 78, "real alphabet\n(100 symbols,\n99.95% of edits)",
         color=BLUE, fontsize=8, va="center", fontweight="bold")
axA.fill_between(m_grid, c["exp_s"] - c["sd"], c["exp_s"] + c["sd"],
                 color=BLUE, alpha=0.18, lw=0, zorder=2)
axA.plot(m_grid, c["exp_s"], color=BLUE, lw=2.0, zorder=3)
axA.axhline(c["M_real"], color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=1)
axA.text(1.15, c["M_real"] - 5, f"$M={c['M_real']}$ real symbols — all of them seen",
         color=INK2, fontsize=8, va="top")
axA.axvline(BIRTHDAY_Q, color=ORANGE, lw=1.4, ls=(0, (5, 2)), zorder=4)
axA.text(BIRTHDAY_Q * 1.18, 8, f"$m^*={BIRTHDAY_Q:.0f}$\nbirthday\nthreshold",
         color=ORANGE, fontsize=8, va="bottom")
for x, lab in MARKS:
    axA.axvline(x, color=AXIS, lw=0.9, ls=":", zorder=1)
    axA.text(x * 1.1, 30, lab, color=INK2, fontsize=7.5, va="top")
axA.set_xscale("log")
axA.set_xlim(1, 1e5); axA.set_ylim(0, 205)
axA.set_xlabel("$m$  — independent write events in the prefix clade")
axA.set_ylabel("$\\mathbb{E}[s\\,|\\,m]$  — distinct symbols seen")
axA.set_title("A.  Against the real alphabet, $s$ saturates inside the subclone range",
              loc="left", fontsize=10, color=INK, pad=8)

# --- Panel B: what m can be recovered
band_hi = np.where(finite, m_hi, m_grid[-1] * 10)
axB.fill_between(m_grid, m_lo, band_hi, color=BLUE, alpha=0.18, lw=0, zorder=2)
axB.plot(m_grid, m_lo, color=BLUE, lw=1.4, zorder=3)
axB.plot(m_grid[finite], m_hi[finite], color=BLUE, lw=1.4, zorder=3)
axB.plot(m_grid, m_grid, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=4)
axB.text(2.2, 1.25, "truth (identity)", color=INK2, fontsize=8, rotation=34)
if np.isfinite(m_break):
    axB.axvline(m_break, color=ORANGE, lw=1.4, ls=(0, (5, 2)), zorder=5)
    axB.text(m_break * 1.25, 2.2, f"above $m\\approx{m_break:,.0f}$\nno upper bound:\n$s$ has saturated",
             color=ORANGE, fontsize=8, va="bottom")
if np.isfinite(m_fac2):
    axB.axvspan(1, m_fac2, color=BLUE, alpha=0.055, lw=0, zorder=0)
    axB.text(1.15, 3.2e5, f"recoverable to\nbetter than 2×\n($m\\lesssim{m_fac2:,.0f}$)",
             color=BLUE, fontsize=8, va="top")
axB.set_xscale("log"); axB.set_yscale("log")
axB.set_xlim(1, 1e5); axB.set_ylim(1, 1e6)
axB.set_xlabel("true $m$")
axB.set_ylabel("$m$ consistent with the observed $s$  (±1 sd)")
axB.set_title("B.  Above the saturation point $m$ is unbounded from above",
              loc="left", fontsize=10, color=INK, pad=8)

fig.suptitle(f"Recovering write-event count $m$ from distinct-symbol count $s$   ·   "
             f"Park {MAIN}, $q={q:.4f}$, effective alphabet $1/q={1/q:.0f}$, real alphabet $M={c['M_real']}$",
             fontsize=10.5, color=INK, x=0.012, ha="left", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.94))
FIG.mkdir(exist_ok=True)
fig.savefig(FIG / "coupon_identifiability.png", dpi=200, facecolor=SURFACE)
print(f"\nwrote {FIG/'coupon_identifiability.png'}")
