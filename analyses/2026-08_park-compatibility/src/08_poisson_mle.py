#!/usr/bin/env python3
r"""
================================================================================
 The set-dependent (Poissonised MLE) estimator for m
================================================================================

WHAT WE ARE ESTIMATING
----------------------
At one node of the prefix trie -- a (clone, tape, prefix p') triple -- some number
m of INDEPENDENT WRITE EVENTS occurred at the next site among the carriers of p'.
Each event drew a symbol independently from the measured insert-probability vector
xi = (xi_1, ..., xi_M).  We cannot see the events.  We see only the SET

    A  =  { symbols that appear at that site among the carriers }        s = |A|

m is not s: whenever two events happened to draw the same symbol (a homoplasy),
two events collapse into one observed symbol.  The whole point is that

    (number of homoplastic recurrences)  =  m - s

so estimating m estimates the homoplasy.  Note carefully that s counts DISTINCT
SYMBOLS, not cells: a hundred cells that inherited one symbol from one ancestor
contribute one event and one symbol, which is why relatedness between cells does
not bias any of this and no tree is needed.


WHY THE COUNT s IS NOT ENOUGH  (the point of this script)
---------------------------------------------------------
Script 06 estimated m by matching the count alone -- solving E[s|m] = s_obs.  That
throws away WHICH symbols were seen, and that is real information:

    * If A is made of RARE symbols, a large m is implausible.  With many draws you
      would almost certainly have hit some common symbol too, and you didn't.
    * If A is made of COMMON symbols, a large m is perfectly consistent -- the
      common symbols soak up the draws.

So two nodes with identical s but different A should get different m-hat.


THE EXACT LIKELIHOOD, AND WHY IT IS UNUSABLE
--------------------------------------------
Write W_B = sum_{i in B} xi_i for the total probability mass of a symbol set B.
For exactly m draws, the distinct-symbol set equals A exactly when every draw
lands inside A AND every member of A is hit at least once.  Inclusion-exclusion
over which members get missed gives

    P(A | m)  =  sum_{B subset of A}  (-1)^(s - |B|)  *  W_B^m                (1)

Check, A = {i}:   only B = {} and B = {i} contribute  ->  xi_i^m.  Correct: all m
                  draws had to be symbol i.
Check, A = {i,j}: (xi_i + xi_j)^m - xi_i^m - xi_j^m.  Correct: all draws inside
                  the pair, minus the all-i case, minus the all-j case.

Equation (1) is exact and useless: it has 2^s terms, and s reaches ~60 here.
The trouble is that with a FIXED number of draws the symbol counts are coupled --
they must sum to m -- which is what forces the inclusion-exclusion.


POISSONISATION REMOVES THE COUPLING
-----------------------------------
Replace "exactly m draws" by "N ~ Poisson(m) draws".  Standard result (Poisson
thinning / the Poissonisation theorem): if N ~ Poisson(m) and each of the N items
is independently assigned to category i with probability xi_i, then the category
counts are INDEPENDENT Poissons,

    N_i ~ Poisson(m * xi_i),   independently across i.                        (2)

Independence is the entire payoff.  Each symbol is now its own coin flip:

    P(symbol i is seen)      = P(N_i >= 1) = 1 - exp(-m * xi_i)
    P(symbol i is not seen)  = P(N_i  = 0) =     exp(-m * xi_i)

and the observed set A is just "which coins came up heads", so the likelihood
factorises completely:

    P(A | m)  =  prod_{i in A} [1 - exp(-m*xi_i)]  *  prod_{i not in A} exp(-m*xi_i)

Taking logs, and using sum_{i not in A} xi_i = 1 - W_A:

                                                                              (3)
    log L(m)  =   sum_{i in A} log( 1 - exp(-m*xi_i) )   -   m * (1 - W_A)
                  \_________________________________/       \____________/
                   every SEEN symbol pushes m UP:             every UNSEEN symbol
                   you needed enough draws to hit it          pushes m DOWN, in
                                                              proportion to its mass

The second term is exactly the "set-dependent" correction.  The penalty on large m
is proportional to the mass NOT observed.  Same s, but an A made of common symbols
=> small (1 - W_A) => weak penalty => larger m-hat.  That is the behaviour we want
and the count-only estimator cannot produce.


SOLVING IT
----------
Differentiate (3).  Using d/dm log(1 - e^{-m*xi}) = xi * e^{-m*xi} / (1 - e^{-m*xi})
                                                  = xi / (e^{m*xi} - 1):

    d log L / dm  =   sum_{i in A}  xi_i / ( exp(m*xi_i) - 1 )   -   (1 - W_A)  (4)

Each term xi/(e^{m*xi} - 1) is strictly DECREASING in m, so (4) is strictly
decreasing, so log L is strictly CONCAVE and the maximum is unique.  Limits:

    m -> 0+ :  xi/(e^{m*xi}-1) -> 1/m, so the sum -> s/m -> +infinity.
    m -> inf:  the sum -> 0, so the derivative -> -(1 - W_A).

Hence a unique finite root exists iff W_A < 1, and it is found by bisection on (4).

SATURATION IS HANDLED HONESTLY.  If A is the entire alphabet then W_A = 1, the
penalty term vanishes, the derivative is positive everywhere, and m-hat = infinity.
That is correct: once every symbol has been seen the data genuinely cannot bound m
from above.  The count-only estimator has the same limit but reaches it silently.


RELATION TO THE COUNT-ONLY ESTIMATOR
------------------------------------
Script 06 solved  sum_i [1 - (1-xi_i)^m]  =  s_obs   (method of moments on s).
Its Poissonised analogue is  sum_i [1 - exp(-m*xi_i)] = s_obs.  Both match only
the EXPECTED count and are blind to which symbols are missing.  Since A is on
average a typical draw, the two estimators agree on average; the MLE wins on
individual nodes, and it wins most when A is atypical -- which is precisely the
large-s, high-homoplasy nodes that carry most of the total.


WHAT THIS SCRIPT DOES
---------------------
 1. Monte-Carlo validation: simulate under the TRUE fixed-m model (not the
    Poissonised one), and compare bias/spread of the MLE against the count-only
    estimator.  This simultaneously measures what Poissonisation costs.
 2. Re-walks the prefix trie and computes both estimators for every node.

Usage:  08_poisson_mle.py --validate            (MC only)
        08_poisson_mle.py Mouse1 Mouse2 ...     (MC, then those tables)
================================================================================
"""
import csv, re, json, gzip, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/data1/choij10/justin/pando/data/cancer_metastasis")
RES = ROOT / "results"
SYM = re.compile(r"^[ACGT]{4}GGA$")
MIN_COUNT, MIN_CARRIERS = 1000, 3
M_CAP = 1e7            # report saturation rather than infinity


# ----------------------------------------------------------------- estimators
def mle_m(xi_A, W_A, cap=M_CAP):
    """Maximise (3) by bisecting its derivative (4).  xi_A = xi values of seen symbols."""
    rest = 1.0 - W_A
    if rest <= 1e-12:                      # A is the whole alphabet: unbounded
        return cap, True
    def dll(m):
        # xi / (exp(m*xi) - 1), guarded against overflow for large m*xi
        z = m * xi_A
        with np.errstate(over="ignore"):
            return float(np.sum(xi_A / np.expm1(np.minimum(z, 700.0)))) - rest
    lo, hi = 1e-9, 1.0
    while dll(hi) > 0:                     # expand until the derivative turns negative
        hi *= 4.0
        if hi > cap:
            return cap, True
    for _ in range(200):                   # bisection: ~1e-60 relative, overkill but cheap
        mid = 0.5 * (lo + hi)
        if dll(mid) > 0: lo = mid
        else:            hi = mid
        if hi - lo < 1e-10 * hi: break
    return 0.5 * (lo + hi), False


def make_moment_inverter(xi):
    """Script 06's estimator: solve sum_i [1-(1-xi_i)^m] = s  (multinomial form)."""
    grid = np.unique(np.round(np.logspace(0, 7, 5000)).astype(int)).astype(float)
    exp_s = (1.0 - np.power(1.0 - xi[None, :], grid[:, None])).sum(axis=1)
    return lambda s: float(np.interp(float(s), exp_s, grid))


# ----------------------------------------------------------------- validation
def validate(xi, reps=3000, seed=0):
    """Simulate under FIXED m (the true model) and score both estimators."""
    rng = np.random.default_rng(seed)
    xi = xi / xi.sum()
    moment = make_moment_inverter(xi)
    cache = {}
    print("\nMONTE-CARLO VALIDATION  (truth = fixed-m multinomial draws)")
    print("  simulating exactly m draws, then estimating m back from the symbol set\n")
    print(f"  {'true m':>7} {'mean s':>7} | {'MLE median':>11} {'MLE bias%':>10} {'MLE IQR%':>9}"
          f" | {'count median':>13} {'count bias%':>12} {'count IQR%':>11}")
    out = {}
    for m_true in [5, 10, 20, 50, 100, 200, 500]:
        r = reps if m_true <= 200 else reps // 3
        mle_v, mom_v, s_v = [], [], []
        for _ in range(r):
            draws = rng.choice(len(xi), size=m_true, p=xi)
            idx = np.unique(draws)
            key = idx.tobytes()
            if key not in cache:
                xa = xi[idx]
                cache[key] = mle_m(xa, float(xa.sum()))[0]
            mle_v.append(cache[key]); s_v.append(len(idx)); mom_v.append(moment(len(idx)))
        mle_v, mom_v = np.array(mle_v), np.array(mom_v)
        f = lambda v: (np.median(v), 100*(np.median(v)-m_true)/m_true,
                       100*(np.percentile(v,75)-np.percentile(v,25))/m_true)
        a, b = f(mle_v), f(mom_v)
        print(f"  {m_true:>7} {np.mean(s_v):>7.1f} | {a[0]:>11.1f} {a[1]:>+9.1f}% {a[2]:>8.1f}%"
              f" | {b[0]:>13.1f} {b[1]:>+11.1f}% {b[2]:>10.1f}%")
        out[m_true] = dict(mean_s=float(np.mean(s_v)), mle=list(map(float, a)),
                           moment=list(map(float, b)))
    return out


# ----------------------------------------------------------------- data pass
def run_table(tbl, clone, xi_all):
    d = xi_all[tbl]
    keep = {k: v for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
    per_site = {}
    for site, sd in d["per_site"].items():
        sub = {k: v for k, v in sd["xi"].items() if k in keep}
        tot = sum(sub.values())
        syms = sorted(sub)
        arr = np.array([sub[k] / tot for k in syms])
        per_site[site] = (dict(zip(syms, range(len(syms)))), arr,
                          make_moment_inverter(arr))
    print(f"\n{tbl}: alphabet {len(keep)} symbols", flush=True)

    nodes = defaultdict(Counter)
    with open(DATA / f"{tbl}_EditTable_filtered.csv", newline="") as fh:
        r = csv.reader(fh); hdr = next(r)[1:]
        blocks = defaultdict(list)
        for i, c in enumerate(hdr):
            t, s = c.rsplit(".", 1); blocks[t].append((int(s[4:]), i))
        blocks = {t: [i for _, i in sorted(v)] for t, v in blocks.items()}
        for row in r:
            cl = clone.get(row[0])
            if cl is None: continue
            vals = row[1:]
            for tape, idxs in blocks.items():
                pref = ()
                for i in idxs:
                    v = vals[i]
                    if v not in keep: break
                    nodes[(cl, tape, pref)][v] += 1
                    pref += (v,)

    cache, out = {}, RES / f"m_mle_{tbl}.tsv.gz"
    n_sat = kept = 0
    with gzip.open(out, "wt") as fo:
        fo.write("clone\ttape\tlevel\tn_next\ts\tW_A\tm_moment\tm_mle\trec_moment\trec_mle\n")
        for (cl, tape, pref), cnt in nodes.items():
            n_next = sum(cnt.values())
            if n_next < MIN_CARRIERS: continue
            site = f"Site{len(pref)+1}"
            index, arr, moment = per_site[site]
            idx = np.sort(np.array([index[k] for k in cnt]))
            key = (site, idx.tobytes())
            if key not in cache:
                xa = arr[idx]
                mm, sat = mle_m(xa, float(xa.sum()))
                cache[key] = (float(xa.sum()), mm, sat)
            W_A, m_mle, sat = cache[key]
            n_sat += sat
            s = len(idx)
            m_mom = max(moment(s), s)
            m_mle_c = min(max(m_mle, s), n_next)
            m_mom_c = min(m_mom, n_next)
            kept += 1
            fo.write(f"{cl}\t{tape}\t{len(pref)}\t{n_next}\t{s}\t{W_A:.6f}\t"
                     f"{m_mom_c:.3f}\t{m_mle_c:.3f}\t{m_mom_c-s:.3f}\t{m_mle_c-s:.3f}\n")
    print(f"  nodes {kept:,} | saturated (W_A=1) {n_sat:,} | -> {out.name}", flush=True)
    return kept


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    xi_all = json.loads((RES / "xi_vectors.json").read_text())
    d = xi_all["Mouse1"]
    keep = {k: v for k, v in d["xi"].items() if v * d["n_edits"] >= MIN_COUNT}
    xi_ref = np.array(list(keep.values()))
    val = validate(xi_ref)
    (RES / "mle_validation.json").write_text(json.dumps(val, indent=2))

    if args:
        clone = {}
        with open(DATA / "clonalbc_percell_hamming1_corrected.csv", newline="") as fh:
            for row in csv.DictReader(fh):
                if row["ClonalBC"] and row["ClonalBC"] != "None":
                    clone[row["CellID"]] = row["ClonalBC"]
        for t in args:
            run_table(t, clone, xi_all)
