# simulator — findings

**State 2026-09-20: the tree layer is BUILT and VALIDATED at both capture fractions. The cost
pilot (`03`) has not been run.** Editing and dropout are separate later layers, not in scope here.

## ✅ Validation result (`02`, 2026-09-20)

**PASS at both $\rho = 0.5$ and $\rho = 8\times10^{-4}$**, 4,000 accepted trees per $n \in \{4,5,8\}$,
turnover 0.3, against closed forms verified independently in §S3.3. **Zero structure faults** in
24,000 trees — every tree had exactly $n-1$ branching times, $2n-1$ nodes, two children per internal
node and all times inside $(0,T)$. The worst statistic in either run sat at **0.64** (`rho50`) and
**0.71** (`rho8e4`) of its own p99 critical value. Cost: 37 s and 25 min, 0.07 / 0.14 GB peak RSS.

⭐ **The two runs really did exercise different regimes.** Nodes traversed per branch point — how
much work degree-2 suppression had to do — was **3.6–3.8 at $\rho=0.5$ against 22.0–25.4 at Park's
$\rho$**, a 6.5× difference. The deep-pruning path that Park adjudication depends on is now tested,
not assumed.

See `CLAUDE.md` in this directory for the full design, and `notes/sciphy_notes.md` §S3 for the
theory record.

## Scope change, recorded (2026-09-20)

**Forward simulation is now the primary route, not the validator.** §S3.3.5 listed direct
conditioned sampling from the BDS density as "the one to build", with forward+reject as validation
only. That ordering assumed the tree generator serves adjudication against Park. On the stated
priority — **purposes I and IV, with IV leading** — it is the other way round: the design sweep has
$n$ and $\rho$ as *design variables*, not quantities to be matched, so it needs no conditioning, no
rejection and no density. The $n/\rho$ blow-up that makes forward simulation expensive is entirely
an artefact of Park's $\rho \approx 8\times10^{-4}$; a design study at $\rho \sim 0.1$–$1$ is cheap.

⇒ build forward first, **measure** what it costs, and reach for the density only where the
measurement says it is needed — which is the sequence `03_cost_curve.py` exists to decide.

## The three scripts

| file | what it does |
|---|---|
| `src/01_bdtree.py` | the simulator. Gillespie loop over the population size → genealogy → $\rho$-sampling → reconstruction. Two implementations: `simulate_tree_naive` (one event at a time, nothing pruned until the end — the oracle) and `simulate_tree` (production) |
| `src/02_validate_tree.py` | five correctness checks against closed forms already verified in §S3.3, reported in **se** not p-values: count law, structure, root split, labelled histories, naive-vs-fast |
| `src/03_cost_curve.py` | the cost measurement and the per-arm projection. Decides whether the BDS density is needed at all |

## Two design points worth keeping

**⭐ The two-pass split, and why it is exact.** At each event the lineage it lands on is chosen
uniformly and independently of everything else, so the law factorises as (count trajectory) ×
(uniform lineage assignments | count trajectory), and the acceptance test reads only $N(T)$ and a
$\mathrm{Binomial}(N(T),\rho)$ draw — never the genealogy. So pass 1 can get $N(T)$ with **no
genealogy at all**, vectorised in numpy, and the $\approx 2.7n$ rejected attempts never allocate a
node; pass 2 replays the identical trajectory from the same seed and builds the genealogy only on
acceptance. Given $\mathrm{Binomial}(N,\rho)=k$ the captured set is a uniform $k$-subset, which is
how pass 2 draws it.

**⭐ The cost study has exactly two free numbers.** $\alpha$ and $\beta$ depend on $(b,\delta,T)$
only through $u=e^{(b-\delta)T}$ and the turnover $\theta=\delta/b$, so at fixed $\rho$ the whole
rejection cost is a function of $(n,\theta)$ and $T$ is a pure scale choice, fixed at 1
(`rate_params`). ⚠ That is **not** true for adjudication against Park, where $T$ is pinned by
$\Lambda=\lambda T$ at the measured $\hat\Lambda$ — it is true for the cost study only.

## ⚠ Corrections to my own proposal of 2026-09-20, before anything ran

1. **Online extinct-pruning was proposed as the memory mitigation; it is not needed.** Working the
   arithmetic: total nodes $\approx 2N(T)/(1-\theta)$, so online pruning saves only ~2× at
   $\theta=0.3$ and ~4× at $\theta=0.75$ — the same order, not the difference between feasible and
   OOM. At Subclone scale the simple "keep everything, prune at the end" path is a few GB.
   ⇒ built the simple version, and `03` records peak RSS so the question is settled by measurement.
2. **Memory is not the binding constraint; wall-clock is.** The genealogy loop (the only sequential
   part) is paid once per acceptance and is a fraction of a percent of the total, so it is
   deliberately left unoptimised.

## What is owed

1. Run `02`. It decides whether the simulator is correct; `03`'s numbers mean nothing until it passes.
2. Run the `03` pilot at small $n$, read $\gamma$, $\kappa$, $\tau$ and `MaxRSS` off a **completed**
   job, then size the large-$n$ points from measurement (CLAUDE.md, "Right-size the request").
3. $R$ — replicates per clone — is deliberately unset. It multiplies the whole cost linearly and is
   set by what the rung 0/1/2 comparison needs. Justin's call, deferred until the timing is known.
4. The forward-vs-BDS equivalence check, **only if the cost measurement says the density is needed**.
   `01_bdtree.py` already carries $h$ and $G$ (fenced, unused by the simulator) and `02` persists
   branching times, so that check needs no re-simulation.
