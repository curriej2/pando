# Analysis: the likelihood engine — exact per-clone likelihood for DNA Typewriter trees

**Question.** Build an exact, differentiable likelihood for one clone's tree under the simulator's
model (editing + heritable silencing + technical dropout + depth-0 censoring), as the core of our own
tree-reconstruction tool, extensible to time-, lineage- and signal-varying editing.

**State (2026-10-05): design settled (decisions 1–7, `notes/sciphy_notes.md` §H.8); proposal parts 1–2 APPROVED 2026-10-06; engine BUILT (`src/01`–`04`), env `envs/pando-jax`
(conda-forge, jax 0.11.2 + CUDA 12). Validation FAILED BY RULE on two CHECK-SET issues: V5 one
finite-difference comparison at a high-curvature boundary point (diagnosed: step-size error, autodiff
correct) and V8 one deterministic library cell; V1–V4 and V7 PASS. Amended V5 + a V8 follow-up
proposed to Justin; level-batched pruning proposed for gradient scaling; V6 (SciPhy Java) not done.
Findings: README.md "Session 22". Write-up: writeup/likelihood_engine.tex.**

Settled design (notes §H.8.6–H.8.9 and the decision records after them):
1. depth-slot state representation (SciPhy prefix chains + missing data; trie later if read errors
   matter); 2. node times stored, ratio-transform coordinates; 3. rate-segment process modules with a
   guardrail, closed-form kernels (matrix exponential only as a background check), parameter classes
   shared with the simulator; 4. lineage-latent hooks only (fixed at "nothing varies"); 5. JAX + a
   NumPy reference twin + SciPhy's Java likelihood as an external base-model check; 6. one clone per
   evaluation, per-tape log-likelihood array, posterior mode on request; 7. one clone at a time, large
   clones first, blank founder, identifiability-by-clone-size study wanted.

---

## PROPOSAL (2026-10-05), PART 1 — THE MATH

### 0. What the engine computes, in one sentence

For one clone, the engine returns, for every tape $z$, the log-probability $\ell_z$ (in nats) of
everything that was read on that tape across the clone's cells, given a tree with node times and the
model's parameters; the clone's log-likelihood is $\ell=\sum_z\ell_z$, and the fitting objective adds
a birth–death prior on the node times and priors on the parameters.

### 1. Notation (grouped; repeated where it is used)

**Tree and time**
| symbol | meaning |
|---|---|
| $t\in[0,1]$ | calendar time as a fraction of the experiment; the clone's founder at $t=0$, harvest at $t=1$ (days $=t\times T$, $T$ known) |
| $h=1-t$ | height: time remaining until harvest |
| $n$ | number of cells (tips) in the clone; the tree is rooted and binary, with $n-1$ internal nodes |
| $o$ | the **origin**: the founder cell at $t=0$, tape blank |
| $r$ | the **root**: the most recent common ancestor of all sampled cells, at time $t_r\in(0,1)$ |
| stem | the branch $o\to r$, spanning $[0,t_r]$ |
| $v,u,w$ | nodes; ${\rm pa}(v)$ = parent, ${\rm ch}(v)$ = children |
| $e$ | a branch, **named by its lower node**: branch $w$ runs from ${\rm pa}(w)$ down to $w$, over $[t^0_w,t^1_w]=[t_{{\rm pa}(w)},t_w]$ |
| $t_v$ | node time; all tips at $t=1$ |
| $x_v$ | the unconstrained coordinate of internal node $v$ (ratio transform, §9) |

**Tapes, cells, reads, hidden states**
| symbol | meaning |
|---|---|
| $k$ | tapes per cell; $z=1..k$ indexes tapes |
| $N$ | slots per tape ($N=6$) |
| $c$ | a cell (a tip) |
| $d_{cz}$ | the **read** of tape $z$ in cell $c$: either an observed string $v_{cz}=(s_1,\dots,s_m)$ with $1\le m\le N$, or "?" (missing). A blank tape never appears as observed — Park's pipeline codes it as missing |
| $X_{vz}$ | the **true (hidden) state** of tape $z$ at node $v$: a string of $\le N$ symbols, or $\mathsf S$ (silenced) |
| $d$ | a **slot** = depth = number of symbols on the tape, $d\in\{0,\dots,N\}$, plus the slot $\mathsf S$ |

**Editing clock**
| symbol | meaning | units |
|---|---|---|
| $\Lambda_T$ | expected edits on an average tape over the whole experiment | edits/tape |
| $\lambda_0(t)$ | time shape of editing, $\ge0$, $\int_0^1\lambda_0=1$, piecewise constant | — |
| $W(t)=\int_0^t\lambda_0$ | share of the experiment's editing spent by time $t$ ($W(0)=0$, $W(1)=1$) | — |
| $r_z$ | speed of tape $z$ relative to the average tape (mean 1) | — |
| $m_w$ | per-branch rate multiplier — **hook for lineage-varying editing; fixed at 1** | — |
| $\lambda_z(t)=r_z\Lambda_T\lambda_0(t)$ | editing rate of tape $z$ at time $t$ | edits/tape/experiment |

**Composition**
| symbol | meaning |
|---|---|
| $\xi_i(t)$ | probability that an edit made at time $t$ writes symbol $i$; $\sum_i\xi_i(t)=1$ |
| $\xi^L$ | the lineage (background) composition — fixed, measured |
| $p_a(t)$, $\xi^{(a)}$ | share and within-channel composition of signal channel $a$ — **hook; no signal channels in the first build**, so $\xi=\xi^L$ |

**Silencing**
| symbol | meaning |
|---|---|
| $\bar H$ | expected fraction of tip × tape entries silenced at harvest |
| $\Gamma_T=-\ln(1-\bar H)$ | total silencing exposure along one root-to-tip lineage (events/tape) |
| $\zeta_0(t)$, $W_\zeta(t)=\int_0^t\zeta_0$ | time shape of silencing ($\int_0^1\zeta_0=1$) |
| $\Gamma(t)=\Gamma_TW_\zeta(t)$ | silencing exposure accumulated by time $t$ |

**Dropout (technical loss at readout)**
| symbol | meaning |
|---|---|
| $\pi_{cz}=\sigma(\alpha_c+\beta_z)$ | probability that tape $z$ of cell $c$ is lost technically; $\sigma(y)=1/(1+e^{-y})$ |
| $\alpha_c$ | cell offset (a poorly captured cell has larger $\alpha_c$) |
| $\beta_z$ | tape offset (a poorly read tape has larger $\beta_z$) |

**Tree prior**: $b$, $\delta$ = division and death rates (per experiment); $\rho$ = fraction of the
final population that was sampled.

**Derived, per branch and segment** (§5): segment $k$ of branch $w$ is $[a_k,b_k]$;
$\mu_{wzk}$ = expected edits on tape $z$ in that segment; $G_{wk}$ = silencing exposure in it;
$\xi^{(k)}$ = composition in it. **Per node and tape** (§4–7): $\mathrm{lcp}_{vz}$, the partial
likelihood $L_{vz}(d)$, the message $M_{w\to v,z}(d)$, the kernel $K_{wz}[d',d]$, the outside
message ${\rm Out}_{vz}(d)$.

### 2. The model the likelihood scores (identical to the simulator's rung 2)

1. **Tree and time.** One founder cell at $t=0$ divides and dies as a birth–death process; at $t=1$
   a fraction $\rho$ of the population is sampled; the tree relating the sampled cells is rooted,
   binary and ultrametric (every tip at $t=1$), with a stem from the founder to the root.
   *[SciPhy's set-up: origin fixed at the known experiment length, stem, birth–death-sampling prior.]*
2. **Editing.** Each tape starts **blank** at the founder. Along every lineage, independently for
   each tape, edits arrive as a Poisson process at rate $\lambda_z(t)$; each edit **appends** one
   symbol, drawn with probabilities $\xi(t)$; once $N$ slots are full, editing stops. At a division
   both daughters inherit the tape.
   *[SciPhy's model, Definitions 1–3, Eqs 3–7; time-varying rate and composition are our §S4.]*
3. **Silencing.** Along every lineage, independently for each tape, silencing strikes at hazard
   $\zeta(t)=\Gamma_T\zeta_0(t)$; a silenced tape stays silenced in every descendant (Dollo).
   *[LAML's absorbing "−1" state; here per tape, in calendar time, with a timing shape — new.]*
4. **Readout.** At harvest each cell's tape is read: a silenced tape is always missing; otherwise it
   is lost technically with probability $\pi_{cz}$; a captured but **blank** tape is also reported
   missing (censoring); otherwise it is read **exactly**.
   *[Two-layer structure from LAML/LAML-Pro; per-(cell,tape) logit dropout and depth-0 censoring —
   new.]*
5. **Independence assumptions.** Given the tree and the parameters: tapes are independent of each
   other (SciPhy; ⚠ rows A5/A8 — cells' tapes share latent factors, so this overstates evidence);
   silencing is independent of editing; technical losses are independent across entries given
   $(\alpha,\beta)$.

### 3. The likelihood is a product over tapes, each a sum over hidden histories

Because tapes are independent given the tree $\mathcal T$ (topology + node times) and the parameters
$\theta$,
$$L(\mathcal T,\theta)=\prod_{z=1}^{k}L_z,\qquad \ell=\log L=\sum_{z=1}^k\ell_z,\qquad \ell_z=\log L_z .$$
For one tape, $L_z$ adds up the probability of every complete story — a true state at every node —
consistent with the reads:
$$L_z=\sum_{\{X_{vz}\}}\underbrace{P_{\rm stem}(\varnothing\to X_{rz})}_{\text{founder blank}}
\prod_{\text{branches }w}\underbrace{P_w(X_{{\rm pa}(w)z}\to X_{wz})}_{\text{what happens on branch }w}
\prod_{\text{cells }c}\underbrace{\mathcal O_{cz}(d_{cz}\mid X_{cz})}_{\text{what was read}} .$$
*Intuition:* each factor answers a local question — given the tape at the top of a branch, how likely
is its state at the bottom; given a cell's true tape, how likely is what the sequencer reported. The
sum is over every way the hidden ancestors could have looked. It is never computed by listing
stories; §7 computes it exactly by passing messages up the tree (Felsenstein pruning).
*[SciPhy Eqs 10–12.]*

### 4. Which states each node can be in: anchored and floating nodes (depth slots)

Write ${\rm Obs}_z(v)$ for the cells below node $v$ whose tape $z$ was read (not missing).
- **Anchored node** (${\rm Obs}_z(v)$ non-empty): the true state must be a prefix of every read tape
  below it, because editing only appends. So it is one of the prefixes of
  $$\mathrm{lcp}_{vz}=\text{longest common prefix of }\{v_{cz}:c\in{\rm Obs}_z(v)\},$$
  and it cannot be $\mathsf S$ (anything readable below means it was never silenced — Dollo).
  Slot $d$ ($0\le d\le|\mathrm{lcp}_{vz}|$) means the specific string $\mathrm{lcp}_{vz}[1..d]$.
- **Floating node** (nothing read below): only the depth matters — a missing read can only say
  whether the tape was blank or silenced, never which symbols it held (verified by exhaustive
  enumeration, notes §H.8.4). Slot $d$ ($0\le d\le N$) means "any string of depth $d$"; slot
  $\mathsf S$ is allowed.

So every (node, tape) carries a vector of $N+2=8$ numbers. Two structural facts used below: if a
child is anchored, so is its parent, and $\mathrm{lcp}_{{\rm pa}(w)z}$ is a prefix of
$\mathrm{lcp}_{wz}$ — so the parent's slot $d'$ and the child's slot $d$ refer to the same string up
to depth $d'$.
*[SciPhy Eqs 8–9 for anchored nodes. The floating-node lumping and the unified 8-slot representation
— new.]*

### 5. Breaking each branch into segments

Each process module (decision 3) reports where it changes: the knots of $\lambda_0$, of $\zeta_0$ and
of the composition. On branch $w$, merge those that fall inside $[t^0_w,t^1_w]$; this cuts the branch
into segments $k=1..K_w$, $[a_k,b_k]$, on each of which the composition is constant. Per segment:
$$\mu_{wzk}=r_z\,\Lambda_T\,m_w\,\big[W(b_k)-W(a_k)\big]\quad(\text{expected edits on tape }z),\qquad
G_{wk}=\Gamma(b_k)-\Gamma(a_k)\quad(\text{silencing exposure}),\qquad \xi^{(k)} .$$
*Intuition:* within a segment where the composition is fixed, only the **total** amount of editing and
of silencing matters, not when inside the segment they happened — speeding the clock up or down plays
the same "film" faster. Segments are needed only where the composition changes, because then which
symbol sits in a slot depends on whether it was written before or after the change. In the first
build (constant composition) every branch is a single segment unless $\lambda_0$ or $\zeta_0$ have
knots — and even then only the integrals matter.
*[Integrated-rate segments: our §S4.4 / §H.8.8. SciPhy uses one constant rate per tape.]*

### 6. The branch kernel — closed form

For branch $w$ (parent $v={\rm pa}(w)$), tape $z$, the kernel $K_{wz}[d',d]$ is the probability that
the tape ends the branch in slot $d$ given it started in slot $d'$ (for an anchored child, *and* wrote
exactly the symbols $\mathrm{lcp}_{wz}[d'+1..d]$ on the way).

**One segment** with editing amount $\mu$, silencing exposure $G$, composition $\xi$:
$$\boxed{\;K[d',d]=\underbrace{e^{-G}}_{\text{not silenced}}\times\underbrace{C_N(d'\to d;\mu)}_{\text{right number of edits}}\times\underbrace{\Xi(d',d)}_{\text{right symbols}}\qquad(d'\le d)\;}$$
- **Right number of edits.** If the tape does not fill ($d<N$), exactly $d-d'$ edits:
  $C_N(d'\to d;\mu)=e^{-\mu}\mu^{\,d-d'}/(d-d')!$. If it fills ($d=N$), at least $N-d'$ edits were
  attempted (extras had nowhere to go): $C_N(d'\to N;\mu)=1-\sum_{j=0}^{N-d'-1}e^{-\mu}\mu^j/j!$. (A full
  tape stays full: $C_N(N\to N)=1$.) *[SciPhy Eq 5.]*
- **Right symbols.** Anchored child: one probability per symbol written,
  $\Xi(d',d)=\prod_{j=d'+1}^{d}\xi_{\mathrm{lcp}_{wz}[j]}$. Floating child: $\Xi=1$ (whatever was written
  is summed over, and each write's probabilities sum to 1). *[SciPhy Eq 6; the floating case — new.]*
- **Not silenced.** $e^{-G}$; silencing is independent of what is on the tape.
  *[LAML's $\Psi$ factorises the same way — verified: LAML's per-site kernel is exactly this at $N=1$.]*
- **Silenced slot** (floating child only): $K[d',\mathsf S]=1-e^{-G}$, $K[\mathsf S,\mathsf S]=1$,
  $K[\mathsf S,d]=0$. For an anchored child the $\mathsf S$ column is not used (Dollo).
- $K[d',d]=0$ for $d<d'$: editing never removes symbols.

**Several segments:** play them in time order,
$$K_{wz}=K^{(1)}K^{(2)}\cdots K^{(K_w)},\qquad\text{i.e.}\quad (AB)[d',d]=\sum_{j}A[d',j]\,B[j,d],$$
where each $K^{(k)}$ is the boxed formula with that segment's $(\mu_{wzk},G_{wk},\xi^{(k)})$, and $j$
runs over the depth the tape had reached at the segment boundary. *Intuition:* the sum enumerates
every way the edits could be split between the segments, each symbol scored by the composition in
force when it was written. *[§S4.4's epoch recursion; checked against direct convolution, §H.8.8.]*

**Written out in full** (added 2026-10-06). With $j_0=d'$, $j_K=d$ and $j_1\le\dots\le j_{K-1}$ the depths
at the $K-1$ segment boundaries,
$$K_{wz}[d',d]=e^{-\sum_kG_{wk}}\sum_{d'\le j_1\le\cdots\le j_{K-1}\le d}\ \prod_{k=1}^{K}C_N\big(j_{k-1}\to j_k;\mu_{wzk}\big)\,\Xi^{(k)}(j_{k-1},j_k),
\qquad K_{wz}[d',\mathsf S]=1-e^{-\sum_kG_{wk}} .$$
**When the composition is the same in every segment** (the first build), the sum collapses to ONE
segment with $\mu=\sum_k\mu_{wzk}$ and $G=\sum_kG_{wk}$ — Poisson counts add, saturation included
(✅ checked 2026-10-06: 3 segments, 200 draws, max $|\Delta|=2.2\times10^{-16}$). Segments then only ever
arise where the composition changes.

**Segment boundaries are not inferred.** They are the knots of $\lambda_0$, $\zeta_0$ and the
composition — fixed inputs of the model (the knot *levels* may be free parameters; their *positions*
are set by the analyst or the experimental design). What is inferred is the node times, and the
kernel is a continuous, differentiable function of them: $\mu$ and $G$ are integrals of piecewise-constant
rates up to the node time, and a node moving across a composition knot just shrinks one segment to
zero length (a zero-length segment is the identity). So gradients exist and no Monte Carlo over
segments is needed. Inferring an unknown switch *time* later would be a latent (decision 4).

### 7. Tips, pruning, and the root

**Tip vectors** (cell $c$, tape $z$):
$$\text{read as }v_{cz}:\quad L_{cz}(d)=(1-\pi_{cz})\,\mathbf 1\big[d=|v_{cz}|\big],\ \ L_{cz}(\mathsf S)=0;\qquad
\text{missing:}\quad L_{cz}(d)=\pi_{cz}+(1-\pi_{cz})\mathbf 1[d=0],\ \ L_{cz}(\mathsf S)=1 .$$
*Intuition:* an observed tape was captured (probability $1-\pi$) and is exactly what we read. A
missing tape was lost technically (probability $\pi$, whatever it held), or captured but blank
(probability $1-\pi$, only if $d=0$), or silenced (probability 1).
*[LAML-Pro's emission framework; per-(cell,tape) dropout and the blank-tape term — new.]*

**Message from child $w$ up to its parent** (for each parent slot $d'$):
$$M_{w\to{\rm pa}(w),z}(d')=\sum_{d}K_{wz}[d',d]\,L_{wz}(d)$$
— "if the parent was in slot $d'$, how probable is everything in $w$'s subtree?"

**Partial likelihood at an internal node** $v$ (the subtrees below $v$ evolve independently given $v$):
$$L_{vz}(d')=\prod_{w\in{\rm ch}(v)}M_{w\to v,z}(d')\quad\text{on }v\text{'s valid slots, and 0 elsewhere.}$$

**Root and stem.** The founder is blank, so the tape enters the stem in slot 0:
$$L_z=\sum_{d}K_{{\rm stem},z}[0,d]\,L_{rz}(d),\qquad \ell_z=\log L_z .$$
(If tape $z$ was missing in every cell, the root is floating and $d$ includes $\mathsf S$.)
Numerically the recursion is carried in log space or with per-node rescaling, since products over
thousands of nodes underflow. *[SciPhy Eq 12 / Felsenstein pruning.]*

### 8. The per-tape array, and tape subsets for support

The engine returns the vector $(\ell_1,\dots,\ell_k)$ in nats; the clone's log-likelihood is
$\ell=\sum_z\ell_z$. Any re-weighting $\ell_{\mathbf w}=\sum_zw_z\ell_z$ is then free:
- tape bootstrap: $\mathbf w$ = counts of a with-replacement draw of $k$ tapes;
- jackknife: $w_z\in\{0,1\}$, one tape left out;
- weighted likelihood bootstrap: $\mathbf w=k\times\mathrm{Dirichlet}(1,\dots,1)$, refit, repeat.
A tape whose $\ell_z$ is persistently out of line flags a problem tape.
*[The tape as the resampling unit: Park's and the Typewriter paper's tape bootstraps; §F.2.
Returning the per-tape likelihood array as a first-class output — new to these methods.]*

### 9. Coordinates, priors, and the objective

**Times from coordinates (ratio transform).** Heights from the root down:
$$h_r=\sigma(x_r),\qquad h_v=\sigma(x_v)\,h_{{\rm pa}(v)}\ \ (v\ne r\text{ internal}),\qquad t_v=1-h_v .$$
Each node sits a fraction $\sigma(x_v)$ of the way from the present up to its parent, so **every real
vector $x$ is a valid ultrametric tree** with the root below the founder. Its log-Jacobian is
$$\log\Big|\frac{\partial h}{\partial x}\Big|=\sum_{v\ \rm internal}\Big[\log h_{{\rm pa}(v)}+\log\sigma(x_v)+\log\big(1-\sigma(x_v)\big)\Big],\qquad h_{{\rm pa}(r)}=1 .$$
(Both properties checked on 200 random trees, notes §H.8.7.) *[Kishino–Thorne–Bruno node-height
ratios — citation from memory, to verify.]*

**Tree prior** (constant-rate birth–death with sampling fraction $\rho$; heights in experiment units;
origin at height 1):
$$p(\mathcal T\mid b,\delta,\rho)\ \propto\ p_1(1)\prod_{v\ \rm internal}b\,p_1(h_v),\qquad
p_1(h)=\frac{\rho\,(b-\delta)^2e^{-(b-\delta)h}}{\big[\rho b+\big(b(1-\rho)-\delta\big)e^{-(b-\delta)h}\big]^2},$$
conditioned on the clone having its observed $n$ cells. $p_1(h)$ is the probability that a lineage
alive at height $h$ leaves exactly one sampled descendant; each internal node contributes "a division
happened here, and each daughter left sampled descendants". This $p_1$ is §S3's (checked there by
finite differences); ⚠ the normalisation for conditioning on $n$, and agreement with the tree library,
are to be verified in part 2. *[SciPhy's / BEAST's birth–death-sampling prior (Stadler). It replaces
ConvexML's minimum branch length and pseudocounts as the regulariser — ConvexML's lesson.]*

**The objective.** With $\theta$ the free parameters of the blocks chosen (decision 7):
$$J(x,\theta)=\underbrace{\sum_z\ell_z\big(t(x),\theta\big)}_{\text{data}}+\underbrace{\log p\big(t(x)\mid b,\delta,\rho\big)}_{\text{tree prior}}+\underbrace{\log p(\theta)}_{\text{parameter priors}}\ \Big[+\log\Big|\tfrac{\partial h}{\partial x}\Big|\ \text{when sampling}\Big].$$
- **Point estimate (penalised ML / MAP in time):** maximise $J$ *without* the Jacobian; $x$ is only
  the optimiser's coordinate system, so the optimum is the MAP in time coordinates.
- **Posterior sampling (HMC) in $x$:** include the Jacobian, so the samples follow the posterior over
  times.
Gradients of $J$ with respect to $x$ and $\theta$ come from automatic differentiation (JAX).

### 10. Posterior mode — ⏸ DEFERRED to a later version (Justin, 2026-10-06)

Not part of the first build. Kept here as the specification for when it is added (it is the bridge to
signal inference and to per-entry missingness attribution).

**Outside messages** — "how probable is everything *outside* $v$'s subtree, if $v$ is in slot $d$" —
computed top-down:
$${\rm Out}_{rz}(d)=K_{{\rm stem},z}[0,d],\qquad
{\rm Out}_{wz}(d)=\sum_{d'}\Big[{\rm Out}_{vz}(d')\prod_{u\in{\rm ch}(v),\,u\ne w}M_{u\to v,z}(d')\Big]K_{wz}[d',d],\quad v={\rm pa}(w).$$
Then, with $L_z$ the tape's likelihood:
- **Node state:** $P(X_{vz}=d\mid\text{reads})={\rm Out}_{vz}(d)\,L_{vz}(d)/L_z$ — the ancestor's depth,
  and at anchored nodes its exact prefix $\mathrm{lcp}_{vz}[1..d]$.
- **Edits on a branch:** $P(X_{vz}=d',X_{wz}=d\mid\text{reads})=\big[{\rm Out}_{vz}(d')\prod_{u\ne w}M_{u\to v,z}(d')\big]K_{wz}[d',d]\,L_{wz}(d)/L_z$
  — the probability that branch $w$ took the tape from depth $d'$ to $d$, i.e. wrote
  $\mathrm{lcp}_{wz}[d'+1..d]$ (anchored child). This is the input to signal inference.
- **Why a missing entry is missing** (cell $c$, tape $z$, with tip posterior $q(d)=P(X_{cz}=d\mid\text{reads})$):
  $$P(\text{silenced})=q(\mathsf S),\qquad P(\text{blank, censored})=q(0)\,(1-\pi_{cz}),\qquad
  P(\text{technical loss})=\sum_{d=1}^{N}q(d)+q(0)\,\pi_{cz},$$
  which sum to 1. *[LAML's per-entry imputation, extended to three causes.]*
(Equivalently, these posteriors are gradients of $\ell_z$ with respect to the tip emissions and branch
kernels — a standard identity, to be checked against brute force in part 2.)

### 11. Cost

Per tape: building lcps is $O(nN)$ (topology only, cached until the topology changes); kernels and
pruning are $O\big((2n-1)\,\bar K\,(N+2)^2\big)$ with $\bar K$ the mean number of segments per branch.
Tapes are independent ⇒ computed in parallel. *[Same order as SciPhy's $O(knN^2)$.]*

### 12. Deliberately not in this version

Read/sequencing error (would need the trie representation); cross-tape coupling (rows A5/A8, the
*trans* limb); lineage-varying rate or composition (hooks fixed); partial silencing; per-division
editing coupling; editing before the founder ($\Lambda_{\rm pre}>0$); pooling across clones; the
site-6 mechanism (§S3 measured slot 6 behaves differently on Park).

### 13. Provenance

| component | taken from / similar to | what is new here |
|---|---|---|
| tape = character; prefix states; kernel Eqs 5–7; ancestral sets Eqs 8–9; pruning Eq 12 | **SciPhy** | — |
| time tree with stem, fixed experiment length, birth–death-sampling prior, per-tape clocks | **SciPhy / BEAST** | — |
| absorbing heritable-silencing state; editing + observation as two layers | **LAML** (its kernel = ours at $N=1$, verified) | silencing per tape, in calendar time, timing shape $\zeta_0$, parameterised by $\bar H$ |
| general per-leaf emission; per-entry posterior attribution | **LAML-Pro / LAML** | per-(cell, tape) logit dropout; **depth-0 censoring**; three-way attribution |
| regularise branch lengths | **ConvexML**'s lesson | done with the tree prior instead of ε / pseudocounts |
| node-height ratio coordinates | Kishino–Thorne–Bruno; BEAST HMC work | — |
| **depth-slot representation** unifying anchored and floating nodes | extends SciPhy Eqs 8–9 | **new** (verified exact) |
| time-varying rate *and* composition via integrated-rate segments | our §S4 | **new** relative to all four methods |
| rate-segment module interface with guardrail; lineage hooks | — | **new** |
| per-tape log-likelihood array as an output; edit-to-branch posterior as the bridge to signal inference | tape bootstrap (Park, Typewriter) | **new** as engine outputs |

---

## PROPOSAL (2026-10-06), PART 2 — BUILD AND VALIDATION. ✅ APPROVED 2026-10-06, BUILT; VALIDATION FAILED BY RULE (check-set issues; README Session 22)

### 1. The question

Does the engine compute exactly the likelihood of part 1 — the probability, under the simulator's
rung-2 model, of the reads on each tape of one clone — with correct gradients, and how fast is it as a
function of clone size? **A failed check is a bug, not a finding.** What depends on it: every later
fit (tree search, support values, the identifiability-by-clone-size study, signal inference) sits on
this function. **No scientific run is included**; the identifiability study returns as its own
proposal, because its numbers will be quoted.

### 2. Deliverables

| step | file | what it is |
|---|---|---|
| 0 | environment | **On Justin's go only:** a separate env `/data1/choij10/justin/envs/pando-jax` (python 3.12, `jax[cuda12]`, numpy, scipy, matplotlib), plus a one-minute GPU smoke job on `-p lesliec,gpu --gres=gpu:1`. The simulator's env is untouched. |
| 1 | `src/01_engine_ref.py` | **NumPy reference twin** — plain loops, written for readability: anchors/lcps, closed-form kernels, pruning, per-tape $\ell_z$, ratio transform, tree prior. Slow and obviously correct. |
| 2 | `src/02_engine.py` | **JAX engine** — the same functions on fixed-shape arrays: topology passed as data (parent index + level order), tapes batched, log-space, float64 throughout, `jit` + `grad`. Returns the per-tape array $(\ell_1..\ell_k)$ and the objective $J$. |
| 3 | `src/03_validate_engine.py` | the validation suite below, one JSON + one printed table |
| 4 | `src/04_time_engine.py` | timing vs clone size (descriptive only) |

Inputs are the simulator's own arrays (true depth, symbols, technical mask, silencing) and parameter
objects (`14_editing.Recorder`, `17_dropout` panel/mask, `22_silencing.SilencingParams`), **imported
as-is — no simulator code is changed** in this build (moving the classes to `src/` stays a later,
separate, bit-identical step). Posterior mode, tree search and optimisers are NOT in this build.

### 3. The checks (pre-registered)

**Hard checks** — exact up to floating-point; tolerance stated; any violation is a defect.
- **V1 exhaustive enumeration.** On small trees ($N=3$, $M=3$, 3–6 cells, all 41 states incl. $\mathsf S$),
  sum the probability of every hidden history directly and compare with the NumPy twin's $\ell_z$.
  Random parameters; observation patterns incl. all-missing subtrees, nothing observed, divergence at
  every slot, homoplasy (same symbol written twice independently), multi-segment branches with a
  composition change. Tolerance: relative error $\le10^{-10}$. *Tests decision 1's representation and
  every formula of part 1 against first principles.*
- **V2 normalisation.** On the same tiny settings, enumerate every possible set of reads (every
  observed string or "?" per cell) and sum $L$ over them: must equal 1 to $10^{-10}$. *A likelihood
  that is a proper probability distribution over what can be observed — tests the emission and
  censoring terms, which V1 alone cannot (V1 checks one dataset at a time).*
- **V3 twins agree.** JAX vs NumPy on simulator data: library trees $n\in\{25,64,210\}$ × 3 capture
  fractions × 2 turnovers × 5 trees, $k=30$, $N=6$, Park $\xi^L$, rung-2 settings (technical dropout,
  censoring, silencing at $\bar H\in\{0,0.1\}$, constant and early timing). Per-tape $\ell_z$ relative
  error $\le10^{-10}$. *Tests the vectorised, padded, topology-as-data implementation.*
- **V4 kernel cross-check.** Closed form vs matrix exponential on random segments and multi-segment
  branches: $\le10^{-12}$ (the background check of decision 3).
- **V5 gradients.** JAX gradients vs central finite differences, w.r.t. every node coordinate $x_v$
  and each parameter ($\Lambda_T$, $\bar H$, $\mu_\beta$, $\alpha_c$, $b$, $\delta$): relative
  $\le10^{-5}$ on 20 trees.
- **V6 SciPhy cross-check — optional, needs its own setup (BEAST 2.7 + SciPhy package + the `java`
  module).** Base model only (no dropout, silencing or censoring; complete tapes): our $\ell_z$ vs
  SciPhy's tree likelihood on the same tree and parameters, $\le10^{-8}$. *The only check written by
  someone else — it catches a misconception the twin and the simulator might share.*

**Statistical checks** — deviation in standard errors across trees (as in `15`/`18`/`23`).
- **V7 score at the truth.** For data simulated at parameters $\theta^*$ on the true tree, the expected
  gradient of the log-likelihood at $\theta^*$ is zero (a property of any correct likelihood). Mean of
  $\partial\ell/\partial\theta$ at $\theta^*$ over 200 simulated datasets, per parameter and per node-time
  block, against 0. *Tests that the engine and the simulator mean the same model — e.g. the same
  normalisation of $\lambda_0$, the same censoring rule.*
- **V8 tree prior.** The owed check of part 1 §9: the birth–death density, conditioned on $n$, against
  the tree library — the distribution of root height and of the internal heights (deciles) implied by
  the density vs the library's realised trees, at the library's $(b,\delta,\rho)$ grid.

### 4. What could make it wrong

**Defeated.** No Park data enters, so the edit-data circularity cannot reach this build. V1/V2 are
first-principles and share no code with the engine; V3/V4 are two implementations each.
**NOT defeated.**
- ⚠ **A misconception shared by the engine and the simulator** (both written by me, with shared
  parameter objects) would pass V3 and V7. V1/V2 reduce this (first principles); only V6 (SciPhy's
  code) removes it, and only for the base model.
- ⚠ **Passing means the engine computes OUR model**, not that our model fits Park: rows A5/A8, partial
  loss, the *trans* limb, slot 6 and read errors are all outside it (part 1 §12).
- ⚠ **Numerical:** float32 underflows on large trees; the build is float64 throughout (on GPU this
  costs speed — measured in step 4).

**Abandon rather than patch if:** V1 fails on floating (all-missing) subtrees while passing on fully
observed ones — that would mean the depth-slot lumping (decision 1) is wrong, and the answer is to
revisit the representation (trie), not to patch the code. If V7 is biased while V1–V5 pass: the engine
and the simulator disagree about the model — find which convention is wrong; never tune the engine to
match the simulator.

### 5. Output — the mock, and the reading rule

`results/validate_engine.json`, printed as one table. ⚠ **Placeholder numbers.**

> *Reading rule.* Hard checks (V1–V6): the number is the worst relative error; PASS iff it is below the
> stated tolerance in every configuration. Statistical checks (V7–V8): $|z|$ against the no-effect
> value 0, se across trees; PASS iff worst $|z|\le3.0$; one trip ⇒ re-run at a new seed, two trips or
> a sign pattern ⇒ defect. **Overall PASS = all hard checks pass AND no statistical defect.**

```
  check                                   compared against                    configs   worst      tol    verdict
  V1  NumPy twin vs exhaustive sum        all hidden histories, 41 states        400     3e-15    1e-10   PASS
  V2  sum of L over all possible reads    1                                      120     1e-14    1e-10   PASS
  V3  JAX vs NumPy, per-tape l_z          library trees n 25-210, rung 2         900     4e-14    1e-10   PASS
  V4  closed form vs matrix exponential   random segments                       1000     5e-15    1e-12   PASS
  V5  autodiff vs finite differences      x_v, Lam_T, Hbar, mu_beta, alpha, b, d  20     3e-07    1e-05   PASS
  V6  ours vs SciPhy (base model)         BEAST 2.7 + SciPhy                      50     2e-10    1e-08   PASS / SKIPPED
  V7  mean score at the truth             0                          per param   12     |z| 2.1   3.0    PASS
  V8  tree prior vs library               root height, height deciles  per cell  60     |z| 2.4   3.0    PASS
```

**Timing (descriptive, not in the verdict)** — seconds per evaluation of all $k$ tapes, and per
gradient, CPU (1 core) and GPU (1× A100), float64, at library sizes $n$ = 25, 210, 1,976, 11,081
(Park's largest clone is 10,997). Reading: the cost of one likelihood / gradient at each size; it
decides whether tree search needs incremental re-evaluation (decision 5, option E).

### 6. Cost and where it runs

V1/V2 are exhaustive and small (minutes, 1 core); V3–V8 minutes to an hour on 1–4 cores; timing one
short GPU job. All via `scripts/submit.sh` on `-p cpu`; GPU only for step 0's smoke job and step 4,
on `-p lesliec,gpu --gres=gpu:1`, small and short. Results JSON committed; nothing large written.

### 7. Decisions for Justin

| decision | proposed | alternative |
|---|---|---|
| JAX environment | separate `pando-jax` with `jax[cuda12]`, created on your go | — |
| simulator code | import the existing parameter classes; change nothing now | graduate to `src/` first (separate bit-identical step) |
| precision | float64 everywhere | float32 on GPU for speed, after V3 shows where it is safe |
| V6 (SciPhy Java) | include, as a separate setup step (download BEAST 2.7 + SciPhy) | defer; rely on V1/V2 for first principles |
| identifiability-by-clone-size study | its own proposal after this PASSes | — |

