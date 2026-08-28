# SciPhy Study Notes

> **ORIENTATION FOR A FRESH SESSION — read this block first.**
>
> Working notes from a slow, equation-by-equation read of the SciPhy paper, plus surrounding
> literature, plus design notes for a planned extension. The reader is **Justin, a student at
> Memorial Sloan Kettering Cancer Center.**
>
> **The project:** extend the SciPhy framework to infer **signalling history** from
> ENGRAM + DNA Typewriter data; secondary goal of making the inference scale.
>
> **Preferred working style (established at the outset):** derive step by step rather than
> summarize; unpack equations term by term; maintain this document as the running record.
>
> **Already settled here — do not re-derive unless asked:**
> - DNA Typewriter mechanism at the sequence level (§0), incl. key-vs-PAM and the *cis*/*trans*
>   decomposition of rate variation.
> - SciPhy's editing model, transition probabilities, likelihood (§1a–1c), incl. the
>   prefix-set / longest-common-prefix structure of Eqs. 8–9 and Felsenstein pruning.
> - Inference background: MCMC, tree space, MAP vs ML, HMC, OU fields, PCM (§2a, §D, §G).
> - Architecture decision: **point-estimate topology + full posterior on continuous parameters** —
>   defensible, and already published in this literature (§F.3, §H.2).
> - **The A1–A9 table** of SciPhy assumptions vs. molecular reality — the master list of extension
>   points, referenced throughout as "row A_n_."
>
> - **Mulberry & Stadler 2026, full close read (§H.6)** — every equation, numerically verified.
>   Three findings not in the paper: the ignored $\mathrm{Cov}(X_\text{in},X_\text{out})$, the
>   CLT-for-a-tail error, and the suboptimal multi-rate tape allocation.
>
> **Live threads, priority order (revised after the Mulberry close read):**
> 1. ⚑⚑ **Read Chen et al. 2024 ENGRAM in full — plus the 2025 *Nature Protocols* "Multichannel
>    genomic recording with ENGRAM."** **This is now the immediate next task** and it displaces
>    Schiffman. Every design conclusion from §H.6 is blocked on molecular facts only ENGRAM
>    supplies. Confirmed already: **ENGRAM as published writes all channels to a *shared* DNA
>    Tape** ⇒ the $q$ problem is unmitigated in the published architecture (§I.6).
> 2. Cheap and decisive, run in parallel: cross-tape character-compatibility fraction on Park et
>    al.'s data (§D.4) — one afternoon, decides whether the perfect-phylogeny route is live.
> 3. **NEW — the corrected reconstruction bound** (§H.6.14): keep the covariance, use exact
>    convolution, add a tape-frailty term, optimise the multi-rate allocation jointly. Now a
>    well-specified, self-contained methods contribution.
> 4. Then **Schiffman (PATH)** and **Zwaans (GABI)**; then **LAML** and **ConvexML**. *(Downstream
>    of the tree — not blocking the design.)*
> 5. Verify the $H_m$ quadrature recursion for time-varying $f$ (§1c.7) numerically.
> 6. Identifiability: can $\xi(t)$ be recovered from $N=6$ ordered positions?
> 7. Experimental design (§I) — three tiers; cancer/HPCS application most developed.
>
> **Also worth adding to the list:** *Liao et al. 2024* (Mulberry's citation for "$k$ is hard to
> increase in practice" — $k$ is the biggest lever, so this is the paper that says why we can't
> pull it); *Loveless et al. 2025, PeChyron* (the other sequential recorder Mulberry models).
>
> **Open lead not yet worked out:** §I.5 (Rudensky / Treg / TGF-β).

Running notes for a slow read of **Seidel et al., "SciPhy: A Bayesian phylogenetic
framework using sequential genetic lineage tracing data"** (*Nat Commun* 2026,
doi:10.1038/s41467-026-73377-6).

**Papers covered in this document:**

*Recorder technology*
- **Choi et al. 2022**, "A time-resolved, multi-symbol molecular recorder via sequential genome
  editing" (*Nature* 608:98–107) — DNA Typewriter.
- **Chen et al. 2024**, ENGRAM (*Nature* 632:1073–1081) — symbolic recording of signalling and
  cis-regulatory activity to DNA. *(cited throughout; not yet read in full — worth doing early)*

*Applications of sequential recorders*
- **Yang, Kim, Seidel et al. 2026**, mouse embryo lineage tracing (bioRxiv 2026.07.29.741625) —
  1.34M-cell E13.5 phylogeny.
- **Park, Chang et al. 2026**, metastasis lineage recording (bioRxiv 2026.08.10.744013) —
  NCI-H1299 xenograft, 166 tapes × 6 sites, 105k cells.

*Inference methods*
- **Pilarski, Stadler & Seidel 2026** (*PLOS Comput Biol* 22:e1014370) — simulation benchmark of
  TiDeTree and SciPhy.
- **Mulberry & Stadler 2026** (*Theor Popul Biol* 168:32–43) — analytic reconstruction bounds for
  sequential recorders.
- **Zwaans, Seidel, Manceau & Stadler 2025** (*Phil Trans R Soc B* 380:20230318) — GABI.
- **Schiffman et al. 2024** (*Nat Genet* 56:2174–2184) — PATH / PATHpro.

*Biology for planned applications*
- **Hamazaki et al. 2024** (*Nat Cell Biol* 26:1790) — RA-gastruloids.
- **Marjanovic et al. 2020** (*Cancer Cell* 38:229) — HPCS discovery.
- **Chan et al. 2026** (*Nature* 651:231) — HPCS lineage tracing and ablation.

**Goal:** understand SciPhy well enough to (a) extend it to inference of *signalling
history* with ENGRAM, and (b) evaluate a move to simulation-based inference for
scalability.

---

## Session 0 — Background: how DNA Typewriter actually works

### 0.1 The division of labor: prime editor = machine, pegRNA = symbol

The **prime editor** is a single generic protein — Cas9 **H840A nickase** fused to an
M-MLV **reverse transcriptase**. PE2 in the original Typewriter work; **PEmax +
hMLH1dn** (dominant-negative MLH1 to suppress mismatch repair reversion) in later work
and in the mouse embryo. On its own it has no address and no message.

All specificity and all *content* come from the **pegRNA**:

| Part of pegRNA | Function |
|---|---|
| 20-nt **spacer** (5′) | base-pairs the protospacer → specifies *where* Cas9n nicks |
| scaffold | binds Cas9 |
| 3′ ext: **PBS** (primer binding site) | anneals the nicked 3′ DNA end |
| 3′ ext: **RTT** (RT template) | encodes *what sequence gets written* |

**Write cycle** (Typewriter Fig. 1b):

1. Cas9n nicks the **PAM-containing strand**, 3 bp 5′ of the PAM.
2. The liberated 3′ DNA end anneals to the pegRNA **PBS**.
3. RT extends that 3′ end, reverse-transcribing the **RTT** into the genome as a
   **3′ flap**.
4. Flap equilibration; 5′ flap excised (FEN1/EXO1); ligation.
5. Mismatch repair resolves the heteroduplex. (hMLH1dn biases this toward *keeping*
   the edit — hence the ~3-fold boost from PEmax+hMLH1dn.)

**Key consequence:** the information written lives in the *pegRNA*, not in the target.
Every DNA Typewriter pegRNA shares the **same 20-nt spacer** and differs only in its
RTT. This is why thousands of concurrent symbols can be written to a **single** target
site — the symbol is a property of the writer, not of the address. It is also why
"epistasis" between symbols is minimal and why symbol identity and target position are
orthogonal degrees of freedom.

Also note: **nickase, not nuclease.** No double-strand break ⇒ no NHEJ, no inter-target
deletion, no wholesale excision of consecutive target sites. This is the single biggest
molecular advantage over GESTALT/scGESTALT/CARLIN-style recorders, and it is why the
"information loss" failure mode is largely absent *within* a tape.

### 0.2 Why sequential editing works — at the level of actual bases

**TAPE-1 architecture:** a 3-bp **key** `GGA`, followed by a tandem array of a 14-bp
**monomer** `TGATGGTGAGCACG` whose positions 4–6 are `TGG` — i.e. **each monomer carries
a PAM**.

```
GGA TGATGGTGAGCACG TGATGGTGAGCACG TGATGGTGAGCACG ...
key      monomer1        monomer2        monomer3
```

A functional target = 20 bp protospacer immediately 5′ of a PAM.

**At monomer2's PAM** (the only intact site):

```
  GGA   +  TGATGGTGAGCACG  +  TGA        PAM = TGG (m2 positions 4-6)
  key         monomer1        m2[1-3]
  = GGATGATGGTGAGCACGTGA   (20 bp)   ✓ ACTIVE
```

**At monomer3's PAM** (and every PAM further downstream):

```
  ACG   +  TGATGGTGAGCACG  +  TGA        PAM = TGG (m3 positions 4-6)
 m1[12-14]    monomer2        m3[1-3]
  = ACGTGATGGTGAGCACGTGA   (20 bp)   ✗ INACTIVE
      ^^^ 3 mismatches vs. spacer, PAM-distal
```

The key `GGA` has been replaced by the tail of the preceding monomer, `ACG`. **This is
what "partial CRISPR–Cas9 target sites, all but the first truncated at their 5′ ends and
therefore inactive" means.** A 3-bp PAM-distal mismatch is empirically sufficient to
block prime editing: only **0.6%** of 2×TAPE-1 reads had site 2 edited without site 1
(Fig. 1d).

**The write, and the carriage return.** Cas9n nicks 3 bp from the PAM, which lands
*exactly at the monomer1/monomer2 junction*. The 5-bp insertion is `NN` (2-bp symbol
barcode, unique per pegRNA) + `GGA` (the constant key):

```
GGA TGATGGTGAGCACG NN GGA TGATGGTGAGCACG TGATGGTGAGCACG ...
key     monomer1   BC key      monomer2        monomer3
```

Two things happen at once:

1. **The old site dies.** Its protospacer becomes `ATGGTGAGCACGNNGGATGA` — the
   PAM-proximal **seed** is destroyed. It can never be re-edited. ⇒ **irreversibility**,
   and it is *structural*, not probabilistic.
2. **A new site is born one monomer downstream.** The inserted `GGA` now sits
   immediately 5′ of monomer2, regenerating `GGA + monomer2 + m3[1-3]` with monomer3's
   `TGG` as PAM — **byte-for-byte the identical 20-bp protospacer.** The same pegRNA
   spacer works again, at position 2.

The insertion carries its own carriage return. The "type guide" (write-head) advances
exactly one monomer per write; the array fills strictly 5′→3′.

#### Key vs. PAM — do not conflate

|  | **PAM** (`TGG`) | **Key** (`GGA`) |
|---|---|---|
| What it is | motif read **directly by the Cas9 protein** (PAM-interacting domain) | 3 bases of the **protospacer** — positions 1–3, the PAM-distal end |
| Position | immediately **3′** of the 20-bp protospacer | immediately **5′** of a monomer; 17–20 bp *upstream* of the PAM |
| Recognized by | protein–DNA contact. **Not present in the pegRNA at all** | **Watson–Crick pairing with nt 1–3 of the pegRNA spacer** |
| Constraint | must be `NGG` — hard SpCas9 requirement | arbitrary; free design choice. TAPE-6 uses a **4-bp** key and has a *lower* sequential error rate |
| Belongs to | **the monomer** (positions 4–6 of `TGATGGTGAGCACG`) | not part of any monomer — it is the separate, transferable token |

The design arithmetic is the whole trick:

$$\underbrace{3}_{\text{key}} \;+\; \underbrace{14}_{\text{monomer}_n} \;+\; \underbrace{3}_{\text{monomer}_{n+1}[1..3]} \;=\; 20\ \text{bp protospacer},\qquad \text{PAM}=\text{monomer}_{n+1}[4..6]$$

A monomer is 14 bp — **exactly 3 bp short** of completing a protospacer for the next
monomer's PAM. The key supplies those 3 bp. Anywhere without a key presents `ACG`
(monomer positions 12–14) in place of `GGA`: a full 3/3 mismatch at protospacer
positions 1–3.

And **the key is written by the edit** — every insert is `[symbol][GGA]`, landing at the
monomer$_n$/monomer$_{n+1}$ junction, so the new `GGA` ends up immediately 5′ of
monomer$_{n+1}$. The key is a token handed one position down the array per write. Same
pegRNA spacer works at every position.

*Why the 3-bp discrimination works so well:* prime editing has a **higher fidelity bar
than nuclease Cas9**, because a productive edit needs not only R-loop formation and
nicking but also stable capture of the nicked 3′ end by the PBS and successful RT
extension. Nuclease Cas9 often tolerates PAM-distal mismatches; prime editing much less
so. The architecture leans on this.

### 0.2b Construct architecture — what is integrated where

**This differs by experiment and it matters.** The three components are *not* three
separate integrations.

| Experiment | TAPE | pegRNA | Prime editor |
|---|---|---|---|
| Typewriter Figs. 1–3 (HEK293T) | piggyBac, stable | **transient plasmid**, re-transfected per "epoch" | **transient plasmid** (PE2, later PEmax + hMLH1dn) |
| Typewriter Figs. 4–5 (lineage tracing) | **lentiviral, co-integrated with its pegRNA**, MOI ≈ 19 | **same cassette as the TAPE** (U6-pegRNA-InsertBC) | Dox-inducible PE2, separately integrated (monoclonal iPE2(+) line) |
| Mouse embryo | **piggyBac, co-integrated with its epegRNA** (`epegRNA::TAPE-BC::circTAPE`) | **same cassette as the TAPE** | PEmax-P2A-dTomato, separate piggyBac construct |

⇒ In both lineage-tracing configurations: **PE separate; pegRNA + TAPE together.**
Choi et al. state the consequence directly: *"pegRNAs expressed by TargetBC-defined
integrants **compete** to mediate insertions at the TAPE-1 arrays **within the same
cell**."*

**Copy numbers, measured:**

- *Prime editor.* Mouse ddPCR: embryo #3 ≈ **7 copies** → extensive editing; embryo #2 ≈
  **2 copies** → modest; embryo #6 → **no detectable PEmax** → editing low-to-absent
  *despite carrying 5 TAPE barcodes*. A clean dose–response, n = 3.
- *TAPE/pegRNA integrations.* Embryo #3: 11, all single-copy. Embryo #2: 17 distinct
  barcodes but **~35 actual integrations** — some barcodes multi-copy from piggyBac
  excision/re-integration. Also seen in high-MOI culture; handled by
  `tape/copynumber.py`.

#### *cis* vs *trans*: which perturbations hit which parameter

**A pegRNA is diffusible and acts in *trans*.** Its copy number and chromatin context
affect the cell-wide pegRNA pool — every tape in that cell equally. This is the key
identifiability boundary.

| Perturbation | Acts in | Data signature | Where it lives in SciPhy |
|---|---|---|---|
| pegRNA cassette copy number | **trans** | changes **which symbol** is written → shifts $\xi_i$, **not** $\lambda$ | insert probabilities $\xi$ ✓ |
| pegRNA cassette chromatin / silencing | **trans** | same — shifts $\xi_i$ for that symbol | $f$ ✓ |
| PE copy number / PE silencing | **trans** | scales **all** tapes in that cell by a common factor | ✗ **not modeled** — would need a per-lineage rate multiplier |
| TAPE locus chromatin / accessibility | ***cis*** | changes that **one** tape's rate | per-tape clock rate $r_z$ ✓ |
| TAPE array length / contraction | ***cis*** | changes that tape's saturation point | $N$ per tape ✓ |

⇒ Differing pegRNA copy number is a **composition** effect, not a rate effect. And the
mouse hematopoietic slowdown is a ***trans*** effect (PEmax silencing / SAMHD1-limited
dNTP supply hitting all 11 tapes in those cells simultaneously) that per-tape $r_z$
structurally **cannot** express. See row A2/A8.

**Nasty coupling:** because pegRNA and TAPE are on the *same* construct, silencing an
integration removes both its symbol from the pool **and** that tape from the readout.
The mouse paper attributes the 20% and 3% recovery of 2 of 11 TAPEs to "epigenetic
silencing of a subset of circTAPE-encoding integrations." **Dropout and rate are
correlated through a shared latent cause** — easy to simulate, awful to write a
likelihood for. (Chromatin dependence of prime editing is directly measured, not
hypothesized: Li et al. 2024, *Cell* 187:2411, "Chromatin context-dependent regulation
and epigenetic manipulation of prime editing" — mouse ref 45.)

#### RNP assembly, mechanistically

PE is translated in the cytoplasm and imported via NLS; pegRNAs are Pol III (U6)
transcripts made in the nucleus and never exported. Loading is competitive and
essentially stoichiometric — any pegRNA can load any PE — so

$$[\text{writer for symbol } i] \;\approx\; [\text{PE}] \times \frac{\text{pegRNA}_i}{\sum_j \text{pegRNA}_j}$$

which is exactly why insertion frequency tracks *relative* pegRNA abundance at
$r = 0.87$–$0.94$. Cas9 is generally limiting relative to guides; apo-Cas9 is unstable
while guide-bound Cas9 is stable. **Two cleanly separated knobs: PE dosage sets the
global rate, guide composition sets the symbol distribution.** Convenient for modeling —
and it is precisely the second knob that ENGRAM hijacks.

### 0.3 Consequences that matter for the phylogenetics

| Property | Why it holds | Why it matters downstream |
|---|---|---|
| **Order is *read*, not inferred** | position in the array = rank order of writing | SciPhy's state space is *ordered tuples*; a huge information gain over unordered scar sets |
| **Irreversible** | edit destroys its own seed region | makes the CTMC upper-triangular; enables the ancestral-state-**intersection** trick in the likelihood |
| **One active site at a time** | only the 5′-most monomer has a key | edits are *serialized* — the process is a walk along a line, not N independent processes |
| **Saturation** | array has finite monomers (5, 6, 12, 20…) | absorbing boundary; the Poisson step in SciPhy needs the "≥ remaining sites" tail mass piled on the last state |
| **Thousands of symbols, one target** | symbol in pegRNA RTT, all spacers identical | rich per-edit information; also the hook for ENGRAM (symbol ↔ signal) |
| **No DSB** | Cas9 H840A nickase | no inter-site deletion within a tape; dropout is a *readout* problem, not a *recording* problem |

### 0.4 Rates — three distinct levels, do not conflate

#### Level 1: symbol-to-symbol (which barcode is inserted)

Two separable components. The paper's `edit score` = log₂(insertion frequency /
pegRNA abundance) is designed exactly to separate them.

- **pegRNA abundance.** Site-1 barcode frequency vs. pegRNA abundance in the plasmid
  pool: **Pearson r = 0.87** (ED Fig. 1d); improves to **r = 0.94** when corrected with
  post-cloning pool abundances. So yes — differential pegRNA supply is a first-order
  driver of *which symbol* gets written. Deliberately exploited: programmes 4 and 5
  encode "signal strength" by mixing pegRNA pairs at 1:3 vs. 1:1 / 1:2 / 1:4 / 1:8 and
  the ratios are recovered from the data (Fig. 2h, ED Fig. 4e,f). **This is the
  conceptual seed of ENGRAM.**
- **Intrinsic sequence bias**, *after* normalizing for abundance. Best `CCGGA` = +0.98,
  worst `TGGGA` = −2.38 → spread of 3.36 log₂ ≈ **10-fold**. Highly reproducible across
  replicates (r = 0.97–0.99) and across positions in the array (r = 0.99), so it is a
  property of the *symbol*, not the site. Mechanism: RTT sequence affects RT
  processivity and 3′-flap stability/annealing; flap sequence affects 5′-flap excision
  and MMR resolution; 3′ extension affects pegRNA folding/stability. epegRNAs
  (structured 3′ motif) narrow the spread: 14/16 NNGGA and 59/64 NNNGGA barcodes within
  a 2-fold range; 1,509/1,908 6N+GGA barcodes within a 4-fold range.

> **Note for SciPhy:** these map onto SciPhy's insert-probability vector
> $f = (f_1, \dots, f_n)$, which the model treats as **shared across all tapes and all
> sites, and constant in time**. Both assumptions are testable and both are things
> ENGRAM would break on purpose.

#### Level 2: tape-to-tape (this is SciPhy's per-tape clock rate)

**Differential pegRNA expression cannot explain this level.** pegRNAs are diffusible and
act in *trans*, so every tape in a given cell sees the same pegRNA pool. Tape-to-tape
rate variation must be *cis*:

- **Chromatin context / integration site.** piggyBac integrates semi-randomly; the
  Typewriter paper attributes uneven TAPE recovery to "site-of-integration effects on
  expression," citing Schep et al. on chromatin context shaping Cas9-induced repair
  outcomes.
- **Copy number.** piggyBac integrations are sometimes duplicated. The mouse paper has a
  dedicated module (`tape/copynumber.py`) inferring integer copy number from read
  abundance.
- **Tape sequence**, if orthogonal tapes are used. In the 48-TAPE screen, efficiency and
  sequential error rate varied substantially by basal spacer / key / monomer. TAPE-27
  had >50% greater efficiency than TAPE-1 but 2× the sequential error rate; TAPE-6
  (same basal spacer as TAPE-1, 4-bp key instead of 3-bp) had a *lower* error rate — a
  longer key buys specificity. *Not a factor in the mouse embryo: all 11 integrations
  are TAPE-1, differing only in integration barcode and genomic locus.*
- **Array length / saturation.** Arrays contract in vivo; observed mouse TargetBCs had
  3, 4, or 5 monomers, and long synthetic arrays (12×, 20×) recovered only 8.4 ± 3.3 and
  12.5 ± 4.3 repeats respectively. Shorter tapes reach the absorbing state sooner.

#### Level 3: cell-to-cell — the messy one

**"Pseudo-processivity."** Unconditional site-1 editing ≈ 6%, but site-2 and site-3
editing *conditional on the preceding site already being edited* ≈ 20% — a ~14% excess of
elongation over initiation, reproduced with episomal tapes and across serial
transfections. Candidate explanations offered: heterogeneous transfection
susceptibility, chromatin context, cell cycle phase. The authors state "the primary
reason remains unclear."

Whatever the cause, the statistical shape is a **latent per-cell (or per-cell × per-tape)
frailty that induces positive correlation between successive edits.** Separately, a
modest 1–10% *decrease* after site 2, possibly because each successive site is "active"
for less time than its predecessor.

**Also observed in vivo (mouse paper), and worse:**
- Editing is **not clock-like in time**: a burst at the 2-cell stage coincident with
  major ZGA (marking the first cleavage unambiguously — a genuine gift for rooting),
  then **quiescence for several divisions**, then resumption around E6 and relative
  stability thereafter.
- Editing rate depends on **cell type**: hematopoietic lineages progressively slow as
  they differentiate. Confounded between PEmax silencing and known dNTP-supply limits on
  prime editing (SAMHD1 / dNTP bottleneck literature, mouse refs 58–60).
- **Missingness is severe**: only ~50% of the full 1.75M cell × 11 TAPE matrix was
  confidently genotyped.

### 0.5 Empirical rate/saturation numbers worth remembering

From the 25-day monoclonal HEK293T expansion (Typewriter Fig. 4):

- 3,257 cells, 13 TargetBC-associated TAPE arrays, **59 editable sites** per cell.
- **39.4 ± 6.3 edits per cell** accumulated over 25 days / ~20 doublings.
- Distribution is **Poisson-like: mean 39.4, variance 40.0.** ⇒ ~2 edits per cell
  division. This is the direct empirical warrant for SciPhy's Poisson step.
- Active type guides declined **13 → 8.6 per cell** over 25 days as tapes saturated ⇒
  the absorbing boundary is *not* a rare corner case; it bites within the experiment.
- Site-specific saturation is strongly graded: site 1 mean 96.8% edited, site 5 mean
  19.7% (Fig. 4c). Note this is a *consequence* of serialization, not of site-specific
  rates.
- Mean pairwise difference between cells: 41.9 ± 5.3 edits; to nearest neighbour
  22.2 ± 6.7.

---

## Running list: DNA Typewriter reality vs. SciPhy model assumptions

Populate as we work through the methods. Each row is a candidate extension point.

| # | SciPhy assumption | Reality (source) | Bearing on ENGRAM / SBI |
|---|---|---|---|
| A1 | Constant editing rate $r_z$ per tape, constant **in time** (strict clock) | ZGA burst → quiescence → resumption (mouse paper, Discussion) | ENGRAM *requires* time-varying rates; this assumption must go |
| A2 | Rate does not depend on cell state | Hematopoietic slowdown (mouse paper) | Rate becomes a function of lineage/state — same machinery as signal-dependence |
| A3 | Insert probabilities $\xi_i$ shared across tapes, sites, and **time** | Symbol frequency ∝ pegRNA abundance × intrinsic bias; abundance is manipulable and, in ENGRAM, signal-driven | **This is the single most important hook.** ENGRAM ⇒ $f_i = \xi_i(t, \text{signal})$ |
| A4 | Equal editing rate across sites within a tape | ~14% elongation-over-initiation excess; 1–10% decline after site 2 | Site-indexed rates, or an explicit latent frailty |
| A5 | Independence across tapes given the tree | Shared per-cell PE/pegRNA supply, shared cell-cycle phase ⇒ frailty couples tapes within a cell | Breaks likelihood factorization over tapes → strong argument for SBI |
| A6 | No dropout / missing-data model (stated limitation) | ~50% of the mouse cell × TAPE matrix missing; circTAPE recovery from nuclei poor | Needs an observation model; easy to *simulate*, hard to write down |
| A7 | No sequential error (edits only ever at the type guide) | 0.6% out-of-order editing in TAPE-1; higher for other TAPEs (TAPE-27 2× worse) | Small, but it is the term that makes exact likelihoods ugly |
| A8 | Rate variation is **per-tape** (*cis*) only | PE copy number / silencing is a ***trans*** effect scaling **all** tapes in a cell at once (embryo #6: no PEmax ⇒ no editing despite 5 TAPEs) | Needs a **per-lineage** rate multiplier orthogonal to $r_z$; also the natural place to put ENGRAM's time-varying rate |
| A9 | Dropout (if added) independent of rate | pegRNA and TAPE are **co-integrated** ⇒ silencing removes the symbol from the pool *and* the tape from the readout — shared latent cause | Correlated missingness; a strong SBI argument since the joint is trivial to simulate, intractable to write down |

---

# Session 1a — The editing model (pp. 10, Definitions 1–3, Eqs. 1–2)

## 1a.1 Notation inventory

| Symbol | Meaning | HEK293T | Mouse |
|---|---|---|---|
| $k$ | **tapes** per cell | 13 | 11 |
| $N$ | **sites** per tape | 5 (three tapes have 4) | 6 |
| $M$ | **alphabet** size | 19 | 8–13 |
| $\Gamma=\{\gamma_1,\dots,\gamma_M\}$ | the alphabet | 19 trinucleotides | 3-nt symbols (+1 4-nt) |
| $\lambda$ | **editing / clock rate**, per tape | inferred | inferred |
| $\xi=(\xi_1,\dots,\xi_M)$ | insert probabilities, $\sum \xi_i = 1$ | inferred (Dirichlet(1.5)) | inferred |

## 1a.2 Definition 1 — sequential susceptibility

> *"Initially, only a tape's first site is susceptible to editing. Subsequent sites are
> edited sequentially: editing at site $i$ is only possible if site $i-1$ is already
> edited."*

Session 0's molecular story as an axiom: only monomer 1 has a `GGA` key 5′ of it; the
edit installs a new key licensing monomer 2 and nothing else.

**Consequence bigger than it looks:** at any instant a tape has **exactly one** editable
site. A tape is not $N$ parallel processes — it is a **single pointer walking a line**.
This is why the state is an ordered *sequence* rather than an $N$-vector, and why the
per-tape rate does not decay as sites fill.

## 1a.3 CTMC primer

### The object

A CTMC is (i) a state space $\Omega$, (ii) a **generator** $Q$ with off-diagonal
$q_{ab}\ge 0$ = rate of $a\to b$ and $q_{aa}=-\sum_{b\ne a}q_{ab}$, (iii) a start state.

**Dynamics:** in state $a$, wait $\text{Exp}(\lambda_a)$ with $\lambda_a=-q_{aa}$, then
jump to $b$ w.p. $q_{ab}/\lambda_a$. Repeat.

Two separable pieces — a **clock** (when) and a **jump chain** (where). SciPhy's entire
factorization follows from this split.

### Why exponential holding times are *forced*

Markov property ⇒ given the present, the future is independent of the past. In
continuous time this **requires** memorylessness,
$P(T>s+t\mid T>s)=P(T>t)$, and the exponential is the unique continuous distribution
with that property. With (say) Gamma holding times, elapsed waiting time would predict
residual waiting time, the state alone would no longer summarise history, and the chain
would not be Markov on $\Omega$.

*Physically reasonable here?* Largely — a well-mixed steady-state RNP pool doing 3D
search with small constant success probability per unit time. No aging, no maturation,
no ratchet; each instant is a fresh independent trial.

*Known failures:* pseudo-processivity (site-2 rate ≈20% conditional on site 1 edited vs.
≈6% unconditional) and the mouse's three-phase burst/quiescence/resumption. Rows **A1**,
**A4**.

### The generator for *this* chain (not written down in the paper)

Write $a\oplus\gamma_i$ for "$a$ with $\gamma_i$ appended". Then

$$q_{a,\,a\oplus\gamma_i} = \lambda \xi_i \quad (|a|<N), \qquad q_{ab}=0 \ \text{otherwise}$$

Total exit rate from any non-absorbing state:
$\lambda_a=\sum_i \lambda \xi_i = \lambda\sum_i \xi_i = \lambda$.
Jump chain picks $\gamma_i$ w.p. $\lambda\xi_i/\lambda=\xi_i$ — a Categorical($\xi$) draw **independent
of the clock**.

> **The whole editing model is: an $\text{Exp}(r)$ clock plus an independent
> Categorical($\xi$) symbol draw, repeated until the tape is full.** The two-step
> factorization of the transition probability is then inevitable, not clever.

Toy chain, $N=2$, $M=2$:

```
                           ┌── rf₁ ──→ (γ₁,γ₁)  ▣
            ┌── rf₁ ──→ (γ₁)
            │              └── rf₂ ──→ (γ₁,γ₂)  ▣
   ∅ ───────┤
            │              ┌── rf₁ ──→ (γ₂,γ₁)  ▣
            └── rf₂ ──→ (γ₂)
                           └── rf₂ ──→ (γ₂,γ₂)  ▣

   exit rate = rf₁ + rf₂ = r   from ∅, (γ₁), (γ₂)
   exit rate = 0               from any ▣  (absorbing)
```

The transition graph is an **$M$-ary tree of depth $N$** rooted at $\varnothing$: one
parent (own prefix), $M$ children, no cycles, no back-edges.

### Why continuous time, not discrete

- **Branch lengths are real numbers** — need $P(a\to b\mid t)$ for arbitrary $t>0$.
- **Editing isn't synchronised to division** — Typewriter notes prime editing continues
  in non-mitotic cells; accumulation is "primarily a function of time."
- **Composability** — $P(t)=e^{Qt}$ gives Chapman–Kolmogorov $P(s+t)=P(s)P(t)$, which is
  what lets you multiply along a root-to-tip path. This is *the* reason phylogenetics
  uses CTMCs.

### How this chain differs from a standard substitution model

SciPhy is deliberately *shaped* like one (so `SciPhySubstitutionModel` can slot in where
Jukes–Cantor goes), but violates the usual assumptions:

| Standard (JC69, GTR…) | SciPhy editing model |
|---|---|
| $|\Omega|=4$ | $|\Omega|\sim10^5$–$10^6$ |
| **Reversible**: $\pi_a q_{ab}=\pi_b q_{ba}$ | **Not reversible**: $q_{ba}=0$ whenever $q_{ab}>0$ |
| Ergodic, stationary $\pi$ | **Absorbing**; all mass ends at length-$N$ states |
| Equilibrium frequencies meaningful | Meaningless — only equilibrium is saturation |
| Root placement irrelevant (pulley principle) | **Root placement matters** ⇒ Eq. 10 pins origin at $\varnothing$ |
| Homoplasy via back-substitution | No back-mutation; only convergent forward edits |

## 1a.4 Definition 2 and the "bounded Poisson process"

$$X(t) = \min\big(N(t),\,N\big), \qquad N(t)\sim\text{Poisson}(\lambda t)$$

**View A — truncated count.** Ordinary Poisson body, piled-up tail:

$$P(X(t)=j)=e^{-\lambda t}\frac{(\lambda t)^j}{j!}\ (j<N),\qquad P(X(t)=N)=1-\sum_{j=0}^{N-1}e^{-\lambda t}\frac{(\lambda t)^j}{j!}$$

All mass that "would have been" at $N,N{+}1,\dots$ collapses onto $N$. → Eq. 5.

**View B — absorbed CTMC.** Pure-birth on $\{0,\dots,N\}$, rate $\lambda$ for $j<N$, rate 0 at
$N$. Same object; A is how to compute, B is what it is.

### Why the count is *exactly* Poisson

**Every non-absorbing state has the identical exit rate $\lambda$.** That is what makes jump
times a genuine Poisson process. A rate depending on tape fullness would give a general
pure-birth process and none of the clean formulas.

What guarantees constant exit rate? **The sequential architecture** — exactly one active
site at all times.

*Contrast, non-sequential recorder:* $N$ independent sites at per-site rate $\rho$ ⇒
total rate $\propto$ (unedited sites remaining), which **declines**. Measured in mouse
Fig. 3D: idealized non-sequential $1\times66$ (red) falls **~6-fold**, 8.8 → 1.5
edits/day; sequential (blue) holds flat at ~4.1/day.

> **"Bounded Poisson" = constant rate $\lambda$ until a hard wall, then zero.**
> Constant-until-saturation is the design goal of DNA Typewriter and the reason the maths
> is this clean.

### Why per-cell totals need not look Poisson

1. **Truncation breaks additivity.** Independent Poissons sum to Poisson; independent
   *truncated* Poissons do not.
2. **Tapes are heterogeneous** — different $r_z$, different $N_z$ (three HEK293T tapes
   have $N=4$).
3. **Cells share ancestry** — a sample of cells is not i.i.d.; sister cells share nearly
   all edit history. Empirical histogram can be multimodal though each marginal is not.

Empirically it often does look Poisson anyway: Typewriter Fig. 4d mean 39.4 / var 40.0
over 3,257 cells; mouse Fig. 2C unimodal at 56/66. **Relevant if you want a Poisson-based
summary statistic for SBI.**

## 1a.5 Why $k \ne M$

|  | HEK293T | Mouse |
|---|---|---|
| $k$ (tapes) | 13 | 11 |
| $M$ (symbols) | **19** | **13** at ≥0.5% (8 at sites 3–6) |

The premise to drop: co-integration of pegRNA and TAPE does **not** pair tapes with
symbols, because **pegRNAs act in *trans***. The alphabet at any tape is the union of all
pegRNAs in the cell. Co-integration is a delivery convenience, not a functional coupling.

- **HEK293T ($M>k$):** MOI ≈19 → ~19 integrations, ~19 InsertBCs; SciPhy retained only
  the **13 most frequent** TargetBCs. Tapes filtered, symbols not.
- **Mouse ($M>k$), lovely reason:** complex library injected, only 11 integrated, yet 12
  symbols at site-1, 9 at site-2, 8 at sites 3–6. *"the additional symbols at sites 1-2
  were likely mediated by epegRNAs transiently expressed shortly after PNI, and are thus
  restricted to the 5′-most sites."* ⇒ **$M$ was larger early and shrank** — a
  time-varying alphabet the model assumes away.

**Different structural roles:**
- $k$ = **replicate recording channels**. Likelihood *factorizes* over $k$ ("embarrassingly
  parallel").
- $M$ = **alphabet size**. Enters only *inside* a single transition probability, via $\xi$.

Capacity: $M^N$ per tape, $M^{kN}$ per cell. Mouse $8^6\approx2.6\times10^5$ per tape,
$8^{66}\approx10^{59}$ per cell — the number the mouse paper quotes.

> **For ENGRAM the decoupling is by design:** $M$ becomes the *signal* axis
> (signal-responsive promoters driving $M$ reporter pegRNAs), $k$ stays the *replication*
> axis (passive substrates). $M$ then has nothing to do with $k$.

## 1a.6 Definition 3 — irreversibility

Molecular basis: the 5-bp insertion lands 3 bp from the PAM and destroys the
PAM-proximal seed. **Structural**, not a low-probability event.

Formally: $q_{ab}=0$ unless $a\sqsubseteq b$ ($a$ a **prefix** of $b$). $Q$ is "upper
triangular" w.r.t. the prefix partial order.

1. **Absorbing state.** $|a|=N \Rightarrow \lambda_a=0$. The "wall".
2. **No back-mutation.** Convergence still possible (two lineages independently writing
   $\gamma_3$ at site 1); reversion is not. Homoplasy is one-directional.
3. **The intersection trick.** Children with ancestor sets $A_k,A_j$ ⇒ parent lies in
   $A_k\cap A_j$, because an edit in one child but not the other *cannot* have been in
   the parent. **Eq. 9 — the load-bearing beam of the likelihood.** → Session 1c.

## 1a.7 The state space (Eq. 2)

$$\Omega=\Big\{(s_i)_{i\le n} : s_i\in\Gamma,\ 0\le n\le N\Big\}$$

$\varnothing$ = unedited start for all $k$ tapes; every length-$N$ tuple absorbing.

**Ordered sequences, not sets or multisets.** $(\gamma_3,\gamma_7)\ne(\gamma_7,\gamma_3)$.
That inequality *is* the information gain over GESTALT-style unordered recorders.

$$|\Omega|=\sum_{n=0}^{N}M^n=\frac{M^{N+1}-1}{M-1}$$

| Dataset | $M$ | $N$ | $\vert\Omega\vert$ |
|---|---|---|---|
| HEK293T | 19 | 5 | **2,613,660** |
| Mouse-like | 8 | 6 | **299,593** |

A naive $e^{Qt}$ on 2.6M states, per branch, per tape, inside MCMC — hopeless.
Tractability comes entirely from (a) prefix structure (never build $Q$) and (b)
intersection pruning (never enumerate $\Omega$).

## 1a.8 The $a$, $b$, $c$ notation

$Y_t=a$, $Y_{t+\Delta t}=b$, $c=b\setminus a$.

**Not ordinary set difference.** Irreversibility forces $a\sqsubseteq b$, so $c$ is the
**suffix** of $b$ beyond position $|a|$. E.g. $a=(\gamma_2)$, $b=(\gamma_2,\gamma_1,\gamma_3)$
⇒ $c=(\gamma_1,\gamma_3)$, $|c|=2$ (Fig. 6b).

Why bother: $(a,b)$ carries the same information as $(a,c)$, and $c$ splits cleanly into
"how many" ($|c|$, the clock) and "which ones in what order" ($c$, the jump chain) —
**independent by construction of the generator**:

$$P_{a,b}(\Delta t;\xi,\lambda)=\underbrace{P(|c|\mid\Delta t,\lambda)}_{\text{bounded Poisson}}\times\underbrace{P(c\mid |c|,\xi)}_{\text{product of }\xi_i}$$

→ Session 1b.

## 1a.9 ⚑ The ENGRAM hook is inside Definition 2b

> *"We additionally assume that insertions operate independently — an insertion at a
> previous site does not influence the insertion probability at the currently active
> site."*

In generator terms: **$\xi_i$ does not depend on $a$.** Symbol draws are i.i.d.
Categorical($\xi$) at every step.

**This is exactly the assumption ENGRAM must break, and the way it breaks *is* the
signal.** If $\xi_i(t)$ varies in time, symbols at successive sites become *marginally*
dependent, because consecutive sites are written at nearby times and $f$ varies smoothly.
**The dependence structure between adjacent positions is the temporal signal trajectory.**

DNA Typewriter already demonstrates this in reverse: **bigram frequency matrices**
(Fig. 2c–g) are strongly non-uniform precisely because the pegRNA pool was changed
between epochs, and off-diagonal asymmetry recovers the true epoch ordering (ED Fig. 4c;
robust down to 2,500 reads). Proof-of-concept that temporal $f$ variation is real and
decodable. ENGRAM swaps the pipette for a biological signal.

### ⚠ Do not assume intractability prematurely

Conditional on $n$ events in $[0,T]$, Poisson event times are distributed as order
statistics of $n$ i.i.d. Uniform$[0,T]$. So

$$P(\text{symbols}\mid n,T)=\int_{0<t_1<\cdots<t_n<T}\prod_{j=1}^{n}f_{s_j}(t_j)\,\frac{n!}{T^n}\,dt_1\cdots dt_n$$

— an $n$-dimensional ordered-simplex integral with $n\le N=6$. **May well admit a
recursion or quadrature.** The likelihood-based route may survive the ENGRAM extension
further than expected. Hold SBI in reserve for what genuinely resists it: correlated
dropout (**A9**), per-lineage rate multipliers (**A8**), cell-state-dependent rates
(**A2**).

---

# Session 1b — Transition probabilities (Eqs. 3–7)

## 1b.1 Step 1 — how many

$$P(|c| = j;\, r, \Delta t) = \begin{cases} \dfrac{(r\Delta t)^j}{j!}e^{-r\Delta t} & |b| < N \quad \text{(Eq. 3)}\\[2ex] 1 - \displaystyle\sum_{m=0}^{j-1}\dfrac{(r\Delta t)^m}{m!}e^{-r\Delta t} & |b| = N \quad \text{(Eq. 4)}\end{cases}$$

Eq. 4 is just $P(N(\Delta t)\ge j)$ — the upper tail piled onto the last state.

**It is a proper pmf** over $j\in\{0,\dots,N-|a|\}$, despite looking like two glued-together
formulas: body terms sum to $\sum_{m=0}^{N-|a|-1}\text{Pois}(m)$, saturating term is
$1-\sum_{m=0}^{N-|a|-1}\text{Pois}(m)$, total exactly 1.

## 1b.2 Step 2 — which ones

$$P(c;\,|c|=j,\xi)=\prod_{m=|a|+1}^{|a|+j}\prod_{l=1}^{M} f_l^{\,\mathbf{1}(\gamma_l,\,s_m)}\quad\text{(Eq. 6)}$$

The inner product over $l$ is **cosmetic** — the one-hot indicator kills all but one factor,
so this is just $\prod_m f_{s_m}$. Written this way for implementation and because
$\partial/\partial f_l$ is then clean for MCMC. **Definition 2b (independence) is what lets
it factorize.**

## 1b.3 Step 3 — the product

$$P_{a,b}(\Delta t;\xi,\lambda)=P(|c|;\Delta t,\lambda)\times P(c;|c|,\xi)\quad\text{(Eq. 7)}$$

and implicitly $P_{a,b}=0$ whenever $a\not\sqsubseteq b$.

**Two observations:**
1. $P_{a,b}$ depends on $a,b$ **only through $|a|$, $|b|$, and the suffix $c$** — never on
   *which* symbols are already in $a$. The generator's homogeneity, surfacing.
2. This factorization *is* the clock/jump-chain split from §1a.3. Nothing new assumed.

---

# Session 1c — The likelihood (Eqs. 8–12)

## 1c.1 The tree

### Why exactly $n-1$ internal nodes

**Edge counting.** Every node but the root has one parent edge ⇒ $n+I-1$ edges. If every
internal node has exactly 2 children, counting the same edges from above gives $2I$. So
$2I = n+I-1 \Rightarrow I = n-1$.

**Coalescence counting (better intuition).** Backward in time you start with $n$ lineages and
end with 1. A **binary** merge reduces the count by exactly **1**. Going $n \to 1$ therefore
needs exactly $n-1$ merges.

⚠ **This depends entirely on the tree being rooted and fully bifurcating.** The paper states
it as if definitional; it is a modeling assumption. With polytomies, a node with $d$ children
reduces the count by $d-1$, so

$$\sum_{i\,\in\,\text{internal}} (d_i - 1) = n-1$$

Any $d_i>2$ forces $I<n-1$ (star tree: $I=1$). Every count in the complexity analysis becomes
an **upper bound** rather than an equality.

### Will the tree be binary? Four different answers

1. **True cell lineage — yes, exactly.** Mitosis gives exactly two daughters. Binary by
   biology, not assumption. Nicer than macroevolution/epidemiology where bifurcation is a
   convenient fiction.
2. **Sampled tree — yes, almost surely.** Internal nodes are divisions where two *sampled*
   sublineages diverge. One division can't yield three daughters, and two divisions at the
   *identical* instant has probability 0 in continuous time. Binarity survives subsampling
   w.p. 1.
3. **Inside BEAST — yes, by construction.** `Tree` is strictly bifurcating; operators
   (subtree slide, Wilson–Balding, narrow/wide exchange) preserve it; the birth–death prior
   is a density on binary ranked trees.
4. **In practice for data like this — NO.** Binary is identifiable *in principle* but not
   *resolvable from the data*. Mouse: query cells share a mean of **0.6 (A) / 1.2 (B)**
   private edits with their anchor; most are **genotype-identical at every recovered site**.
   Hence: *"Query cells with identical genotypes are represented as a single polytomy rather
   than being arranged into an arbitrary ladder."* 65,879 terminal polytomies, **max fan-out
   168**; LSD2 run with *"near-zero input branches collapsed to polytomies at `-l 0.01`."*

> **SciPhy assumes binary; the mouse tree explicitly is not.** Real friction point for
> applying SciPhy machinery to that dataset or benchmarking against it.

### Aside: what a polytomy is, and what "collapsing" means

**Polytomy** (= multifurcation): an internal node with **≥3 children**.

- **Hard polytomy** — a genuine simultaneous multi-way split.
- **Soft polytomy** — a binary split really happened; the data can't order it. An admission
  of ignorance, not a claim about biology.

> **In a cell lineage tree every polytomy is soft by construction**, since mitosis yields
> exactly two daughters. There is always a true binary answer; the polytomy declines to guess.

**Mechanically, collapsing = edge contraction.** Four cells with identical genotypes; a binary
tree is forced to commit:

```
        ┌── A                              ┌── A
     ┌──┤                                  ├── B
  ┌──┤  └── B      ──contract──▶           ├── C
  │  └───── C      zero-length             └── D
  └──────── D      internal edges
                                        degree-4 polytomy
  arbitrary ladder:
  all 15 rooted binary topologies on 4 tips fit equally well;
  internal branches carry ZERO edits
```

Contracting an edge deletes it and merges its endpoints, promoting the child's children to the
parent. Same information, no invented structure.

**Mouse paper does this at two stages:**
1. **Placement** — *"Query cells with identical genotypes are represented as a single polytomy
   rather than being arranged into an arbitrary ladder."*
2. **Dating** — LSD2 with `-l 0.01`: *"near-zero input branches collapsed to polytomies."*

**Scale:** 361,162 cherries vs. **65,879 polytomies (up to 168 tips)**. Since a polytomy of $m$
tips contributes $\binom{m}{2}$ sibling pairs, polytomies generate **598,053** pairs — more
than all cherries combined.

**Why it matters for SciPhy:** BEAST cannot represent a polytomy (strictly binary tree object).
SciPhy encodes the identical ignorance as a **posterior spread over many binary resolutions
with near-zero branch lengths** — more principled, since it quantifies rather than erases the
uncertainty, but exactly what makes MCMC crawl. The 15 equally-good topologies on 4 identical
cells become $10^{2860}$ near-equally-good topologies on 1,000 mostly-identical cells, with no
gradient to follow.

> Same biological fact, two encodings: **collapse it** (cheap, honest, lossy) or **integrate
> over it** (expensive, honest, complete).

### How the two frameworks handle unresolvability

| | Mouse (NJ + LSD2) | SciPhy (Bayesian) |
|---|---|---|
| Unresolved region | **Collapsed** to explicit polytomy | Kept binary; uncertainty in the **posterior** |
| Reported | Polytomous point estimate | Distribution over binary topologies, summarized by CCD0 |
| Honesty | Explicit about the unknown | Also honest, and quantitative |
| Cost | Cheap | **Enormous** — MCMC must explore a plateau of near-equiprobable topologies |

> ⚑ **This is the real scalability bottleneck.** $(2n-3)!!$ rooted binary trees on $n$ labeled
> leaves ≈ $10^{2860}$ at $n=1000$. When many cells are genotype-identical, vast regions have
> *nearly identical* posterior density and the chain has nothing to climb. SciPhy caps out
> ~1,000 cells at $10^9$ iterations. **The bottleneck is topology mixing, not likelihood
> evaluation** (which is already linear in $n$). "Make the likelihood faster" attacks the
> wrong thing.

### Node numbering

Paper says $n-1$ internal nodes but numbers them $n{+}1,\dots,2n$ — $n$ labels. Resolves as:

$$\underbrace{n}_{\text{tips }1..n}+\underbrace{(n-1)}_{\text{internal }n+1..2n-1}+\underbrace{1}_{\text{origin }2n}=2n$$

Consistent with the equations: Eq. 11 sums $a_{n+1},\dots,a_{2n-1}$ (the $n-1$ genuine internal
nodes); node $2n$ is fixed at $\varnothing$ by Eq. 10; the product $\prod_{i=1}^{2n-1}$ runs
over every node that *has* a parent edge (the origin has none).

*Hidden implementation detail:* ordering internal nodes by **increasing age** guarantees every
node's children have smaller indices. So iterating $i=n{+}1,n{+}2,\dots$ **is** a valid
post-order traversal — no recursion or sort needed for the pruning pass.

### Origin vs. root

```
   ORIGIN  (node 2n)     t = 0    experiment starts; state = ∅ by definition
     │                            NOT a branching event — a boundary condition
     │
     │  stem branch (τ_stem)      edits accumulate here
     │                            → shared by EVERY sampled cell
     │
   ROOT  (node 2n−1)     t = t_MRCA    first split among SAMPLED lineages
    ╱ ╲                                IS a branching event (a real cell division)
   ╱   ╲
  ⋯     ⋯
  │ │   │ │
 tips (1..n)            t = T     all sampled simultaneously → ultrametric
```

**Root** = MRCA of the *sampled* cells; an actual cell division; its time is a random variable
inferred jointly with everything else.
**Origin** = start of the *process*; one founding lineage (monoclonal founder / zygote) with an
unedited tape; not a division — it's where the CTMC is initialized.

**Why both are needed:**

1. **Different times.** The founder predates the first divergence *among sampled cells*.
   Sparse sampling, lineage extinction, and the waiting time to the first informative split
   all push $t_\text{MRCA}$ later than $t_\text{origin}$.
2. **The stem carries the shared edits.** The root is **not** unedited; whatever accrued on
   the stem is in every sampled cell. $P_{\varnothing,\,a_\text{root}}(\tau_\text{stem})$ is
   exactly that term. *Mouse illustration:* the ZGA burst wrote ~20 edits per blastomere
   shared by all descendants; they date A/B subtree roots at **E1.5** and *"joined [them] onto
   a shared day-0 zygote root by a 1.5-day stem."* That stem **is** the origin→root branch.
3. **The tree prior conditions on one or the other, and the densities differ.** Origin ⇒
   process began with **one** lineage at $t_\text{or}$; root ⇒ began with **two** at
   $t_\text{MRCA}$. SciPhy conditions on the origin — and then **fixes** it: Tables 3–4 list
   *"Origin or experiment duration: 25 [fixed]"* (HEK293T), *"11 [fixed]"* (gastruloid).
4. **Fixing the origin is what makes $\lambda$ identifiable in absolute units.** Rate and time enter
   the likelihood almost entirely through their product; classical phylogenetics breaks this
   with fossil calibrations or serially sampled tips. Here you break it by *knowing the
   experiment ran 25 days*. That converts the clock rate into **edits per day** and is why the
   prior reads *"between one and ten edits per tape over 25 days."*

> Structural advantage of lineage recording over epidemiology/macroevolution: **the origin is
> normally a diffuse nuisance parameter; in a designed experiment it is a known constant.**

**Ultrametricity.** All cells harvested at one time point (day 25, E13.5) ⇒ all tips
contemporaneous ⇒ **ultrametric tree**. In the birth–death prior this is $\rho$-sampling at the
present, sampling proportion fixed (0.0008 HEK293T, 0.0195 gastruloid). No serial sampling, no
tip dates to infer — another simplification the experimental design buys, inherited by any
ENGRAM extension.

### Factorization over tapes

Tapes independent given the tree (**A5**) ⇒ likelihood = product of $k$ per-tape densities.
Everything below is for one tape.

## 1c.2 Eqs. 8–9 — the state-set reduction

### What problem they solve

Eq. 11 sums over every joint assignment of states to internal nodes: $|\Omega|^{n-1}$ terms
with $|\Omega|=2{,}613{,}660$. Eqs. 8–9 shrink the per-node candidate set to $\le N+1$.

Character of the step, worth noting up front:
- **Purely combinatorial** — no $\lambda$, no $\xi$, no branch lengths. Depends only on tip data
  and topology.
- Therefore a **preprocessing pass**, once per (tape, topology).
- ⇒ **MCMC moves on $\lambda$, $\xi$, branch lengths leave all $A_i$ untouched; only topology moves
  force recomputation.** A meaningful caching boundary.

### Eq. 8 — the leaf sets

$$A_m=\big\{\varnothing,\ (s_1),\ (s_1,s_2),\ \dots,\ (s_1,\dots,s_{|v^m|})\big\}$$

**Necessity (nothing else belongs).** If $a$ with $|a|=p$ was an ancestor of $v^m$: Def. 1
says edits only *append* at the single active site; Def. 3 says they are never removed or
overwritten. So positions $1..p$ of $v^m$ still hold exactly what $a$ held, and $p{+}1..L$
were appended after. Hence $a=(s_1,\dots,s_p)$. Non-prefixes are **impossible**, not unlikely.

**Sufficiency (nothing can be dropped).** For every $p$, $P_{(s_1..s_p),v^m}(\Delta t)>0$ for
any $\Delta t>0$ (positive Poisson term × $\prod f$, all strictly positive). $|A_m|=L+1$.

⚠ **"Power set" is loose wording.** A literal power set of the edit multiset has $2^L$
elements (32 vs. 6 prefixes at $L=5$), and most members — e.g. $\{s_2,s_4\}$ — aren't valid
chain states at all. **The equation is correct and unambiguous; read the equation.**

**Concrete ($N=6$):** $v^m=(\text{AAG},\text{CAC},\text{GCC},\text{CTG})$, sites 5–6 unwritten:

$$A_m=\{\varnothing,\ (\text{AAG}),\ (\text{AAG},\text{CAC}),\ (\text{AAG},\text{CAC},\text{GCC}),\ (\text{AAG},\text{CAC},\text{GCC},\text{CTG})\}$$

The leaf's **own state is included** — its parent may have been identical, probability
$e^{-r\tau}$ of no edits on the branch.

> **⚠ Partial editing is NOT missing data.** 4 of 6 sites written means sites 5–6 are
> *observed to be unwritten* — a positive datum contributing likelihood. Categorically
> different from a tape that failed to sequence. SciPhy handles the first natively, the second
> not at all (**A6**). Easy way to misread the model.

> **The whole ballgame:** $|\Omega|=2{,}613{,}660$; each leaf collapses its ancestors to
> **≤ 6 states**. Five orders of magnitude, bought entirely with Definitions 1 and 3.

### Eq. 9 — the intersection

$$A_i = A_k \cap A_j$$

**Forward.** If $a$ is the state at $i$, the state at child $k$ is some $b\sqsupseteq a$ with
$b\in A_k$. Since $A_k$ is **prefix-closed**, $a\sqsubseteq b\in A_k \Rightarrow a\in A_k$.
Same for $j$. (Prefix-closure does real work — the argument fails for arbitrary candidate sets.)

**Reverse.** If $a\in A_k\cap A_j$, then $a$ is a consistent ancestor of everything below both
children, hence a possible state at $i$. Equality.

**The two verbal clauses:**
- *Pruning direction* — an edit in $k$ but not $j$ cannot have been in the parent;
  irreversibility would have forced $j$ to inherit it.
- *Anti-pruning direction* (**more important**) — shared edits do **not** license concluding
  the parent had them. Convergent independent writes give the same pattern.

> **Convergence is common, not a corner case.** P(two independent lineages write the same
> symbol at a site) $=\sum_i \xi_i^2 \ge 1/M$ (Cauchy–Schwarz, equality iff uniform $\xi$). Mouse
> $M=8$ ⇒ **≥12.5% per site**, and higher since $f$ is measurably non-uniform (10-fold edit-score
> spread). Homoplasy at the first shared site is routine.

**Parsimony vs. likelihood.** The mouse branch-support statistic takes the *maximum* element of
$A_i$: *"If two descendant lineages differ at a site, that site and all downstream sites are
inferred to be unedited in their common ancestor."* SciPhy keeps the whole set and lets the
likelihood weight every member. **Parsimony picks the top of the set; likelihood integrates
over it.** At a 12.5% collision rate this is not cosmetic.

```
child k:  (AAG, CAC, GCC, CTG)   A_k = {∅, (AAG), (AAG CAC), (AAG CAC GCC), (AAG CAC GCC CTG)}
child j:  (AAG, CAC, CCC)        A_j = {∅, (AAG), (AAG CAC), (AAG CAC CCC)}
          ───────────
   lcp =  (AAG, CAC)             A_i = {∅, (AAG), (AAG CAC)}        |A_i| = 3
```

Diverged at site 3 (GCC vs CCC) ⇒ parent certainly had not written site 3. Whether it wrote
0, 1, or 2 sites is genuinely undetermined — AAG and CAC could each be inherited *or* written
twice independently.

### The lcp reformulation (not in the paper)

Every $A$ is $\text{pre}(u)$ for a single string $u$, and

$$\text{pre}(u)\cap\text{pre}(w)=\text{pre}\big(\mathrm{lcp}(u,w)\big)$$

By induction, $u_i=\mathrm{lcp}(u_k,u_j)$. **The entire Eq. 8/9 computation is: compute the
longest common prefix of descendant leaves at every internal node.** One post-order pass.

**Representational payoff:** never materialize a set of tuples. Store **one integer per node**
(the lcp length) plus a pointer to any descendant leaf, since $u_i$ is a prefix of that leaf's
tape. $O(1)$ memory/node, $O(N)$ time/intersection.

**Two properties:**
1. **The restriction is EXACT, not an approximation.** States outside $A_i$ contribute *exactly
   zero* — some $P_{a,b}$ in the product is 0 by irreversibility. Not truncating a sum;
   skipping proven-vanishing terms.
2. **Monotone rootward.** lcp lengths non-increasing. $u_\text{root}$ = lcp of **all** sampled
   cells. *Mouse:* founder signatures of **23 edits (A) / 19 (B)** across 11 tapes ≈ 2 per tape
   ⇒ $|A_\text{root}|\approx 3$ per tape. That shared prefix is what let them partition A from B
   directly instead of trusting a distance method with the first split.

$\varnothing$ is in every prefix set, so Eq. 10's $A_{2n}=\{\varnothing\}$ is always consistent
with $A_\text{root}$. It is a genuine modeling claim, though — **it asserts the founder was
unedited.** True for a zygote; verify before assuming it for a founder isolated mid-experiment.

Cost: $n-1$ intersections × $O(N)$, per tape → $O(knN)$.

### ⚠ Two failure modes

**1. A single site-1 error is catastrophic and propagates to the root.** One miscalled symbol
at position 1 in leaf $x$ ⇒ $\mathrm{lcp}(x,\text{anything})=\varnothing$ ⇒ since lcp lengths
are non-increasing and $\varnothing$ is absorbing, **$u=\varnothing$ at every node on the path
from $x$ to the root.** The whole ancestral-state structure along that path collapses.

The likelihood doesn't go to zero — it goes *confidently wrong*. A cell sharing no edits with
anyone looks early-diverging, gets attached near the root, and drags topology and node dates
with it. **The worst place for an error is site 1, and it produces a spurious deep branch** —
the most damaging error type in a rooted tree.

Hence the mouse pipeline's paranoia: chimera removal, cross-barcode swap filter, ≥2 molecules
per UMI, ≥3 molecules per locus call, ≥0.9 dominance. 0.6% sequential error enters through a
**hard constraint**, not a soft weight.

**2. Missing tapes break the representation, not the mathematics.** Correct handling is
$A_m=\Omega$ with uninformative tip likelihood $L_m(a)=1$ — standard phylogenetic
marginalization; the intersection passes through ($A_k\cap\Omega=A_k$). But $\Omega$ is
prefix-closed *without* being the prefix set of any single string, so the one-integer
representation needs a wildcard extension. Nuisance, not barrier.

The real reason they punt: **dropout is non-random.** *"Ignoring such sources of loss may still
bias tree or population parameter inference by violating a common assumption of birth-death
sampling models regarding uniform sampling of the population of interest, e.g., if DNA Tapes
are simultaneously lost for groups of related cells."* Correlated loss (**A9**, heritable
silencing) breaks the **tree prior**, not the pruning. Naive marginalization would be worse
than dropping cells, because it would look correct.

### ⚑ Implications for the extensions

**Eqs. 8–9 are invariant under the ENGRAM extension.** They invoke only Definitions 1 and 3 and
never touch $\lambda$ or $\xi$. Time-varying $\xi(t)$ leaves the state sets, the lcp recursion, the
$O(1)$ representation, and the $O(knN^2)$ bound entirely alone. **Third independent
confirmation** that the signal-history extension changes exactly one function, $P_{a,b}$.

**But a probabilistic error model directly conflicts with this machinery.** If out-of-order
editing has probability $\varepsilon>0$, non-prefix states acquire positive probability, $A_m$
inflates toward $\Omega$, and the five-orders-of-magnitude reduction evaporates. **Efficiency
here rests on a *hard* combinatorial constraint; softening it destroys it.**

| Extension | Compatible with the Eq. 8/9 reduction? |
|---|---|
| Time-varying $\xi(t)$ (ENGRAM) | **Yes** — untouched |
| Time-varying $r(t)$ (non-clock recording) | **Yes** — untouched |
| Per-lineage rate multipliers (**A8**) | **Yes** — untouched |
| Missing-tape marginalization (**A6**) | Yes, with a wildcard representation |
| Sequential error $\varepsilon>0$ (**A7**) | **No** — reduction collapses |
| Correlated dropout (**A9**) | Reduction survives; the **tree prior** breaks |

⇒ The likelihood route survives everything in the top half. **SBI earns its keep specifically
on the bottom two rows** — where simulation is indifferent to exactly the constraints that make
the exact likelihood tractable.

## 1c.3½ Refresher — likelihood vs. probability

**One function, $p(D\mid\theta)$, read two ways:**

| Fix | Vary | Called | Normalizes? |
|---|---|---|---|
| $\theta$ | $D$ | **probability (density)** | Yes — sums/integrates to 1 over $D$ |
| $D$ | $\theta$ | **likelihood** $L(\theta)$ | **No** — nothing forces $\int L(\theta)d\theta=1$ |

Probability asks *"if $\theta$ were true, how often would I see each dataset?"*; likelihood asks
*"given the one dataset I have, how well does each $\theta$ account for it?"*

1. **Likelihood is not a distribution over $\theta$.** Probability statements about parameters
   require Bayes.
2. **Only ratios matter** — defined up to a multiplicative constant.
3. **Notation is routinely abused**, including here.

### Applied to SciPhy

$D$ = observed tape states at the tips (cell × tape × site matrix).
$\theta$ = $T$ (tree: topology **and** node times), $\lambda$ (per-tape clock rate), $\xi$ (insert
probabilities).

> **The tree is a parameter** — high-dimensional, discrete *and* continuous at once. This is
> why the bottleneck is topology mixing, not likelihood evaluation.

Eq. 10 is itself an instance of the ambiguity:

$$\underbrace{\text{Lik}(T,\xi\mid D)}_{\text{likelihood reading: fn of }T,\xi}=\underbrace{P(D\mid T,\lambda,\xi,A_{2n})}_{\text{probability reading: fn of }D}$$

The `|` on the left is borrowed notation — $T,\xi$ are not random variables being conditioned.
Cleaner: $L(T,\xi\,;D)$.

### Where it sits in the full model

$$p(T,\lambda,\xi,b,\delta\mid D)\ \propto\ \underbrace{P(D\mid T,\lambda,\xi)}_{\textbf{Eqs. 10–12}}\times\underbrace{p(T\mid b,\delta,\rho)}_{\text{tree prior}}\times\underbrace{p(\lambda)p(\xi)p(b)p(\delta)}_{\text{Tables 2–4}}$$

- **Editing model → the likelihood.** Says nothing about how likely a *tree shape* is.
- **Birth–death-sampling → the prior on trees.** Its parameters are biologically meaningful
  (birth = division rate, death = death rate, $\rho$ = sampling proportion), so inferring
  $b,\delta$ **is** the phylodynamic result. In ordinary phylogenetics the tree prior is a
  nuisance; here it is half the point.

**Marginalized vs. inferred:** ancestral tape states $a_{n+1},\dots,a_{2n-1}$ are **summed out**
(latent nuisance variables). Tree, $\lambda$, $\xi$ are **inferred**. This distinction becomes
load-bearing in §1c.5.

## 1c.4 Eqs. 10–12 — Felsenstein pruning

### Eq. 10 — the statement

$$\text{Lik}(T,\xi\mid D)=P(D\mid T,\lambda,\xi,A_{2n}),\qquad A_{2n}=\{\varnothing\}$$

Everything conditional on the origin being unedited. **Imposed, not derived.**

### Eq. 11 — the honest, unusable form

$$\text{Lik}=\sum_{a_{2n-1}\in A_{2n-1}}\!\!\cdots\!\!\sum_{a_{n+1}\in A_{n+1}}\ \prod_{i=1}^{2n-1}P_{a_{\pi_i},\,a_i}(\tau_i,\xi)$$

- **Sums** range over internal nodes $n{+}1..2n{-}1$ ($n-1$ of them). Tips fixed at observed
  states; node $2n$ fixed at $\varnothing$.
- **Product index $i=1..2n-1$** covers every node that *has* a parent = every branch. $2n$ nodes
  ⇒ $2n-1$ edges; only the origin lacks an incoming branch. Branches indexed by their child.
- $\pi_i$ = parent of $i$; $\tau_i$ = branch length above $i$; $P_{a_{\pi_i},a_i}$ = Eq. 7.

**Why a product?** For a fixed complete assignment,

$$P(\text{all states})=\underbrace{P(a_{2n})}_{=1,\ \text{fixed}}\prod_{i=1}^{2n-1}P(a_i\mid a_{\pi_i})$$

by the Markov property on the tree. **Each term of Eq. 11 is the probability of one complete
evolutionary history**; the sum is over all histories consistent with the data.

**Unusable:** $\prod_i|A_i|\le(N+1)^{n-1}$. At $N=5,n=1000$: $6^{999}\approx10^{777}$.

### Eq. 12 — pruning

**Key observation:** in that product, an internal node's state $a_j$ appears in **exactly three
factors** — once as a child ($P_{a_{\pi_j},a_j}$), twice as a parent ($P_{a_j,a_{c_1}}$,
$P_{a_j,a_{c_2}}$). Everything else is constant w.r.t. $a_j$. By the distributive law,
$\sum_{a_j}$ pushes inside to act on those three alone. Bottom-up, each node's sum is a small
local computation on numbers already available below it.

> **This is the sum-product algorithm.** Felsenstein pruning (1973/1981), the HMM forward
> algorithm, and variable elimination in graphical models are the same algorithm, discovered
> three times. If HMMs are familiar: **the forward algorithm generalized from a chain to a tree.**

**Partial likelihood** — a vector of length $|A_i|\le N+1$ per node:

$$L_i(a)=P\big(\text{data at all tips below } i \mid \text{node } i \text{ in state } a\big)$$

$$\textbf{base (tips):}\quad L_m(a)=\mathbf{1}[a=v^m]$$

> ⚑ This is an **emission probability** — $P(\text{observation}\mid\text{true state})$ —
> degenerate (0/1) because observation is assumed perfect. **Replace this single indicator with
> a proper $P(\text{observed}\mid\text{true})$ and you have a dropout/error model. That line is
> where A6 and A7 enter.**

$$\textbf{recursion:}\quad L_i(a)=\underbrace{\Bigg[\sum_{b\in A_k}P_{a,b}(\tau_k)L_k(b)\Bigg]}_{\text{left subtree}}\times\underbrace{\Bigg[\sum_{b\in A_j}P_{a,b}(\tau_j)L_j(b)\Bigg]}_{\text{right subtree}}$$

*Product of brackets* = conditional independence of subtrees given $i$'s state. *Sum inside each
bracket* = marginalize the child's unknown state.

$$\textbf{termination:}\quad \text{Lik}=\sum_{a\in A_\text{root}}P_{\varnothing,a}(\tau_\text{stem})\,L_\text{root}(a)$$

One post-order traversal — and per §1c.1 the node numbering **is** a valid post-order, so it's a
plain loop.

### Worked example

$N=2$, $M=2$ ($\alpha,\beta$), $f_\alpha=f_\beta=0.5$, $r=1$, all branches length 1.
Tips: $v^1=(\alpha)$, $v^2=(\alpha,\beta)$.

$A_1=\{\varnothing,(\alpha)\}$, $A_2=\{\varnothing,(\alpha),(\alpha\beta)\}$,
$A_\text{root}=A_1\cap A_2=\{\varnothing,(\alpha)\}$.

| Transition | Form | Value |
|---|---|---|
| $P_{\varnothing,\varnothing}$ | $e^{-1}$ | 0.3679 |
| $P_{\varnothing,(\alpha)}$ | $1\cdot e^{-1}\cdot f_\alpha$ | 0.1839 |
| $P_{(\alpha),(\alpha)}$ | $e^{-1}$ | 0.3679 |
| $P_{(\alpha),(\alpha\beta)}$ | $[1-e^{-1}]\cdot f_\beta$ | 0.3161 |
| $P_{\varnothing,(\alpha\beta)}$ | $[1-e^{-1}(1+1)]\cdot f_\alpha f_\beta$ | 0.0661 |

($|b|=2=N$ triggers the saturating branch of Eq. 5.)

$L_\text{root}(\varnothing)=0.1839\times0.0661=0.01215$
$L_\text{root}((\alpha))=0.3679\times0.3161=0.11627$

$$\text{Lik}=0.3679\times0.01215+0.1839\times0.11627=0.00447+0.02139=\mathbf{0.02586}$$

Brute-force Eq. 11 cross-check: $a_\text{root}=\varnothing$ gives
$0.3679\times0.1839\times0.0661=0.00447$; $a_\text{root}=(\alpha)$ gives
$0.1839\times0.3679\times0.3161=0.02139$. **Identical ✓**

Posterior weight on the root state is already visible: $0.02139/0.02586=83\%$ on $(\alpha)$ —
sensible, since both tips share $\alpha$ at site 1.

*(At $n=2$ pruning saves nothing — one internal node. The saving is
$(N+1)^{n-1}\to n(N+1)^2$, so it only bites as $n$ grows.)*

## 1c.5 Complexity and practicalities

$\le(N+1)^2$ transition probabilities per node per tape ⇒

$$O(\underbrace{knN}_{\text{Eqs. 8–9}}+\underbrace{knN^2}_{\text{pruning}})=O(knN^2)$$

**Pruning converts exponential in $n$ into linear in $n$.** That is the entire content of Eq. 12.

*Structural caveat:* this works because the phylogeny is a **tree**. On a network with
reticulation, variable elimination gains a treewidth-dependent cost. Cell lineages are genuine
trees (no recombination, no cell fusion), so unlike viral/bacterial phylogenetics this is free
here.

**Two practicalities the paper skips:**

- **Underflow.** $L_i(a)$ is a product of many small probabilities; at $n=1000$ it underflows
  double precision well before the root. Needs per-node scaling factors or log-space arithmetic
  — inherited from BEAST 2's `GenericTreeLikelihood`.
- **Caching.** *"the likelihoods for subtrees are saved at each internal node."* Combined with
  §1c.2 ($A_i$ is topology-only), the caching boundary is clean: **moves on $\lambda$, $\xi$, branch
  lengths reuse $A_i$; topology moves invalidate both $A_i$ and $L_i$ above the change.**
  3.5× on a laptop, 8× overall with tape-level parallelism.

## 1c.5½ ⚑ What pruning does NOT give you — and why it matters for ENGRAM

Pruning returns $P(D\mid\theta)$ with ancestral states **marginalized out**. It does not tell you
what those states *were*. Recovering them needs a second root-to-tip pass — the backward half of
forward–backward, or stochastic mapping.

**"Inference of signalling history" is an ancestral-state problem**, and which kind determines
how much machinery you must build:

| | Global signal trajectory | Lineage-specific signal |
|---|---|---|
| Model | one shared $\xi(t)$, a few knots | each branch has its own signal exposure |
| Latent dimension | $O(1)$ parameters | $O(n)$ latent variables |
| Pruning suffices? | **Yes** — $\xi(t)$ is just a parameter; Eqs. 8–12 unchanged | **No** — needs per-branch latents + reconstruction pass |
| Biological question | *"when did the embryo experience Wnt?"* | *"which lineages experienced Wnt?"* |

The first drops straight in via the $H_m$ recursion (§1c.7). The second is what you probably want
biologically, and is substantially harder.

> **Structural analogue already solved in this literature: branch-specific $\xi$ is to SciPhy
> what a relaxed clock is to standard phylogenetics** — branch-specific rates drawn from a shared
> distribution, sampled by MCMC alongside the tree. SciPhy already cites the uncorrelated relaxed
> clock (Drummond et al. 2006, ref. 23) and **random local clocks** (Drummond & Suchard, ref. 24).
> Borrow rather than invent — the random-local-clock formulation (signal switches on at a node and
> persists through a clade) maps unusually well onto how developmental signalling behaves.

## 1c.6 ⚠ What this costs — three liabilities

1. **The intersection is brittle.** $A_i$ is a **hard** constraint from exact sequence
   identity. One miscalled symbol at site 1 in one cell forces $\text{lcp}=\varnothing$ for
   every node above it, annihilating all ancestral information in that clade. Distance methods
   degrade gracefully; this does not. Sequential error (**A7**) is ~0.6% for TAPE-1 but enters
   **multiplicatively across nodes**. This is why the mouse de-noising is so elaborate.
2. **Missing data handled by deletion, not modeling.** Validation *"generate[s] input
   alignments by constructing the largest complete set of tapes without missing data that
   includes at least a third of all cells."* Subset to a complete submatrix, discard the rest.
   With ~50% of the mouse matrix missing this is **the binding constraint** on applying SciPhy
   there. Row **A6**; acknowledged: *"Follow-up work should integrate a dropout model."*
3. **Root placement is load-bearing.** No reversibility ⇒ no pulley principle ⇒
   $A_{2n}=\{\varnothing\}$ is what makes the likelihood well-defined, not a convenience.

## 1c.7 ⚑ ENGRAM: the extension appears to stay analytically tractable

**Felsenstein pruning does not require time-homogeneity** — only the Markov property and
conditional independence of subtrees given the node state. So with time-varying $f$, the
*entire algorithm above survives unchanged*: state sets, intersection, recursion, $O(knN^2)$.
**Only $P_{a,b}$ changes**, from a function of $\Delta t$ to a function of $(t,t+\Delta t)$.
The tree is time-calibrated, so every node already carries an absolute date.

Write the joint density over event times and symbols directly ($j=|c|$, $T=t+\Delta t$):

$$P_{a,b}(t,T)=\int_{t<u_1<\cdots<u_j<T} w(u_j)\; r^j\prod_{m=1}^{j}f_{c_m}(u_m)\;du_1\cdots du_j$$

$$w(u_j)=\begin{cases}e^{-r(T-t)} & \text{non-saturating } (|b|<N)\\ e^{-r(u_j-t)} & \text{saturating } (|b|=N)\end{cases}$$

The saturating weight drops $e^{-r(T-u_j)}$ because once the tape is full, further Cas9
activity is irrelevant.

**Sanity check — constant $f$:**
- *Non-saturating:* simplex volume $=(\Delta t)^j/j!$ ⇒ $\frac{(r\Delta t)^j}{j!}e^{-r\Delta t}\prod_m f_{c_m}$ = **Eq. 3 × Eq. 6** ✓
- *Saturating:* integral $=P(T_j\le\Delta t)=1-\sum_{m<j}\text{Pois}(m;r\Delta t)$ ⇒ **Eq. 4 × Eq. 6** ✓

**Computes by nested 1-D quadratures:**

$$H_0(u)=1,\qquad H_m(u)=\int_t^u f_{c_m}(s)\,H_{m-1}(s)\,ds$$

$$P_{a,b}=\begin{cases} r^j e^{-r\Delta t}\,H_j(T) & \text{non-saturating}\\[1ex] r^j\displaystyle\int_t^T f_{c_j}(u)\,e^{-r(u-t)}\,H_{j-1}(u)\,du & \text{saturating}\end{cases}$$

At most $N=6$ nested 1-D integrals, each cheap; tabulable once per $(t,T)$ grid rather than
per branch.

> **Cost of the ENGRAM extension: $O(knN^2\cdot Q)$ — still linear in the number of cells.**
> The scalability argument for abandoning likelihood-based inference does **not** obviously
> apply to time-varying $f$.

**⇒ Sharpens the SBI case.** Reserve SBI for what genuinely resists a written-down
likelihood — correlated dropout (**A9**), per-lineage rate multipliers from PE silencing
(**A8**), cell-state-dependent rates (**A2**), within-cell frailty coupling tapes (**A5**).
Those are where the likelihood does not factorize at all. The signal-history part looks like
it stays in reach.

*To verify:* (i) numerical check of the recursion against simulation; (ii) identifiability —
can $\xi(t)$ actually be recovered, given only $N=6$ ordered positions per tape? (iii) does
$\lambda$ also need to vary, and does the factorization survive if it does?

---

## Open empirical question: is there a *cis*-preference of a pegRNA for its own tape?

**Question.** Since pegRNA and TAPE are co-integrated on one construct, does a pegRNA
preferentially edit the tape it sits next to, by proximity?

**Status: not reported in either paper.** Closest panel is Typewriter **ED Fig. 6c**
("Read counts of InsertBCs observed in TAPE-1 arrays") — a rank-abundance plot of 19
InsertBCs *pooled across all tapes*, i.e. a **marginal**, not the TargetBC × InsertBC
contingency table. No such analysis in main text, extended data, or methods of either
paper. *(Not checked: the ENGRAM paper, mouse supplementary tables.)*

**Testability differs by dataset:**

- **Typewriter lineage data — not testable.** They identify 19 TargetBCs and 19
  InsertBCs but **never link which InsertBC is encoded on which integration** (would
  require sequencing the intact provirus). The diagonal is undefined.
- **Mouse data — testable.** *"we performed amplicon sequencing of the full
  `epegRNA::TAPE-BC::circTAPE` integrations, which recovered **10 of the 11** TAPE
  integrations and showed them to collectively encode all 8 insertions present at ≥0.5%
  at each of the six sites."* Pairing known for 10/11. Data public: GEO **GSE341627**,
  `github.com/seidels/dtt-mouse-analysis`.

### Prior: expect essentially zero, on four grounds

1. **Timescale separation.** Cas9 RNP nucleoplasmic $D \sim 0.1$–$3\ \mu\text{m}^2\,
   \text{s}^{-1}$; nuclear radius ~5 µm ⇒ mixing time $t \sim L^2/6D$ = seconds to tens
   of seconds. Editing runs at ~0.8 edits/tape/day (mouse peak) or ~0.12 (Typewriter).
   **4–6 orders of magnitude** of separation ⇒ complete mixing.
2. **Low per-encounter efficiency.** Cas9 target search involves extensive futile PAM
   interrogation; prime editing additionally needs PBS annealing + RT extension + flap
   resolution (hence single-digit % per-site efficiency). A productive RNP has sampled
   the genome enormously many times; positional bias is washed out.
3. **Steady state, not pulse.** Constitutive U6 expression ⇒ stationary, near-uniform
   nuclear pegRNA distribution. Cis-bias would need consumption to outrun diffusion —
   the opposite of the actual regime.
4. **Decisive: the pegRNA is inert until it loads a PE, and PE is expressed from a
   *different* locus** (separate piggyBac integration; ~7 copies in mouse embryo #3).
   The pegRNA **must** diffuse away from its birthplace to meet a prime editor. By the
   time the RNP is competent to target, it has no memory of where its guide was made.
   **Locality is destroyed before targeting begins.**

*Contrast:* RNAs that genuinely act in *cis* (canonically **Xist**) require dedicated
retention machinery (A-repeat, SPEN/SHARP) precisely to resist diffusive
homogenization. Pol III sgRNAs have none. 1D sliding is short-range (tens–hundreds of
bp) while the pegRNA cassette sits kilobases from the TAPE across vector backbone.

### The shared-cause confound — and why it does *not* break the test

Integration $i$'s chromatin openness $\theta_i$ raises pegRNA$_i$ abundance (lifting
$\xi_i$ **globally**, on every tape) *and* raises tape $i$'s accessibility (lifting
$r_i$). Since the cassette and TAPE are adjacent on one construct, this coupling should
be strong. **But** under multiplicative independence,

$$\mathbb{E}[\text{edits of symbol } i \text{ on tape } j] \;\propto\; \xi_i \times \lambda_j$$

so $\mathbb{E}[\text{symbol } i \text{ on tape } i] \propto \xi_i \lambda_i$ is elevated
*exactly as independence predicts*. The shared cause inflates the **marginals**, not the
interaction. A test for **diagonal enrichment conditional on both marginals**
(quasi-independence / log-linear model with a diagonal indicator; row sums = per-tape
edit totals, column sums = per-symbol totals) absorbs it. Same for copy number. **The
test is clean.**

Two things that *do* bite:

- **Phylogenetic non-independence.** An edit inherited by a 50,000-cell clade is **one**
  event, not 50,000. A naive cell-level table is catastrophically anticonservative.
  Count each edit once at the branch where it first appears — the same reconstruction
  the mouse paper already implements for its edit-count branch-support statistic.
- **Symbol degeneracy.** 10 integrations encode 8 distinct symbols ⇒ "the diagonal" is
  the set of (tape, symbol) pairs where that tape's construct encodes that symbol.

Power: ~55 edits/cell × 1.34M cells; even collapsed to independent events on ~10⁶
internal branches, a few-percent effect should be detectable.

### Why it matters for the modeling

| Outcome | Consequence |
|---|---|
| **No cis-preference** (predicted) | Clean empirical license for SciPhy's shared-$f$-across-tapes assumption — currently asserted, never tested |
| **Cis-preference exists** | $\xi$ becomes tape-indexed $\xi_i^{(z)}$ (parameter count × #tapes) and introduces **structured** diagonal enrichment any simulator must reproduce. Better to know now than to meet it as unexplained residual structure in posterior predictive checks |

*Action:* ask Sophie Seidel (author on both papers) whether this has been looked at
before spending time on it.

---

## Open methodological question: transferring SBI calibration across experiments

**The worry (correct).** SBI is markedly less robust to model misspecification than
likelihood-based inference: a density estimator trained on simulations can be
*confidently* wrong on data off the simulated support. If chromatin context differs
between the calibration dataset and a new experiment, calibrated nuisance parameters
will not transfer.

**The reframe.** Reject "calibrate on X, transfer to Y." Chromatin context is **fixed
within an experiment** — integrations happen once, in the zygote or founder cell, and
are clonally inherited by every assayed cell. So the per-tape rate is *not* a universal
biological constant to be calibrated once; it is an **experiment-specific nuisance
parameter inferred from the data at hand.** SciPhy already gets this right by inferring
$r_z$ per tape. Preserve it.

What transfers is **structure + hyperpriors**, not point estimates:

1. **Amortize over the nuisance.** Train the posterior estimator on simulations with
   $r_z$, dropout, and $f$ drawn from broad priors spanning plausible chromatin effects.
   Then you have marginalized over contexts rather than calibrated to one. The failure
   mode bites only if the real experiment falls outside the prior's support — a
   *checkable* condition, not an invisible one.
2. **Go hierarchical.** Estimate the *distribution* of per-tape rates across many
   experiments (e.g. lognormal with fitted spread) and use it as a hyperprior. That
   object genuinely transfers; individual $r_z$ values do not.
3. **Keep the validation machinery.** SciPhy already runs simulation-based calibration —
   "well-calibrated simulations," coverage of the 95% HPD, citing Dawid (1982) and
   Mendes et al. (2024). These checks port to SBI unchanged. Add prior-predictive checks
   against held-out real data to test support coverage.

**The honest tradeoff.** SBI buys model richness (dropout, frailty, cell-type-dependent
rates, non-clock recording — all trivial to simulate, all intractable to write down) at
the cost of transferability. So spend that richness on making the nuisance structure
**explicit and inferrable**, not on hard-coding calibrated constants. Hard-coding
chromatin effects yields a machine that works on exactly one dataset.

---

---

# Session 2a — Inference at the theoretical level

## 2a.1 MCMC refresher

**The problem:** $p(\theta\mid D)=p(D\mid\theta)p(\theta)/p(D)$ with
$p(D)=\int p(D\mid\theta)p(\theta)d\theta$ — an integral over every topology × every node-time
assignment × every $\lambda,\xi,b,\delta$. Hopeless.

**The escape:** never compute $p(D)$. Build a Markov chain whose stationary distribution *is*
the posterior; evaluate only ratios, in which $p(D)$ cancels.

### Metropolis–Hastings

From current $\theta$: propose $\theta'\sim q(\theta'\mid\theta)$, accept w.p.

$$\alpha=\min\left(1,\ \underbrace{\frac{p(D\mid\theta')}{p(D\mid\theta)}}_{\text{likelihood ratio}}\cdot\underbrace{\frac{p(\theta')}{p(\theta)}}_{\text{prior ratio}}\cdot\underbrace{\frac{q(\theta\mid\theta')}{q(\theta'\mid\theta)}}_{\text{Hastings ratio}}\right)$$

otherwise **stay put** (and record the current state again).

The **Hastings ratio** corrects proposal asymmetry; =1 for symmetric proposals. ⚠ Where most
implementation bugs live — a subtly wrong Jacobian gives a chain that converges *confidently to
the wrong distribution*.

**Why it works:** detailed balance $\pi(\theta)P(\theta\to\theta')=\pi(\theta')P(\theta'\to\theta)$
+ irreducibility + aperiodicity ⇒ convergence to $\pi$ from any start.

### Diagnostics, and what SciPhy's numbers tell you

| Concept | What it is | SciPhy |
|---|---|---|
| **Burn-in** | discard pre-stationary prefix | first **10%** |
| **ESS** | how many *independent* draws the correlated chain is worth | threshold **200**, every parameter |
| **$\hat R$** | Gelman–Rubin: do independent chains agree? ≈1 | **3 chains**, LogCombiner |
| **Chain length** | | **$10^9$ iterations** |

> ⚑ Typical BEAST analyses run $10^7$–$10^8$. **Ten to a hundred times that, for ~1,000 cells, is
> a direct measurement of how badly this posterior mixes.** They even *extrapolate* runtime to
> convergence from achieved ESS — what you do when reaching ESS 200 is impractical.

## 2a.2 What "a distribution over trees" means

### A tree is a pair

$$T=(\underbrace{\tau}_{\text{topology}},\ \underbrace{\mathbf t}_{\text{node times}})$$

- **Topology** — discrete, from a finite set of size $(2n-3)!!$.
- **Node times** — continuous, $\mathbb R^{n-1}_+$, constrained so every parent is older than its
  children.

⇒ The parameter space is **a finite stack of continuous spaces**. A book where each *page* is one
topology and the page's surface is the space of compatible node-time assignments.

```
   page 1: ((A,B),C)     page 2: ((A,C),B)     page 3: ((B,C),A)
   ┌───────────────┐     ┌───────────────┐     ┌───────────────┐
   │   ·  ·   ·    │     │        ·      │     │   ·           │
   │ ·   ···  ·    │     │    ·  ·       │     │        ·      │
   │  ·  ·· ·      │     │      ·        │     │   ·           │
   └───────────────┘     └───────────────┘     └───────────────┘
   each dot = one (topology, node-times) state visited by the chain

   WITHIN a page  = change node times   (continuous move)
   BETWEEN pages  = change topology     (discrete jump)
```

Page counts: $n=4\to15$; $n=10\to34{,}459{,}425$; $n=1000\to\sim10^{2860}$.

### Operationally, it's a bag of trees

You never write the distribution down. MCMC hands you a **sample**; you answer questions by
**counting**. With $n=3$ and 1,000 samples:

| Topology | Count | Posterior probability |
|---|---|---|
| `((A,B),C)` | 700 | **0.70** |
| `((A,C),B)` | 200 | 0.20 |
| `((B,C),A)` | 100 | 0.10 |

"Posterior probability that A and B are sisters = 0.70" is literally the fraction of sampled
trees containing that clade — the Bayesian counterpart of a bootstrap value. **Exactly the
quantity the mouse paper could not compute** (*"resampling-based bootstrap support is
computationally infeasible at this scale"*), hence their edit-count support statistic.

Within those 700, the A–B ancestor's age varies sample to sample; that spread **is** the posterior
for that node's date. Every marginal is a histogram over the bag.

**Summary trees.** MCC (maximum clade credibility) picks the best *sampled* tree — but when tree
space dwarfs the sample, no sampled tree may be good. SciPhy uses the **conditional clade
distribution, CCD0 parameterization** (Berling et al., ref. 19), which builds a distribution from
conditional clade frequencies and can therefore score and select trees that **never appeared in
the sample**. Right tool when the sample is sparse relative to the space.

### ⚠ What the tree prior actually contributes

**Under a constant-rate birth–death-sampling process, all labeled topologies are equally
probable** — the process is exchangeable in the tips. So the tree prior is essentially
**uninformative about topology**.

It is highly informative about **timing** — the joint distribution of node ages given
$b,\delta,\rho$ and the origin. Three jobs, none of them "prefer certain shapes":

1. **Regularizes node times** toward what a branching process would produce.
2. **Makes $b,\delta$ inferrable** — the phylodynamic payoff; the prior is a *scientific model*,
   not a nuisance.
3. **Makes the posterior proper.**

> **Instructive contrast.** The mouse paper hit the same disease (implausibly early node dates) and
> treated it with a hard external constraint — *"a literature-derived embryo cell-count ceiling"*
> halved per blastomere, imposed as an order-statistic ladder. SciPhy gets equivalent
> regularization **endogenously** from the birth–death prior, with population parameters inferred
> rather than looked up. Same disease; the difference in treatment is a fair summary of what the
> Bayesian machinery buys.

## 2a.3 Branch lengths

**They are a co-equal component of the tree parameter, not an aspect of the topology.** The same
topology supports infinitely many time assignments; the tree is the pair, both halves inferred.
And the $\Delta t$ of Eqs. 3–7 **is** the $\tau_i$ of Eqs. 11–12 — same object.

| Component | Role in the likelihood |
|---|---|
| **Topology** | parent–child relations → structure of the pruning recursion, which leaves feed each $A_i$ |
| **Branch lengths** | enter *numerically* via $P_{a,b}(\tau_i)$, always as the product $r\tau_i$ |

### Parameterization is node HEIGHTS, not lengths

BEAST stores node **heights** (ages); $\tau_i=\text{height}(\pi_i)-\text{height}(i)$ is derived.

*Why:* the tree is **ultrametric** (all cells harvested at one instant). With heights that's
automatic — fix tips at 0, constrain parent older than child. With branch lengths as free
parameters you'd have to enforce that **every** root-to-tip path sums to the same total, a coupled
constraint almost no proposal would satisfy.

### Rate–time confounding, and what breaks it

The likelihood sees $\tau_i$ **only** through $r\tau_i$. Scale every branch by $c$, divide $\lambda$ by
$c$ ⇒ identical likelihood. Three things break it, and **your extension inherits all three**:

1. **The origin is fixed** at the known experiment duration (25 days; 11 days) — not inferred, not
   given a prior; *fixed*.
2. **The tree is ultrametric** ⇒ total origin-to-tip time is that same fixed number for every cell.
3. Together ⇒ **you cannot rescale**; tree height is pinned at both ends.

Hence $\lambda$ is identified in genuine edits/tape/day, and the prior reads *"between one and ten edits
per tape over 25 days."* Classical phylogenetics needs fossil calibrations or serially sampled tips
for this; a designed experiment hands it over.

*MCMC consequence:* $\lambda$ and node heights are strongly **correlated** in the posterior (only the
product is locally constrained). Remedy is an **up-down operator** scaling $\lambda$ down while scaling
heights up — moving *along* the ridge instead of across it. Without one, mixing degrades badly.

## 2a.4 What the chain moves, and why it struggles

| Kind | Operator | Effect |
|---|---|---|
| **Between pages** | Narrow exchange (NNI) | swap adjacent subtrees — local, high acceptance |
| | Wide exchange | swap random subtrees — big jump, low acceptance |
| | Subtree slide | detach/reattach nearby, adjusting heights |
| | Wilson–Balding | detach/reattach anywhere — largest move |
| **Within a page** | Uniform node height | move one node between child-max and parent height |
| | Scale / root-height scale | rescale subtree or whole tree |
| **Parameters** | Scale | $r,b,\delta$ |
| | Delta-exchange / Dirichlet | $\xi$ — must stay on the simplex |
| **Joint** | Up-down | $\lambda$ vs. heights, along the confounding ridge |

Matches the paper: *"tree operators modify the tree to explore all possible topologies that could
have produced the data. As some of these changes are local, the likelihood for subtrees that are
not modified remains unchanged. We make use of that property by implementing caching."*

**Why mixing is hard here — four compounding reasons:**

1. **The posterior is genuinely flat over huge regions.** Genotype-identical cells ⇒ thousands of
   topologies with *literally equal* likelihood. Not hill-climbing — **diffusing on a plateau**.
   Diffusion time scales quadratically in distance; the plateau is enormous.
2. **Local moves, distant targets.** NNI changes one bipartition; reaching a very different topology
   needs many consecutive accepted moves, and any intervening dip stalls you.
3. **Rate–time correlation** — a diagonal ridge that axis-aligned proposals cross rather than follow.
4. **Dimension grows with $n$** — $n-1$ heights plus topology; each sweep perturbs a vanishing
   fraction.

> ⚑⚑ **The $O(knN^2)$ likelihood is NOT the bottleneck; the random walk over tree space is.**
> Faster pruning buys more steps of a chain that needs exponentially many. Leverage is in items
> 1–3: better proposals, or a formulation that doesn't visit topologies one at a time.
>
> **This reframes the SBI question.** "Amortize the editing-model likelihood" attacks nothing.
> "Amortize over tree space" — learn a map from data directly to a posterior over trees, or to a
> distribution you can sample *independently* rather than by random walk — attacks the actual
> problem. Very different research programs; the second is far more ambitious.

---

# Comparison: three tree-building approaches on DNA Typewriter data

Third paper added: **Park, Chang et al. 2026**, *"Time-resolved lineage recording reveals a
pre-existing, heritable cell state underlying metastatic potential"* (bioRxiv
2026.08.10.744013) — NCI-H1299 lung cancer, orthotopic xenograft, 3 mice, 105,260 cells.

| | **Mouse embryo** (Yang/Seidel) | **Metastasis** (Park/Choi) | **SciPhy** |
|---|---|---|---|
| Tips | 1.34M | ~74 cells/clone × 75 clones; subclones 313–10,506 | ≲1,000 |
| Recorder | 11 tapes × 6 = 66 sites | **166 tapes × 6 = 996 sites** | 13×5 = 65 (HEK293T) |
| Tree search | NJ on typewriter distance + nearest-anchor placement | NJ start → **NNI hill climb** + guided SPR | **MCMC over topology + times** |
| Criterion | distance | **sequential (prefix) parsimony** (explicit Camin–Sokal generalization) + $\lambda\cdot$Fitch(organ) | **posterior under the mechanistic CTMC** |
| Ancestral states | most-parsimonious | *"each internal node is labeled by the longest common prefix of its non-missing descendants"* | **integrated over $A_i$** |
| Clock | strict, LSD2 | **saturating** $E(t)=A(1-e^{-t/\tau})$, calibrated externally, held fixed | strict |
| Date regularization | **external cell-count ceiling** | **external in-vitro calibration** | **birth–death prior, inferred** |
| Uncertainty | edit-count branch support (bootstrap *"computationally infeasible at this scale"*) | tape bootstrap, Kishino–Hasegawa test, $\lambda$-frontier, model-free distance-ratio | **full posterior** |
| Population params | none | per-clone $\lambda$ by OLS on in vitro frequencies | **$b,\delta$ from tree shape** |

> ⚑ **Park et al.'s ancestral reconstruction is literally SciPhy's $u_i=\mathrm{lcp}$.** They
> independently derived the same combinatorial object (§1c.2) and then took its **maximum
> element** as a point estimate rather than integrating over the set. The parsimony/likelihood
> distinction, now confirmed in two independent papers.

## What SciPhy buys — read off what the others had to bolt on

1. **Both import their clock from outside the data.** Mouse: literature-derived cell-count
   ceiling as an order-statistic ladder. Park: $A,\tau$ calibrated from *separate* in vitro time
   courses, held fixed — and the headline result (liver seeding at day 21 / day 34) rests
   entirely on it. Bayesian infers the clock from the same data and propagates uncertainty into
   the dissemination date.
2. **Both build bespoke per-analysis uncertainty machinery.** Park runs four separate devices
   (tape bootstrap, KH test, $\lambda$-ladder in three passes, distance-ratio cross-check). A
   posterior answers all uniformly.
3. **Neither handles homoplasy properly.** Park argues prefix parsimony is *"insensitive to
   homoplasy at the later sites"* — correct, since spurious lcp extension to depth $k$ needs
   matches at sites $1..k$, decaying like $(\sum \xi_i^2)^k$. **But site-1 homoplasy is untouched
   by that argument, and site-1 errors corrupt deep structure** (§1c.2 failure mode 1).
4. **The signalling extension has no parsimony analogue.** Inferring $\xi(t)$ is a latent-variable
   problem; there is no parsimony criterion for it.

## The scale question — sharpened

**The target "scale Bayesian inference to 1.3M cells" is wrong, for three reasons.**

**1. The plateau does not yield to compute.** 0.6 private edits/cell ⇒ the posterior genuinely
*is* flat over enormous regions. Scaling MCMC there spends vast compute characterizing a
distribution whose honest content is "I don't know." An **information** problem, not a compute
problem.

**2. ⚑ The complementary lever is recording capacity — and Park et al. just demonstrated it.**
996 sites vs. 66. ~254 edits/cell by day 7, several hundred by harvest, ~150,000 unique
insertion patterns, ~15 bits entropy, *"recording diversity grew rather than plateaued."* With
NNNNGGA, $M\approx256$ ⇒ $\sum \xi_i^2\sim0.4\%$ — **near homoplasy-free**.  ⚠ **Wrong — measured $q=0.0170$, $1/q\approx58$; see the correction in §D.4b.**

> **Better recorders don't just give better trees — they make Bayesian inference tractable by
> removing the flatness.** Sharp posterior ⇒ MCMC has a gradient to climb instead of a plateau to
> diffuse across. Note large $M$ is **free** in SciPhy: $|\Omega|$ never enters the computation,
> only $|A_i|\le N+1$.

**3. Subsampling is statistically principled and underused.** The birth–death-sampling prior
models $\rho$-sampling explicitly; SciPhy fixes $\rho=0.0008$ for HEK293T — i.e. it *already*
correctly accounts for analysing 1,000 of 1.2M cells.

| Question | Cells needed |
|---|---|
| Population dynamics ($b,\delta$, growth phases) | correct **random subsample** — $\rho$ handles it |
| **Global signal trajectory $\xi(t)$** | same — **the tree is a nuisance parameter** |
| Editing-model parameters ($\lambda,\xi$, dropout) | same |
| *Which* lineage founded this metastasis | **that specific cell** — no subsampling |
| Fine-grained coupling between rare cell types | the rare cells, at minimum |

⇒ **For the ENGRAM goal, most value sits in the top three rows, where the scale problem
substantially dissolves.**

## Where SBI actually helps — three distinct programs

**Program A — marginalize the tree, infer the parameters. Achievable, high value.**
Simulate birth–death tree → editing → dropout → cell×tape×site matrix. *Learn* summaries via a
permutation-invariant embedding over cells. Train NPE for
$(b,\delta,r,\xi,\xi(t)\text{-knots})$. The tree is never inferred — marginalized by
simulation. **Handles exactly the A-table rows where the exact likelihood fails: A9, A8, A2,
A7.** Well-posed; startable now.

**Program B — amortize over tree space itself. Hard research program.**
Relevant prior art: variational Bayesian phylogenetic inference with **subsplit Bayesian
networks** (Zhang & Matsen). ⚑ **SBNs and the conditional clade distributions SciPhy already uses
for CCD0 are the same family of object** — both factorize a topology distribution via conditional
clade probabilities. SciPhy uses it post hoc for summarization; VBPI uses it as the variational
family. Shorter bridge than it looks.

**Program C — ML tree search under the SciPhy likelihood. The missing rung; nobody has taken it.**
The mouse authors' own wish was *"an eventual transition to **maximum likelihood** methods that
integrate our prior understanding of the molecular biology of sequential genome editing"* — ML,
not Bayesian. ML tree search scales to $10^5$–$10^6$ tips routinely.

> **Park et al. have already built 90% of the machinery** — NJ start, NNI hill climb, guided SPR,
> penalized objective, bootstrap for comparing adjacent solutions. **Swapping
> `sequential_parsimony(T)` for `−log Lik(T)` from the pruning recursion is close to a drop-in.**
> The likelihood is $O(knN^2)$ — cheap. You get the mechanistic model, correct homoplasy
> treatment, jointly estimated $\lambda$, at ~100× the scale.
>
> The field currently jumps from parsimony straight to full Bayesian MCMC, skipping the rung that
> would work today.

## Recommended sequence

1. **Immediate — run SciPhy as-is on Park et al.'s per-clone trees.** ~74 cells/clone, 75 clones,
   $k=166$, $N=6$, $M\approx256$, near-zero homoplasy, tapes far from saturating. Squarely in
   range, embarrassingly parallel across clones. Replaces *"day 21/34 with a bootstrap CI
   conditional on a fixed external clock"* with genuine posterior dissemination dates + jointly
   inferred editing rates. Cleanest possible test of whether the Bayesian treatment changes
   conclusions — or doesn't, which is itself worth knowing.
2. **Then — Program C** (ML tree search under the SciPhy likelihood): the scalability answer for
   topology.
3. **In parallel — Program A** for the ENGRAM extension, where the tree is a nuisance and the
   misspecification terms are the whole difficulty.

> **Framing to resist:** "Bayesian inference is more precise, so scale it to everything."
> **Sharper:** Bayesian inference is worth its cost **specifically where uncertainty propagates
> into the scientific claim** — dating, migration counts, signal trajectories, population
> parameters — and those questions mostly *don't* need 1.3M cells. Where you need every cell,
> you're usually asking a question a point estimate answers adequately.

---

# Design notes: MAP tree search, and the tape-independence question

## D.1 MAP vs. ML — definitions

$$\hat\theta_{\text{ML}}=\arg\max_\theta p(D\mid\theta)\qquad\qquad \hat\theta_{\text{MAP}}=\arg\max_\theta \underbrace{p(D\mid\theta)p(\theta)}_{\propto\,p(\theta\mid D)}$$

**MAP = Maximum A Posteriori = the *mode* of the posterior.** ⚠ **Not the posterior mean** — the
mean is a different estimator (Bayes estimator under squared-error loss; MAP is under 0–1 loss).
They coincide only for symmetric unimodal posteriors.

> **The useful framing: MAP = penalized maximum likelihood.** $\log p(\theta)$ *is* the penalty.
> Gaussian prior ↔ L2; Laplace prior ↔ L1. Here the "penalty" is $\log p(T\mid b,\delta,\rho)$ —
> the birth–death density on node times, which is what keeps node ages plausible. Pure ML has no
> such term and drifts into the mouse paper's pathology.

⚠ **Parameterization invariance.** The MLE **is** invariant under bijective reparameterization;
**the MAP is not**, because densities pick up a Jacobian that shifts the mode. This bites directly
— the node-height ratio transform is exactly such a reparameterization, so the MAP in ratio space
≠ image of the MAP in height space. Optimize in the unconstrained space, add $\log|J|$ if you want
the mode in original coordinates, and state which you report.

**Gradients:** yes for **continuous** parameters (node heights, $\lambda$, $\xi$, $b$, $\delta$,
$\xi(t)$ knots) — autodiff through pruning, then L-BFGS/Adam. **No** for topology — discrete,
searched combinatorially. Hybrid: discrete search outside, gradient ascent inside.

## D.2 Assessment of the direction

**Why it's promising.** Rigorous + scalable + extensible-to-signalling are jointly served by
MAP-under-a-mechanistic-likelihood and by nothing else on offer. Parsimony: scalable, not rigorous,
**no parsimony criterion exists for $\xi(t)$**. MCMC: rigorous and extensible, not scalable. And the
no-matrix-exponential gradient advantage is real and unusual.

**Caveat 1 — the plateau changes character, it doesn't vanish.** Where topologies are near-equally
likely, MCMC diffuses forever; hill climbing stops arbitrarily **and reports it confidently** —
arguably worse. Build bootstrap/SH-aLRT support in from the start, not as an afterthought.

**Caveat 2 — framing.** "MAP instead of MCMC" is incremental. "Signalling history inferred from
lineage recordings at scale, which required a new inference framework" is new. **Let the ENGRAM
application drive the project; MAP is the enabling technique, not the goal.**

## D.3 The tape-independence question

All three papers assume **conditional independence of tapes given the tree**:
SciPhy multiplies per-tape densities; Park sums prefix-parsimony *"over tapes"*; the mouse pools
shared edits into one distance.

### ⚑ Key correction: factorization is NOT information loss

$$P(D\mid T,\theta)=\prod_{z=1}^{k}P(D_z\mid T,\theta)$$

If tapes are conditionally independent given the tree, this **is the exact joint likelihood** — not
an approximation. Every symbol on every tape enters; the likelihood is a sufficient statistic for
the parameters. **There is no extra information in "tapes considered jointly" that the product
misses, because the product *is* the joint distribution.** Any other combination rule could only
lose information.

### ⚠ What the assumption actually costs: overconfidence

Tapes in one cell share the PE pool, pegRNA pool, cell-cycle phase, dNTP supply, global chromatin
state (**A5**, **A8**). Multiplying likelihoods of *correlated* observations as if independent
doesn't lose information — **it overstates the evidence.** Classic pseudo-replication: intervals
shrink too fast, clade support looks stronger than it is. With $k=166$, treating 166 correlated
observations as independent could inflate confidence a lot.

**The fix is a shared latent, not a different combination rule:**

$$P(D\mid T,\theta)=\int\Big[\prod_{z=1}^{k}P(D_z\mid T,\theta,\epsilon)\Big]p(\epsilon)\,d\epsilon$$

Conditionally independent given $\epsilon$; marginally dependent after integration. That is row A8
done properly.

### Where joint treatment genuinely adds

**1. ⚑ Cross-tape correlation IS the signal for per-lineage rate variation.** An independence model
literally cannot see it: a cell whose tapes are *all* unusually far along reads as 166 separate
coincidences rather than one fast lineage. With $k=166$ this is **strongly powered** (166
near-replicate observations of one per-cell rate); with the mouse's $k=11$ it is noisy.

> **Same structure as the ENGRAM problem.** If signalling drives $f$, then *correlated deviation in
> symbol composition across a cell's tapes* is the signal. Per-tape-independent modelling averages
> it away. **Your extension requires this machinery — build it rather than bolt it on.**

**2. Per-tape likelihood contributions are a free diagnostic.** A tape systematically disagreeing
with the consensus flags copy-number ambiguity, a silenced integration, or a genotyping artifact —
exactly the mouse's 2-of-11 tapes recovering at 20% and 3%. Free once you compute the product;
invisible to pooled distance methods.

**3. Batching.** $k$ is the GPU batch dimension. Statistically neutral, practically large.

## D.4 ⚑⚑ The idea to chase: irreversibility ⇒ near **perfect phylogeny**

Treat each edit event (tape $z$, position $i$, symbol $\gamma$) as a **character** whose clade is
the cell set carrying it. Two characters are **compatible** if their cell sets are nested or
disjoint. Classical perfect-phylogeny theory: with irreversible characters and no homoplasy, all
pairs compatible ⇒ the clades form a **laminar family**, and **that laminar family IS the tree** —
readable off in $O(kn)$ by hashing, **no search at all**.

Two DNA-Typewriter-specific features make this nearly true:

1. **Within a tape, compatibility is automatic.** The characters "has prefix $p$" are nested by
   construction (length-3 prefix ⊂ its length-2 prefix). Each tape's characters form a **trie**,
   laminar for free. ⇒ **Incompatibility can only arise BETWEEN tapes.** Large reduction in what
   must be checked. *(Specific to sequential recorders — this is the part that makes it cheap, and
   I'm not aware of anyone exploiting it.)*
2. **Cross-tape incompatibility is rare when $M$ is large.** Arises from homoplasy at rate
   $\approx\sum_i \xi_i^2$. Park: per-site Shannon entropy ~6 bits at early sites ⇒ effective alphabet
   ~64 ⇒ low-single-digit-% collisions. Mouse ($M=8$): $\ge12.5\%$ — much worse.

```
1. Build k tries (one per tape)                          O(kn)
2. Find maximal cross-tape-compatible character set      (NP-hard in general;
                                                          easy when conflicts sparse)
3. Read off the perfect-phylogeny skeleton               O(kn)   ← the tree, mostly done
4. Resolve remaining conflicts LOCALLY by likelihood     small, bounded regions
5. Optimize node heights + r + η + f(t) by gradient      differentiable, GPU
```

**Why it could matter:** the current pipelines' costs are $O(n^2)$ distance computation (NJ at
$n=10^5$ is $10^{10}$ ops — the mouse needed DecentTree precisely for this) plus generic NNI/SPR
over an unstructured space. This replaces both with a near-linear construction that **uses all tapes
jointly by construction**, reserving likelihood search for the small subset where characters
actually conflict.

> **The real answer to "can joint treatment improve search efficiency":** not by changing how the
> likelihood combines information (already optimal) — **but by exploiting the laminar structure
> irreversibility imposes, to avoid searching most of tree space at all.**

**Caveats to check first:**
- Missing data makes perfect phylogeny NP-hard in general; needs heuristics for partial characters,
  and ~50% recovery is not benign.
- **Measure the actual compatibility fraction on Park's data before betting on it** — a
  one-afternoon computation that settles whether the idea is live.
- Cassiopeia's greedy solver is a distant cousin (splits on the most frequent character) but does
  not exploit the ordered/prefix structure.

---

# Design notes II: representing conditional-CRE activity

## E.1 ⚑ ENGRAM writes COMPOSITION, not rate

From Session 0: $[\text{writer for symbol }i]\approx[\text{PE}]\times\frac{\text{pegRNA}_i}{\sum_j \text{pegRNA}_j}$.
**PE dosage sets the rate; pegRNA composition sets which symbol.** PE is almost certainly limiting
(Pol II/III transcripts abundant, Cas9-fusion protein scarce), so CRE activity moves the second
factor only:

$$\lambda(t,\ell)\approx \lambda_0 \qquad\qquad \xi_i(t,\ell)=\frac{\beta_i\,a_i(t,\ell)}{\sum_j\beta_j\,a_j(t,\ell)}$$

$a_i$ = CRE activity (the target); $\beta_i$ = intrinsic per-symbol editing efficiency (Session 0's
"edit score," 10-fold sequence bias) — **calibrate once from a constitutive experiment and fix**.

⇒ Sharpen "editing **rate** over time and lineage" to "**symbol composition** over time and
lineage." Two consequences:

1. **$f$ lives on a simplex ⇒ only *relative* activity is recoverable.** Scaling all $a_i$ by a
   constant leaves $f$ unchanged; a uniform doubling of CRE activity is invisible.
   **⚠ A constitutively-driven reference symbol is a hard requirement**, not a nicety — then
   $\xi_i/\xi_0\propto\beta_i a_i/\beta_0 a_0$ recovers $a_i$ up to a fixed constant.
2. **Experimental design lever:** run *pegRNA*-limiting instead (low CRE output, abundant PE) and
   total rate scales with total pegRNA ⇒ absolute activity via the rate **and** relative activity
   via composition. Strictly more identifiable, at the cost of a slower, noisier recorder.
   **Decide this deliberately; it determines what is estimable.**

**Three drivers are not the same problem:**

| Driver | Latent structure | Key feature |
|---|---|---|
| **Signal** (Wnt, NF-κB) | dynamic trajectory over time and lineage | fully latent |
| **Cell type** | function of cell state | **cell state is OBSERVED at the tips** (scRNA-seq) |
| **TF activity** | similar | partly inferrable from expression |

## E.2 What resolution is achievable

Each edit = **one categorical draw** from $f$ at its latent write time. The count that matters is
**independent events in the tree**, not observed edits (an edit inherited by a 50,000-cell clade is
one draw).

$$\#\text{events}\approx k\cdot r_\text{cell}\cdot\underbrace{\textstyle\int L(t)dt}_{\text{total tree length}},\qquad \text{distributed}\propto L(t)$$

Rates: **Park ≈ 12 edits/cell/day** (166 tapes, ~500–600 edits/45 d, decelerating);
**mouse ≈ 4.1 edits/cell/day**.

1. **⚠ Temporal resolution improves ~exponentially toward the present**, since $L(t)$ grows
   exponentially. Near the root the limit is the *number of events*, not observation noise.
   **Early developmental signals are the hardest case** — awkward, since that's where embryologists
   want them. Metastasis is a good match (dissemination days 21–34 of 45, many lineages alive).
2. **Per-branch estimation is marginal even at $k=166$.** A branch of duration $\tau$ carries
   ~$12\tau$ events (Park) / ~$4\tau$ (mouse). At $\tau=2$ d: 24 vs 8 events, spread over $M=8$
   symbols ⇒ ~3/symbol vs ~1/symbol. **Park: noisy but nonzero. Mouse: hopeless.**
   ⇒ **$k$ is the lever — the quantitative case for high-capacity recorders in signal recording.**
3. ⇒ **Natural resolution is per-clade / per-regime, not per-branch.** Statistical reason, not just
   computational.

## E.3 The shared latent — structure, and where it bites

$$\log a_i(t,\ell)=\underbrace{\mu_i}_{\text{baseline}}+\underbrace{g_i(t)}_{\text{global trajectory}}+\underbrace{u_i(\ell)}_{\text{lineage effect}}$$

$g_i(t)$ = spline, few knots — the "global trajectory" case (§1c.5½), needs no new machinery. The
design decision is $u_i(\ell)$:

| | Structure | Fit to biology | Fit to gradient inference |
|---|---|---|---|
| (a) Switch-point / random local clock | constant within clades, sparse switches | **good** (signalling is on/off) | poor — discrete, RJMCMC/spike-slab |
| (b) **OU / Brownian on the tree** | continuous diffusion along branches | fair (smooth, not switch-like) | **excellent** |
| (c) Discrete regime CTMC on tree | $R$ regimes, transitions along branches | **best** | poor — see below |

### ⚠ Why (c) fails — the general statement of the coupling problem

Hoped-for: expand the pruning state space to $|A_i|\times R\le 7\times5=35$, pay $9\times$. **It
doesn't factorize.** Conditional on the regime tapes are independent, but pruning must marginalize
the regime at each node, and

$$\sum_{s'}\prod_z g_z(a^{(z)},s')\ \neq\ \prod_z\sum_{s'}g_z(a^{(z)},s')$$

The sum over the **shared** regime sits *outside* the product over tapes ⇒ the partial likelihood is
not separable in $\mathbf a=(a^{(1)},\dots,a^{(k)})$. State space is genuinely $(N{+}1)^k\times R$ —
**exponential in $k$**. At $k=166$, a dead end.

> **This is the precise statement of the intuition: a shared latent couples the tapes, and coupling
> destroys the factorization that makes SciPhy cheap.**

### The fix: augment, don't marginalize

$$\log p(\text{data},u)=\underbrace{\sum_{z=1}^{k}\log P(D_z\mid T,r,\beta,u)}_{k\text{ independent prunings — unchanged, }O(knN^2)}+\underbrace{\log p(u\mid \text{tree},\theta_\text{OU})}_{\text{tree-structured Gaussian, }O(Mn)}$$

Treat the latent as something you **optimize or sample**, not integrate out inside the likelihood.
Conditional on $u$, the tape product is restored. With a Brownian/OU prior, $u$ is a multivariate
Gaussian with covariance from shared branch lengths, and its log-density is **linear-time** via
Felsenstein's independent contrasts — the same pruning trick on a continuous trait. Fully
differentiable ⇒ $u$ joins node heights, $\lambda$, and spline knots in one gradient optimization.

**Structurally identical to a relaxed clock** — latent per-branch value with a shared-distribution
prior. SciPhy already cites the uncorrelated relaxed clock (ref 23) and random local clocks (ref
24). Established machinery, applied to $\xi$ instead of $\lambda$.

**Regularization / effective resolution.** $M\times2n$ latents vs. single-digit events per latent at
Park scale ⇒ badly under-determined alone. The OU prior saves it; its **length-scale is estimable
from the data**.

> **Honest characterization: you do not get reliable per-branch estimates. You get reliable
> estimates at whatever temporal/phylogenetic scale the data support — and the fitted OU
> length-scale TELLS you what that scale is.** A feature: the model reports its own resolution.

### ⚑ Build this version first: cell-type CREs jointly with the transcriptome

If the CRE is cell-type-specific, $u(\ell)$ is a function of cell state — **and cell state is
measured at every tip.** Both papers profile transcriptome and lineage from the same cells. The
latent is *observed at the leaves*, latent only internally ⇒ unsupervised becomes semi-supervised,
with a large identifiability gain.

Joint model: latent cell-state process diffuses on the tree; tapes record its CRE-driven output
along branches; scRNA-seq anchors it at the tips. **Tapes date the tree and record history; the
transcriptome pins down what the states are.**

Both papers do a weak post-hoc version (Park: "does pre-transplant state predict colonization,"
via barcode matching across separate analyses). A joint model does it properly and propagates
uncertainty. **More tractable than fully-latent signal inference, and a real paper.**

### Uncertainty

- **Laplace** at the MAP over $(u,\text{heights},r,\text{knots})$ — Hessian free from autodiff.
- **Variational** with structured Gaussian $q(u)$ — better calibrated under skew, still gradient-based.
- **Tape bootstrap** — non-parametric, expensive, honest.
- **Full MCMC on a small subset** as calibration. ⇒ **SciPhy's lasting role: not the production
  tool, but the reference implementation you validate against.**

⚠ Laplace on an $Mn$-dimensional latent field with a hierarchical prior is exactly where variance
is commonly under-estimated. **Check against MCMC on a small case before trusting it.**

## E.4 Perfect phylogeny, explained properly

**Core idea:** each edit is a character some cells have and others don't. Irreversible + never twice
independently ⇒ **the set of cells carrying an edit is exactly a clade.** Each edit *names* a clade;
collect all edits ⇒ collect all clades ⇒ that is the tree.

```
Cells:   A  B  C  D  E
edit-1:  ✓  ✓  ✓  ✓  ✓   → {A,B,C,D,E}         ┌── A
edit-2:  ✓  ✓  ✓  ·  ·   → {A,B,C}         ┌───┤
edit-3:  ✓  ✓  ·  ·  ·   → {A,B}         ┌─┤   └── B
edit-4:  ·  ·  ·  ✓  ✓   → {D,E}         │ └────── C
                                       ──┤
  every pair nested or disjoint          │   ┌──── D
  = LAMINAR family = a tree              └───┤
                                             └──── E
```

Read off directly, **no search**, linear in (cells × characters) with hashing.

**Incompatibility** = partial overlap:
```
edit-5:  ✓  ✓  ·  ·  ·   → {A,B}
edit-6:  ·  ✓  ✓  ·  ·   → {B,C}    ← overlap, neither contains the other
```
No tree has both. One is homoplastic or erroneous.

**Two Typewriter-specific facts:**

1. ⚑ **Within a tape, characters are automatically nested.** "carries $(\gamma_3,\gamma_7,\gamma_1)$"
   $\subset$ "carries $(\gamma_3,\gamma_7)$" $\subset$ "carries $(\gamma_3)$" — **prefix containment
   IS set containment.** Each tape's characters form a **trie**, already laminar. **A single tape can
   never be internally incompatible.** *Specific to sequential recorders; unordered recorders
   (GESTALT, Cas9 arrays) get no such guarantee.*
2. ⇒ **Conflict is purely cross-tape**, from homoplasy at rate $\approx\sum_i \xi_i^2$ — low single
   digits for Park's alphabet, $\ge12.5\%$ for the mouse's.

### Why it reduces the search space

**Reframe: you have $k$ trees, one per tape, each a coarse partial view. Merge them.** Any clade
surviving the merge is **fixed** — never searched. What remains is exactly the polytomies of that
skeleton: where tapes conflict, or where no tape distinguishes the cells.

This is a **constraint tree** — standard in ML phylogenetics (RAxML, IQ-TREE both support them),
and it cuts search cost dramatically. Resolve $C$ of $n-1$ internal nodes from compatible
characters ⇒ search only $n-1-C$. With $k=166$ and a large alphabet, $C$ should be most of them.

> **Yes: using all tapes at once shrinks the search space, by converting agreement across tapes from
> something you search into something you assert.** Pleasing symmetry with Eq. 9 — SciPhy intersects
> candidate sets *across children of a node*; this intersects clade families *across tapes*. Same
> operation, orthogonal axis.

### ⚑ And it gives a clean entry point for signal inference

With a dated skeleton, each edit event is assigned to **the branch where it first appears** — which
is exactly what the character's clade tells you. Free output:

| edit event | symbol | branch | approx. time |
|---|---|---|---|
| … | $\gamma_i$ | $\ell$ | uniform on $\ell$'s interval |

**Almost literally the raw material for estimating $f(t,\ell)$.** Plug-in estimator: bin by time
window and lineage group, tabulate symbol frequencies, normalize by the constitutive reference.

Biased and noisy (write times smeared along the branch; homoplasy misassigns; dropout thins). But:
**fast** (one pass, no optimization), an **excellent initialization** for the full latent model, a
**sanity check** before committing to heavy machinery, and an **immediate read on whether the signal
is visible at all**.

> **First experiment on real data:** run exactly this on Park et al.'s per-clone trees. Build the
> compatible skeleton, assign edits to branches, plot symbol composition vs. time. **Their pegRNAs
> aren't CRE-driven, so composition SHOULD be flat — a clean negative control for the whole
> estimator.** Structure there means homoplasy, dropout, or a bug. Know that before pointing this at
> ENGRAM data.

---

# Design notes III: bootstrap, and the fixed-topology / full-posterior hybrid

## F.1 Bootstrap — what it is

One dataset $D$, statistic $\hat\theta=s(D)$; want the sampling distribution of $\hat\theta$. Efron:
**use the data as a stand-in for $F$** and resample from it. Draw $n$ observations **with
replacement** → $D^*$; compute $\hat\theta^*$; repeat $B$ times; the spread estimates the sampling
distribution.

> **Everything hinges on the resampling unit.** It must be an independent, exchangeable draw. Get
> that wrong and the intervals are meaningless.

## F.2 ⚑ Here the unit is the TAPE, not the site

Felsenstein (1985; mouse ref. 55) resamples **alignment columns**, valid because sites are ~i.i.d.
given the tree.

**Sites are the WRONG unit for DNA Typewriter** — sites within a tape are sequentially coupled
(site 3 needs sites 1–2 edited). Independent resampling shreds the prefix structure and produces
probability-zero states. **Tapes are exactly the units the likelihood already treats as
conditionally independent.**

> **The factorization assumption and the bootstrap resampling unit are the same assumption.**

| Paper | Bootstrap | Unit | Measures |
|---|---|---|---|
| Typewriter | branch support | *"resampling with replacement of the intact TAPE-1 arrays"* | finite **recording capacity** |
| Park | parsimony diff $\delta=P_{\le k}-P_{\le k+1}$, 2,000 resamples | tapes | which topology the edits support |
| Park | distance-ratio CI, 2,000 resamples | **cells** | finite **cell sampling** |
| Mouse | — | — | *"computationally infeasible at this scale"* → edit-count support |

Keep the tape/cell distinction: *"would more tapes change my answer?"* vs *"would more cells?"*

### What it does and doesn't capture

**Captures:** variability from the finite number of tapes (or cells).
**Does NOT capture:** model misspecification, dropout uncertainty, clock-calibration error, tree
search failing to find the optimum. **A bootstrap can be beautifully tight around a wrong answer.**

- ⚠ **Resolution scales with $k$.** Park ($k=166$): fine-grained. Mouse ($k=11$): resampling eleven
  things ⇒ lumpy, coarse support values.
- ⚠ **Correlated tapes make it anti-conservative** — same root cause as the product likelihood's
  overconfidence (A5/A8). With a shared-latent model, consider bootstrapping at cell or clade level.
- ⚠ **Bootstrap support ≠ posterior probability.** Different objects; in phylogenetics they disagree
  systematically (bootstrap conservative, posterior clade probabilities anti-conservative).
  **Never report one as the other.**

### Variants for the toolkit

- **UFBoot** — resample per-unit log-likelihoods and reweight instead of refitting (RELL; the same
  Kishino–Hasegawa machinery Park already uses). Here: per-**tape** log-likelihoods. Orders of
  magnitude cheaper.
- **Jackknife** — leave-one-tape-out. With $k=166$, a useful influence diagnostic.
- **Parametric bootstrap** — simulate fresh data from the fitted model, refit, measure spread.
  ⚑ **This is literally your SBI simulator. Build one, get both.**

## F.3 The hybrid: point-estimate topology, full posterior on everything else

$$\underbrace{p(\theta\mid D)}_{\text{marginal}}=\sum_\tau p(\theta\mid D,\tau)p(\tau\mid D)\quad\longrightarrow\quad\underbrace{p(\theta\mid D,\hat\tau)}_{\text{conditional}}$$

Error depends on **how much $p(\theta\mid D,\tau)$ varies with $\tau$.**

**Why it's small for the parameters that matter.** Global parameters are functions of **coarse tree
summaries**: $\lambda$ from total edits over 45 days; $\xi$ from pooled symbol frequencies; $b,\delta$
from the LTT curve; $g_i(t)$ from symbol composition binned by time. If $\theta \perp \tau \mid
S(\tau)$ for a low-dimensional $S$ stable across the high-density region, conditioning ≈
marginalizing. **The topology uncertainty that matters is uncertainty in $S(\tau)$, not in $\tau$.**

> ⚑ **Self-consistency that helps:** the topology is most uncertain exactly where cells share no
> distinguishing edits — which is also where there is least information about $f$. **Topology
> uncertainty and signal information are anti-correlated.** Ambiguity concentrates where it costs
> least. *(Not airtight — a wrong DEEP split misassigns whole clades. Local tip unresolvability,
> which is nearly all of it, is benign.)*

| Parameter | Conditioning safe? | Why |
|---|---|---|
| $r,\xi,\beta$ | ✅ very | aggregate counts over the whole tree |
| $b,\delta$ | ✅ yes | function of the LTT curve |
| **Global signal $g_i(t)$** | ✅ yes | edits binned by time; binning is coarse |
| Root height, LTT, ages of **well-supported** clades | ✅ yes | topology-robust quantities |
| **Individual** node heights | ⚠️ conditional | "node 47" isn't the same node under another topology |
| **Fine-grained $u_i(\ell)$** | ❌ no | change the topology and "lineage $\ell$" changes |
| Clade regimes on **well-supported** clades | ✅ mostly | the clade is stable |

⇒ Same split as before: **global trajectories safe; per-lineage signal at fine resolution not.**

### ⚑ Precedent: this IS phylogenetic comparative methods

PGLS, ancestral state reconstruction, OU/Brownian trait models — **essentially all of PCM conditions
on a fixed tree** and infers trait-evolution parameters with full uncertainty. **Your latent
$u_i(\ell)$ with an OU prior *is* a phylogenetic comparative method.** And that field's standard fix
for the caveat is exactly option 2 below: repeat over a sample of plausible trees and pool.

### Menu, increasing rigor

1. **Plug-in.** Fix $\hat\tau$, infer the rest. Cheapest; understates uncertainty; honest if stated.
2. **Multiple imputation over topologies.** $B$ plausible topologies (tape bootstrap / top-$B$ from
   search / CCD-SBN samples), conditional posterior on each, pool by Rubin's rules:
   $$\text{Var}_\text{total}=\overline{\text{Var}_\text{within}}+\left(1+\tfrac1B\right)\text{Var}_\text{between}$$
   ⚑ **The two variance components separate** — you can report exactly how much uncertainty comes
   from not knowing the tree.
3. **Weighted likelihood bootstrap** (Newton & Raftery 1994; cf. Rubin's Bayesian bootstrap 1981).
   Draw Dirichlet$(1,\dots,1)$ weights over the $k$ tapes, maximize the **weighted** likelihood
   (topology included), repeat. Approximates posterior sampling. Each replicate = one MAP fit ⇒
   $B=100$ gives an approximate **joint** posterior *including topology* for 100 optimizations
   instead of $10^9$ MCMC steps. **Most elegant option; look hard at it.**
4. **Full MCMC** on a small subset as the calibration reference for 1–3.

### ⚑⚑ Why this is the right architecture

**Conditional on topology, the remaining posterior is continuous, smooth, probably unimodal** — node
heights (ratio transform), $\lambda$, $\xi$, the OU field $u$, $b$, $\delta$, all with autodiff
gradients. ⇒ **HMC/NUTS, not random-walk Metropolis.** Efficiency gap is enormous: random walk needs
~$O(d)$ steps per effective sample, HMC ~$O(d^{1/4})$.

> **SciPhy's $10^9$ iterations are not about the continuous parameters — they're the
> trans-dimensional random walk over discrete topology space.**
>
> **Fixing the topology converts an intractable discrete–continuous problem into a standard smooth
> continuous one, where HMC gives honest, well-calibrated posteriors in ~$10^4$ gradient
> evaluations.**

Near-ideal for the signalling goal: the thing you most want a posterior over (the latent field $u$,
the trajectory $g(t)$) is the continuous part HMC handles well; the thing you give up (topology)
contributes least to it.

### Reporting structure

- **Topology:** point estimate + per-branch support (UFBoot / tape bootstrap).
- **Node dates:** credible intervals conditional on topology, or pooled via Rubin's rules.
- **Signal, $\lambda$, $\xi$, $b$, $\delta$:** full posteriors.

Strictly more informative than either paper — Park's dates come from a bootstrap conditional on a
**fixed externally calibrated clock**, a considerably stronger conditioning than this.

### Risks

- **Multimodality in node heights** — saturated tapes flatten the likelihood in recent times. HMC
  handles curvature, not multimodality. Check with multiple initializations.
- **Laplace/VI understating variance** on the high-dimensional $u$ field.
- **Deep topology errors** — benign case is unresolved tips; dangerous case is a misplaced deep
  split. Watch bootstrap support on deep branches; it's cheap.
- ⚑ **The validation you owe:** on ≥1 small dataset, run full MCMC and show the conditional
  posterior matches the marginal. **This is the experiment that licenses the whole architecture —
  and it's a simulation study, doable before any ENGRAM data exists.**

---

# Method primers: HMC, OU fields, PCM

> **They layer:** PCM is the *framework*, OU is the *model* inside it, HMC is the *inference engine*.

## G.1 Hamiltonian Monte Carlo

**Problem it solves.** RWM proposes $\theta'=\theta+\epsilon$; step size must shrink like $d^{-1/2}$
for decent acceptance, and being a *random walk*, covering distance 1 takes $(1/\epsilon)^2\sim d$
steps ⇒ $O(d)$ likelihood evals per effective sample. The chain diffuses; it doesn't travel.

**Construction.** Potential energy $U(\theta)=-\log p(\theta\mid D)$. Add auxiliary **momentum**
$\rho$ with $K(\rho)=\tfrac12\rho^\top M^{-1}\rho$; Hamiltonian $H=U+K$. Then

$$p(\theta,\rho)\propto e^{-H(\theta,\rho)}=e^{-U(\theta)}\cdot e^{-K(\rho)}$$

**factorizes** ⇒ $\rho$ is independent Gaussian noise; marginalizing recovers the target exactly.
Momentum costs nothing statistically — pure algorithmic scaffolding.

**Algorithm.** (1) draw $\rho\sim N(0,M)$; (2) $L$ **leapfrog** steps of size $\epsilon$:
$$\rho\leftarrow\rho-\tfrac\epsilon2\nabla U(\theta);\quad \theta\leftarrow\theta+\epsilon M^{-1}\rho;\quad \rho\leftarrow\rho-\tfrac\epsilon2\nabla U(\theta)$$
(3) accept endpoint w.p. $\min(1,e^{H(\theta,\rho)-H(\theta',\rho')})$.

**Why it works.** Hamiltonian dynamics conserve energy ⇒ perfect integrator would give acceptance 1.
Leapfrog is **symplectic**: conserves phase-space volume *exactly* (needed for detailed balance) and
energy *approximately*, with error that oscillates rather than drifting ⇒ high acceptance even for
long trajectories. The particle rolls **along** posterior contours instead of random-walking across.

| | Step size | Evals per effective sample |
|---|---|---|
| RWM | $\epsilon\sim d^{-1/2}$ | $O(d)$ |
| HMC (tuned) | $\epsilon\sim d^{-1/4}$, $L\sim d^{1/4}$ | $O(d^{1/4})$ |

Optimal acceptance: 0.234 (RWM) vs ~0.65 (HMC). **NUTS** removes tuning — extends each trajectory
until it doubles back, adapts $\epsilon$ by dual averaging in warmup. Stan / NumPyro.

**Our dimension** after fixing topology: ~999 heights + 166 rates + knots + ~16,000 OU field
$\approx d\sim1.7\times10^4$. $d^{1/4}\approx11$ vs $d\approx17{,}000$ ⇒ **~1,500×**. This is why
$10^4$ gradient evals is realistic where $10^9$ random-walk steps was not.

**Limitations:**
- ⚠ **Requires continuous, differentiable parameters. Cannot touch discrete ones.** *This is why
  topology must leave the sampler — a hard constraint, not a preference.*
- Explores one mode superbly; jumps between separated modes poorly.
- ⚠⚠ **Funnels — you WILL hit this.** Hierarchical models produce Neal's funnel: small variance
  hyperparameter ⇒ the latent field must concentrate into a narrow neck; no single step size works
  for both neck and mouth. **The OU field with estimated variance is exactly this geometry.**
  **Fix (build it in from the start): non-centered parameterization.** Sample standardized
  $z\sim N(0,I)$ and set $u=\sigma Lz$ with $LL^\top=\Sigma$. For a tree-structured OU the natural
  whitening is via **standardized independent contrasts** — the same Felsenstein recursion you
  already have.

## G.2 Ornstein–Uhlenbeck processes and OU fields

$$dX_t=\underbrace{\theta(\mu-X_t)dt}_{\text{pull toward }\mu}+\underbrace{\sigma\,dW_t}_{\text{wobble}}$$

$$X_t\mid X_0\sim N\!\left(\mu+(X_0-\mu)e^{-\theta t},\ \tfrac{\sigma^2}{2\theta}(1-e^{-2\theta t})\right)$$

Stationary: $N(\mu,\sigma^2/2\theta)$. Autocorrelation $e^{-\theta|t-s|}$.

- **Mean-reverting** — unlike BM, stays bounded. Right for log CRE activity.
- **One interpretable knob:** $1/\theta$ = **characteristic memory time**. $\theta\to0$ ⇒ Brownian
  motion (perfect heritability); $\theta\to\infty$ ⇒ white noise (branches independent). **Smoothly
  interpolates "fully heritable" ↔ "not heritable," and is estimated from data.**
- **Gaussian throughout** ⇒ tractable and differentiable.

**On a tree:** run down every branch from the root; lineages share trajectory until divergence. With
$t_{ij}$ = time since MRCA (ultrametric):

$$\text{Cov}(X_i,X_j)=\frac{\sigma^2}{2\theta}e^{-2\theta\,t_{ij}}$$

**"Field"** = a value at every branch/node rather than one scalar — a random field indexed by the
tree.

**Computation.** Naive: $2n\times2n$ covariance, $O(n^3)$ — dead. But OU on a tree is
**Gauss–Markov** ⇒ log-density in **$O(n)$** by post-order recursion (Felsenstein's independent
contrasts for BM; OU generalization).

> ⚑ **Pleasing symmetry: two pruning recursions on the same tree** — discrete CTMC for tape data
> (Eqs. 8–12), Gaussian for the latent field. Both linear time, both differentiable, both feeding
> the same HMC gradient.

### ⚑⚑ The payoff

$$\boxed{\ 1/\theta\ =\ \text{the heritability timescale of the signalling state}\ }$$

**Not a nuisance parameter — the scientific answer to Park et al.'s central question.** Their paper
argues metastatic potential is *"a pre-existing, heritable cell state,"* established via clonal
barcodes and qualitative comparison. **An OU field turns that into a number with a credible
interval:** this state persists X days of lineage time. Long relative to the experiment ⇒ heritable;
short ⇒ transient and niche-acquired.

**Extension:** **multi-regime OU** (OUCH/SURFACE family) lets $\mu$ shift at certain nodes — the
continuous bridge to the discrete-regime idea, and probably a better fit to switch-like signalling
than plain OU.

## G.3 Phylogenetic comparative methods (PCM)

**Founding problem (Felsenstein 1985).** Regress one trait on another across related taxa and the
residuals are correlated — close relatives resemble each other through **shared ancestry**, not
causation. $n$ observations are treated as $n$ independent points when the effective sample size is
far smaller. Significance wildly overstated. **Phylogenetic non-independence.**

| Method | What it does |
|---|---|
| **Independent contrasts** | transform $n$ correlated tips into $n-1$ contrasts independent under BM; $O(n)$ |
| **PGLS** | regression with tree-derived covariance; generalizes contrasts |
| **Ancestral state reconstruction** | infer trait values at internal nodes |
| **Trait-evolution models** (BM, OU, early-burst) | compare to ask *how* a trait evolves |
| **Phylogenetic signal** (Pagel's $\lambda$, Blomberg's $K$) | how much trait variance does the tree explain? |

**Four connections:**

1. **Your latent $u_i(\ell)$ with an OU prior literally IS a PCM** — ancestral state reconstruction
   of a continuous trait under OU. Mature, solved, decades of theory and software. **You are not
   inventing this part.**
2. **PCMs condition on a fixed tree — essentially universally.** The precedent for the hybrid
   architecture; and their standard remedy (repeat over a tree sample, pool) is exactly the
   multiple-imputation option (§F.3).
3. ⚑ **All three papers already do informal PCM without the framework.**
   - Mouse: *"A null distribution was generated by permuting cell-type labels across tips within
     each blastomere subtree"* — a phylogenetic permutation test, the nonparametric cousin of PCM.
   - Park: *"349 genes recurrently heritable in vitro"*; *"phylogenetic fitness ... correlated with
     a hypoxia-, glycolysis- and partial-EMT-centered axis."* **Heritability along lineage IS
     phylogenetic signal**, computed with bespoke clone-level statistics.

   ⇒ **Adopting the framework properly — model-based rather than permutation-based, parameters with
   units and credible intervals — is a clean contribution on its own.**
4. **What's genuinely novel in your version.** Standard PCM assumes the trait is **observed at the
   tips**. Here CRE activity is **never directly observed** — inferred only from its effect on a
   time-integrated recorder. **PCM with a latent trait seen through a noisy, temporally-smeared
   readout.** A real extension (related to threshold/liability and phylogenetic latent-variable
   models, but not identical). *For cell-type CREs the trait is partly observed via scRNA-seq at the
   tips — better identified, and the reason to start there.*

```
  PCM          framework:  traits evolving on a tree, non-independence handled
    │
    ├── OU     model:      mean-reverting, heritable, θ = memory timescale
    │                      log-density in O(n) by pruning
    │
    └── HMC    engine:     gradients → O(d^{1/4}) instead of O(d)
                           requires everything continuous ⇒ topology fixed
```

---

# Pilarski, Stadler & Seidel 2026 (PLOS Comp Biol 22:e1014370)

*"Assessing the inference of single-cell phylogenies and population dynamics from CRISPR lineage
recordings."* A systematic simulation benchmark of the exact framework we plan to extend, **by the
group that built it.**

**Design.** 100 time-scaled phylogenies (20–200 tips, $T=40$) under five tree-generating processes,
from **fully synchronous divisions** (early embryogenesis) to **stochastic birth–death**. Barcodes
simulated under both recorders — **TiDeTree** (non-sequential Cas9 scarring) and **SciPhy**
(sequential tapes). Joint inference in BEAST 2; compared to truth. $10^8$ steps or ESS>200, 10%
burn-in. Tree priors: BDSKY (homogeneous), BDMM-Prime (multi-type).

## ⚠⚠ Two setup details that matter enormously

> **They FIXED the insert probabilities $\xi$ to their true values.** *"We inferred the editing rate
> and fixed the scarring and insert probabilities to true values, because we reasoned that the
> occurrence and relative frequency of each editing outcome can be quantified in real experiments
> using sequencing data."*

**Every performance number assumes $\xi$ is known.** Our entire ENGRAM extension is about inferring
a *time- and lineage-varying* $\xi$. This study says **nothing** about whether that is feasible —
simultaneously the opportunity and a caution against reading its optimism across.

Also: **$\rho$ fixed to truth** *"due to identifiability reasons"* ⇒ sampling proportion and
birth/death rates are confounded. Origin fixed at 40 (consistent with §2a.3).

## Headline results

**1. Barcode diversity is the master variable.** Kendall's $\tau=-0.68$ between proportion of unique
barcodes and wRF distance. wRF **halved** from 5→20 targets, improved further at 40 (matches the
$\sqrt n$ expectation).

**2. ⚠ But it saturates hard.** *"Once nearly all cells acquired unique barcodes, the lineage
relationships could be resolved almost perfectly. Further accumulation of edits resulted in only
minimal improvements."* Sequential at $r=0.05$ → 40% of sites filled, wRF 0.10; at $r=0.15$ → 90%
filled (2× informative sites), wRF 0.09. **Steep diminishing returns — there is a capacity target,
and exceeding it buys almost nothing.**

**3. ⚑⚑ Sequential editing improves TOPOLOGY but NOT TIMING.** Matched comparison: 20 independent
sites vs. **one 20-site tape**, editing rate tuned to 0.45 so barcode diversity and edits/cell match.

| Metric | Sequential vs. non-sequential |
|---|---|
| wRF distance | **significantly better** (Wilcoxon, $p<0.001$) |
| Topological RF, Shared PI | **significantly better** |
| Branch-length distributions (Wasserstein) | slightly better ($p<0.01$); **KS not significant** |
| Division / death / growth rate **bias** | **not significant** |
| HPD widths | **not significant** (except death rate, $p<0.05$) |

*"Sequential editing primarily improved the reconstruction of lineage relationships, while the
temporal aspect of the trees remained largely unaffected. Consequently, sequential editing did not
consistently increase the accuracy nor reduce the uncertainty in the estimated cell division and
death rates."*

> **Order information buys topology, not dates.** The specific advantage of DNA Typewriter over Cas9
> recorders is narrower than the field's rhetoric. **If the scientific question is a rate or a
> timing, the sequential architecture is not where the gain comes from.**

**4. Division rates robustly inferrable; death rates essentially not.** Division: high coverage, low
bias, HPD shrinking >90% vs prior across all scenarios. Death: **coverage 0%** in several scenarios,
relative bias up to +1. Principled reason — death signs only via the LTT curve, very noisy in small
populations, worse under partial sampling.

**5. Misspecification is real and hits timing.** Memoryless birth–death fitted to synchronous
divisions ⇒ overdispersed branch lengths, death rate under-estimated (with death) / over-estimated
(without). But: *"the inferred branch length distributions matched the truth better for sequential
recordings, indicating that more informative data yielded estimates of cell division timings that
were more robust to model misspecification."* ⇒ **More information buys robustness to a wrong prior.**

**6. The $\ln 2$ correction.** Synchronous growth is exponential base 2; birth–death assumes base $e$
⇒ per-lineage and population-level rates differ by exactly $\ln2$. *"Birth rates inferred under the
birth-death model can be divided by $\ln 2$ to estimate the doubling time of cell populations growing
by synchronous and regular cell divisions."* Directly actionable for early embryonic dating.

## ✅ Three results validating our plan

**Filtering to a complete submatrix works — for parameters.** Simulated silencing (0.01) + dropout
(0.2) → 40–60% missing; filtered exactly as SciPhy does (→ subsets of 35–250 cells, ~5 of 20
targets/tapes retained). **Division rate, editing rate, tree height and length recovered in ≥90% of
simulations while retaining only ~2.5% of the data** — *"albeit at the cost of reconstructing lineage
trees for only a small subset of cells."*

> **Direct simulation evidence for the hybrid architecture: population-level and editing parameters
> survive drastic subsetting; only the fine-grained tree is lost.** The "which questions need which
> cells" table (§ comparison) now has empirical backing.

**Latent discrete states on trees already work.** Multi-type analysis (4 types, BDMM-Prime, terminal
vs chain-like transitions): type-specific division/death and transition rates recovered with >80%
coverage; **ancestral cell types correct at a median 97.8% (terminal) / 87.4% (chain-like) of
internal nodes.** Workflow: infer **tip-typed** trees, then **stochastic mapping** for ancestral
states — exactly the two-pass structure of §1c.5½.

> ⚑ **But read the difference — it sharpens our novelty claim.** In BDMM-Prime the type affects the
> **tree prior** (per-type birth/death/transition rates), **not the editing model**. Types never
> enter the tape likelihood ⇒ tapes stay independent ⇒ **the coupling problem of §E.3 never arises.**
> Our ENGRAM model puts the latent state *inside* the editing likelihood, shared across tapes.
> **The field has solved the uncoupled version; the coupled version is genuinely ours.**

**They name our problem as the open one.** *"It remains an open question how many cells are needed to
get accurate estimates of cell division, differentiation, and death rates, and thus, to what number
of cells the methods need to scale."*

## ⚠ Two results that complicate it

**Runtime is worse than we thought.** *"Each MCMC chain required at least a few hours... Inference
with the multi-type phylodynamic model lasted up to three weeks."* On trees of **20–200 tips** — not
1,000, a couple of hundred. **Adding a latent state on the tree explodes the cost.** Direct warning
about the latent-field plan under MCMC; strongest argument yet for the gradient route.

**More data does not fix misspecification.** *"If the model fails to capture the true dynamics, more
data may lead to over-confident but biased posteriors."* Their own biases from synchronous-division
misspecification **persisted and became more evident at increased sample size.**

⚠ **Cautionary subtlety:** filtering *improved* branch-length recovery on synchronous trees, because
*"filtering effectively induced non-uniform sampling of lineages, making branches in the synchronous
trees more irregular and more compatible with the birth-death model used for inference."*
**Two errors partially cancelling** — exactly what makes a pipeline look validated when it isn't.

## ⚠ Corrections to earlier claims in these notes

**1. "Nobody has done ML tree search on lineage recordings" was too strong.** Their ref. 46 is
**Chu, Mai, Schmidt & Raphael, "Maximum likelihood inference of time-scaled cell lineage trees with
mixed-type missing data using LAML," Genome Biology 2025** — cited precisely for the point about
explicit error models: *"Explicitly incorporating error processes into the mechanistic, inferential
models of CRISPR editing (as demonstrated in [46]) would circumvent the need for extensive
filtering."*

> **Accurate version: LAML handles Cas9-style *unordered* recorders with mixed-type missing data.
> Nobody has done it for *sequential* recorders under the SciPhy likelihood.**
> **⇒ Read LAML before building anything — it may be the right thing to extend rather than
> reimplement.**

**2. The authors' own scalability bet is not ours.** Their *"recent computational advancements"* are
refs 53–54: **Bouckaert et al., "Improving the scalability of Bayesian phylodynamic inference through
efficient MCMC proposals"** (bioRxiv 2025) and **Delphy, "scalable, near-real-time Bayesian
phylogenetics for outbreaks"** (bioRxiv 2025). **Better proposals + faster engineering, staying
inside MCMC.** Not VI, not ML, not SBI. An unoccupied niche for us — *and* a reason to understand why
they didn't go there before assuming they overlooked it.

## Reading list this surfaces

- ⚑ **Schiffman et al., "Defining heritability, plasticity, and transition dynamics of cellular
  phenotypes in somatic evolution," Nat Genet 2024** (ref. 33). Title suggests the closest existing
  work to the **OU heritability-timescale idea** (§G.2). **Joshua Schiffman is a coauthor on Park et
  al.** ⇒ the metastasis paper's heritability analysis very likely uses this framework.
  **Read before proposing our own version.**
- **Zwaans, Seidel, Manceau & Stadler, "A Bayesian phylodynamic inference framework for single-cell
  CRISPR/Cas9 lineage tracing barcode data with *dependent target sites*," Phil Trans R Soc B 2025**
  (ref. 24). **Dependent target sites = the coupling problem (§D.3, §E.3).**
- **Mulberry & Stadler, "Strategies for resolving cellular phylogenies from sequential lineage
  tracing data," Theor Popul Biol 2026** (ref. 22). Theory specific to sequential recorders.
- **Wang, Zhang, Khodaverdian & Yosef, "Theoretical guarantees for phylogeny inference from
  single-cell lineage tracing," PNAS 2023** (ref. 21). Identifiability theory.
- **Chadly et al., "Regenerative base editing enables deep lineage recording," bioRxiv 2026**
  (ref. 47) — the "hypercascade" recorder; sequential-adjacent, higher capacity.

---

# Literature deep-dive: Mulberry, Zwaans, Schiffman

## H.1 ⚑⚑ Mulberry & Stadler 2026 (Theor Popul Biol 168:32–43) — the foundation

*"Strategies for resolving cellular phylogenies from sequential lineage tracing data."*
**Analytic bounds on P(exact tree reconstruction) for SEQUENTIAL recorders.** Built on the same
generative model as SciPhy (censored Poisson, citing Seidel et al. explicitly).

### Their parameterization — adopt it

| | Meaning |
|---|---|
| $m$ | tape copy number per cell |
| $k$ | target sites per tape |
| $j$ | number of unique characters |
| $\lambda$ | editing rate |
| $\xi_i$ | insertion probabilities |
| $q=\sum_i\xi_i^2$ | **collision probability** (two independent edits → same character) |
| $\ell$ | lower bound on minimal internal branch length = **the resolution you demand** |

$q$ is exactly the $\sum \xi_i^2$ of §1c.2. Two bounds: $B_\infty$ (assumes $q=0$),
$B_q$ (general). Regimes $q=1/64$ (high diversity) vs $q=1/4$ (low). $B_\infty$ *"significantly
overestimates the accuracy under low rates when $q=1/4$"* ⇒ **homoplasy at the mouse's alphabet
size is not a rounding error.**

**Their distance is again summed per-tape shared prefix:** $D_{a,b}=\sum_i(k-k_i')$.
⇒ **Third independent appearance** (SciPhy's lcp, Park's "prefix divergence", this).
**Clearly the canonical statistic for sequential recorders.**

Reconstruction criterion: triplet condition $D_{a,b}<\min(D_{b,c},D_{a,c})$ for all triplets
(Aho et al. 1981; Gascuel & Steel 2016). Distance methods (NJ, UPGMA) **outperform** the triplet
score ⇒ bounds are conservative.

### The four results that matter most

**⚑ (a) Multiple editing rates dramatically improve reconstruction.** §3.4 / Fig. 7: with fast
early and slow late divisions, **two** editing rates collapse RF distance from ~1,900 to ~0 at
$\alpha=0.1$. Qualitative, not marginal.

> **Reframes row A1: rate variation is not just a nuisance to model — it is a design feature that
> buys resolution.** If you can *engineer* rate diversity, do. And where it happens naturally
> (mouse ZGA burst → quiescence → resumption), a model capturing it gains a lot.

> ⚠ **CORRECTION (see §H.6.13, after the close read).** The mechanism in §3.4 is **not** a rate
> that varies *in time* on one tape. It is **two populations of tapes with different constant
> rates in the same cell** ($\lambda_1$ on $m_1$ tapes, $\lambda_2$ on $m_2$, $m_1+m_2\le m$).
> Keep these separate: **A1** = temporal, uncontrolled, shared across tapes (must be *modelled*);
> **§3.4** = across-tape, engineered, constant per tape (a *design lever*). They are different
> objects that happen to buy the same thing.

**⚑ (b) Accuracy depends on $\xi$ ONLY through $q$.** Fig. 8: skewed over 21 characters vs uniform
over 12, matched at $q=1/12$ ⇒ statistically indistinguishable accuracy. **The shape of $\xi$
carries no topological information beyond its second moment.**

⇒ **Good news: composition information and topology information are largely orthogonal.** Varying
$\xi(t)$ for signal recording does not degrade tree reconstruction *provided $q$ stays comparable*.

**⚑⚑ (c) But this implies a design tension nobody has stated.** Strong ENGRAM signal ⇒ one symbol
dominates ⇒ $f$ peaked ⇒ **$q$ rises** ⇒ homoplasy rises ⇒ tree degrades.

> **A STRONG SIGNAL DEGRADES THE PHYLOGENY.** The more informative the ENGRAM readout at a given
> moment, the less resolvable the lineage at that moment.
>
> **Mitigation:** include enough always-on background symbols / a high-diversity constitutive
> pegRNA pool to keep $q$ low even when signal symbols spike.

**Their $q$-estimator makes this checkable:** $\hat q=\sum c_{ij}/\sum n_{ij}$ (independent
co-occurrences / sites where both sequences have an edit). Validated Fig. 9. **Use this.**

**(d) Sensitivity is to resolution $\ell$, NOT sample size $n$.** For $k=5$: ~30 tapes for a
512-tip tree at $\ell\approx0.1$; ~24 for 100 tips (barely different). Halve resolution to
$\ell=0.05$ ⇒ **50** tapes for the same 100 cells. *"These results are much more sensitive to
$\ell$ rather than $n$."*

> ⇒ **Reframes the scale question: the binding constraint is not how many cells you analyze, it's
> how finely you want to resolve time.** Exactly the currency the signal-history project trades in.

**Bottom line on current tech:** Typewriter as published ($k=5$, $j=64$, $m\approx20$) is
*"generally not yet in a regime where exact reconstruction is possible for trees on the order of
1000 cells."*

**Their own stated future work:** *"a probabilistic approach to a relaxed problem"*; *"an
alternative approach [that] could take into account both the underlying cellular dynamics and the
effect of sampling on the branch length distribution."* **An open invitation.**

## H.2 Zwaans et al. 2025 (Phil Trans R Soc B 380:20230318) — GABI

GESTALT → Bayesian in BEAST 2, extending **GAPML** (Feng et al. 2021, penalized ML) by adding a
molecular clock rate $\lambda$ for **absolute** time scaling.

### ⚠ Correction 1: "dependent target sites" ≠ our coupling problem

It means GESTALT's **inter-target deletions** — one Cas9 double-cut spans and overwrites adjacent
sites, so an edit at one site destroys information at others. A *within-barcode structural*
dependence, **not** the shared-latent coupling across tapes (§E.3). Instructive as an example of
handling non-independence; different problem.

### ⚑⚑ Correction 2: the fixed-topology hybrid is ALREADY PUBLISHED, by this group

**Pipeline 3** takes a maximum-likelihood tree from GAPML, **fixes the topology**, and runs GABI for
branch lengths + all parameters. Their justification is nearly verbatim §F.3:

> *"This is consistent with many pipelines for phylodynamic inference in both fields of
> macroevolution and epidemiology, where pre-existing, or separately inferred, fixed tree
> topologies are used to infer population dynamic parameters. This is particularly relevant in the
> context of Bayesian inference, where the exploration of the exponentially growing lineage tree
> space quickly becomes intractable."*

And it works: fixed-topology division rate **2.6 h⁻¹** vs **2.2 h⁻¹** from full joint inference —
consistent, and handling 100 sequences where full inference was limited to 20.

⇒ **The architecture is de-risked and accepted practice in this exact literature.** Good (no need to
justify from first principles); cautionary (**not the novel part**).

### ⚑ A fourth advantage of sequential recorders I'd missed

GABI's stated bottleneck: *"the ancestral state pruning algorithm, which scales poorly in the
presence of large inter-target deletions."* **SciPhy's ancestral state sets are bounded by $N+1$
because editing is sequential and irreversible** — the prefix structure caps them.

> **Sequential recorders are computationally better behaved, not just more informative.** Fourth
> advantage over GESTALT, alongside no DSBs, order information, and constant hazard.

**Also:** pooling replicates by multiplying likelihoods gave **26% narrower HPDs**. Park's 75 clones
share editing parameters ⇒ same trick applies, 75 independent replicates of one editing process.

**Their other limitations:** hundreds of sequences max, days-to-weeks runtime; call explicitly for
relaxed/time-dependent clocks; note birth–death and coalescent both fail on synchronous divisions.

## H.3 Schiffman et al. 2024 (Nat Genet 56:2174–2184) — PATH

**PATH** = Phylogenetic Analysis of Trait Heritability. Adapts **Moran's $I$** into *phylogenetic
autocorrelation* — how much close relatives resemble each other phenotypically — as a quantitative
**heritability vs plasticity** measure. Plus *cross-correlation* between states. Critically:

> *"we revealed a direct link between phylogenetic correlations and phenotypic state transitions,
> and thus can transform phylogenetic correlations to state transition inferences with high
> precision."*

**PATHpro** adds state-specific proliferation rates (multitype branching process). Benchmarked vs
SSE-MLE: comparable accuracy, **dramatically faster** — PATHpro <40 s vs SSE-MLE >5 h; scales far
better in both #states and #cells (Fig. 2e,f). Can even **impute branch lengths** from the sampling
rate when barcodes are too short.

### ⚠⚠ Honest reassessment of the OU heritability-timescale idea (§G.2)

Transition rates have units of inverse time ⇒ $1/\text{rate}$ **is** a heritability timescale.
**PATH got there first**, with discrete states + Markov transition rates rather than a continuous
OU. And **Joshua Schiffman is a coauthor on Park et al.**, so their heritability analysis very
likely *is* PATH. Their EMT result is directly parallel: terminal states heritable, intermediate
hybrids plastic, with Simeonov's peak-metastatic EMT scores 20–22 in the **most plastic** bin.

> **Do NOT propose "quantify heritability timescale" as the contribution.** It exists, it's fast,
> it's in Nature Genetics.

### Where our version IS genuinely different — a clean line

> **PATH requires the phenotype to be OBSERVED at the tips.** Autocorrelation is computed on
> measured leaf values.
>
> **For ENGRAM, CRE activity is NEVER observed.** It is recorded in the tape, integrated over time,
> read out indirectly. **There is no leaf value to correlate.** PATH cannot be applied to ENGRAM
> data at all.

**Corollary to act on: for the cell-type-driven case, USE PATH, don't rebuild it.** "Is this
transcriptional state heritable" → 40 seconds. Reserve our machinery for what PATH structurally
cannot reach.

**Their limitations we'd inherit:** Markovian transitions; multitype branching process;
**near-equilibrium state frequencies**; and *"incomplete sampling leads to a decline in cellular
relatedness... underestimating the true heritability"* — bites hard at $\rho\sim10^{-3}$.

## H.4 ⚠ Consolidated ML-landscape correction

| Method | What it does |
|---|---|
| **GAPML** (Feng et al. 2021, AoAS) | penalized ML for GESTALT trees; relative time only |
| **ConvexML** (Prillo et al. 2023) | *"scalable and accurate inference of single-cell chronograms"* — convex optimization for dated trees |
| **LAML** (Chu et al. 2025, Genome Biol) | ML **time-scaled** trees with **mixed-type missing data** |
| **GABI Pipeline 3** (Zwaans et al. 2025) | fixed ML topology + Bayesian branch lengths/parameters |
| **PATH / PATHpro** (Schiffman et al. 2024) | fast correlation-based alternative to MLE for state dynamics |

> **The gap, stated accurately: ML and fixed-topology-Bayesian pipelines for lineage tracing exist
> and are mature. NONE handles a sequential/ordered recorder under the SciPhy likelihood, and NONE
> handles a latent, time-varying, tape-recorded covariate. Those two together are the
> contribution.** Narrower than "make Bayesian lineage tracing scale," and far more defensible.

## H.5 Reading priority

1. **Mulberry & Stadler** — foundational for our exact model. Adopt their notation; take the
   $q$-estimator; take the multi-rate result seriously (both an argument for the extension and a
   design recommendation). **§3.4 and Fig. 7 are the highest-value pages.**
2. **Schiffman et al. (PATH)** — establish what's already claimed about heritability, to position
   against it and to know when to just use PATH.
3. **Zwaans et al.** — skim the model, read §3 pipelines carefully; **Pipeline 3 is our
   architecture, already validated.**

**Then, before writing code:** **LAML** (closest existing ML implementation, with the error model we
need) and **ConvexML** (convex-optimization framing of dating — possibly directly reusable for the
node-height step).

---

# Experimental applications: what data to generate

Papers: **Marjanovic et al. 2020** (Cancer Cell 38:229) HPCS discovery; **Chan et al. 2026**
(Nature 651:231) HPCS lineage tracing + ablation; **Hamazaki et al. 2024** (Nat Cell Biol 26:1790)
RA-gastruloids (Shendure lab).

## I.1 What each tier actually tests

| Tier | Signal structure | Computational demand |
|---|---|---|
| **1. Culture, defined program** | one global $\xi(t)$, shared by all lineages | **Easiest.** Pruning suffices; no latent field; tree is a nuisance |
| **2. Gastruloid** | lineage-specific exposure, **known answer** | Latent field, but shallow tree (~6 divisions) → few latents |
| **3. Cancer** | lineage-specific, unknown, deep tree, dropout | **Everything in the A-table at once** |

### ⚠ Tier 1 as described tests only the TEMPORAL axis

A well-mixed culture ⇒ **every lineage sees the same signal**. You recover $\xi(t)$ — the global case
needing none of the lineage machinery. Necessary positive control; won't test whether you can
identify *which lineages responded*.

**Cheap fix: make the response heterogeneous.** Mix two clones of differing reporter sensitivity
(or responsive + pathway-null), pool, apply the program, ask whether the model partitions the tree
by response. **Ground truth on both axes**, repeatable ×50 — the dataset you want for
simulation-based calibration anyway.

## I.2 Tier 2 — why RA-gastruloids fit unusually well

1. **The RA protocol is discontinuous and the discontinuity is functionally essential.** RA 0–24 h,
   withdrawn, again 48–120 h. *"the first pulse of RA (0–24 h), together with Matrigel starting at
   48 h, was sufficient to induce structures resembling a neural tube flanked by somites, the
   second pulse of RA (48–120 h) was not."* ⇒ **The two pulses are not interchangeable** — a
   temporal-order question with a known, non-trivial answer.
2. **Controlled AND endogenous signals in one system.** RA/CHIR applied on a known schedule (ground
   truth) while the gastruloid generates its own gradients — *"operational signalling gradients
   along both the dorsoventral (neural crest) and mediolateral (intermediate mesoderm) axes."*
3. ⚑ **The question to build it around:**

> **Did the neural-tube vs. somite fate split occur BEFORE or AFTER the first RA pulse?**

NMPs are bipotent; RA at 0–24 h maintains bipotency; the neural/mesodermal decision follows.
ENGRAM+Typewriter gives both halves in one molecule: **divergence time** (tape position) and **RA
exposure** (symbol identity). Currently answerable only from population-level marker dynamics.

⚠ 96–120 h ≈ 5–7 divisions ⇒ shallow tree, few internal branches. Fine for divergence-ordering;
poor for deep phylogenetic structure. *SciPhy's own gastruloid analysis is precedent here.*

## I.3 Tier 3 — the cancer application

### ⚑ The gap in Chan et al. is unusually clean

Their tracing is **Cre-based**: a tamoxifen pulse marks Slc4a11⁺ cells *at that instant*; follow
descendants. A prospective fate map from **one time point**.

**Structurally cannot ask:**
- Which cells passed *through* the HPCS at any past point? (HPCS at wk 8 → AT1-like at wk 14 is
  invisible unless you pulsed at wk 8)
- How long did a lineage spend in the HPCS?
- Did it enter more than once?
- What was the signalling environment at the moment of entry?

**Their own stated open question:** *"Future work will identify specific niche-derived inducers of
the HPCS in vivo."* Plus: *"As its physiological correlate is induced by epithelial injury,
regenerative niche signals may similarly induce the HPCS in a subset of cancer cells."*

> **A continuous cell-autonomous recorder is the instrument that question is missing. Cre gives one
> timestamp per animal; the tape gives every transition, ordered, in every lineage.**

### Which pathway to record

Their HPCS-enriched gene list is the shortlist: *"Tnf, Il23, Vegf, Cxcl2, Wnt7a, Wnt7b, Wnt10a,
Areg, Ereg and Tgfb1."*

**⚑ Strongest: NF-κB / inflammatory signalling.** The HPCS's physiological correlate is the
injury-induced DATP/PATS/ADI state, and the DATP literature they cite (ref. 52, Choi et al.,
*"Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors that mediate
alveolar regeneration"*) establishes **IL-1β/NF-κB as the inducer**. Mechanistic hypothesis + named
pathway + physiological precedent + (believed) existing ENGRAM NF-κB reporter. TNF/IL23/CXCL2 being
HPCS-enriched closes the loop — the state may both respond to and propagate inflammatory signalling.

**Runners-up:** TGF-β (Tgfb1 enriched; drives KRT8⁺ transitional states, EMT); WNT (three ligands
enriched; central to alveolar regeneration); hypoxia/HIF (Park's phylogenetic fitness tracked a
*"hypoxia-, glycolysis- and partial-EMT-centered axis"*).

### ⚑⚑ The two-channel design — record STATE as well as signal

**Drive an ENGRAM pegRNA off a cell-state-specific promoter — the *Slc4a11* promoter itself.**
Converts a transient state into a permanent, time-stamped record. Every HPCS episode written to the
tape, in every lineage, regardless of tamoxifen timing.

Enables: how many distinct HPCS episodes? how long each? do EMT-fated lineages differ from
alveolar-fated ones? **do ALL therapy-surviving lineages have an HPCS episode in their past?**

Pair with signal channel + constitutive reference ⇒ the question becomes an **ordering** question
*within a single tape*:

$$\underbrace{\text{NF-}\kappa\text{B symbol}}_{\text{position } i}\ \longrightarrow\ \underbrace{\textit{Slc4a11}\text{ symbol}}_{\text{position } j>i}\ ?$$

> **Was the niche signal on BEFORE the cell entered the high-plasticity state — in the same
> lineage?** As close to a causal ordering claim as a recorder can get, and it answers Chan et al.'s
> stated open question directly.

### ⚑ Highest-impact version: therapy resistance

Chan et al.: *"therapeutic responses primarily eliminate non-HPCS-derived cells, whereas
HPCS-derived cells become strongly enriched in MRD"* (cisplatin, MRTX1133).

> **Were surviving lineages already different before treatment, or did treatment induce the
> difference?** The pre-existing-vs-induced resistance question — the same one Park et al. built
> their paper on ("a pre-existing, heritable cell state"), which they could answer only with clonal
> barcodes + phylogeny, **with no access to signalling history. ENGRAM adds exactly that axis.**

### Practical: cell line, not GEMM

Chan et al. generated *"a clonal KP;Slc4a11^MCD/+ LUAD cell line with complete excision of the FSF
cassette"* in which the reporter *"faithfully marks the HPCS in subcutaneous allografts."*
Transduce with Typewriter tapes + ENGRAM cassettes; orthotopic allografts; harvest across
timepoints.

Exactly Park's proven route (NCI-H1299 xenograft); skips years of mouse engineering; and crucially
**the mScarlet reporter gives an independent ground-truth readout of HPCS identity at the tips** —
validate recorded history against directly measured state rather than going in blind. Move to the
autochthonous GEMM once the readout is trusted. ⚠ Use a **syngeneic** host, not NSG, if NF-κB /
inflammation is the signal — the immune compartment is part of the biology.

### ⚠ Three risks to design around

1. **⚠⚠ The HPCS is quiescent, and quiescence may slow editing.** *"the HPCS is highly quiescent"*
   (Ki-67 low) while derivatives proliferate rapidly. Prime editing is dNTP/cell-cycle sensitive
   (mouse embryo hematopoietic slowdown ← SAMHD1-limited dNTP supply). ⇒ **The state you most want
   to time may record most slowly**, systematically under-representing HPCS residency. **Row A2
   landing precisely on the variable of interest.** Measurable and correctable via the constitutive
   reference channel — but design for it.
2. **Timescale/saturation.** KP tumours take 12–30 weeks vs Park's 45 days. 4–6× longer ⇒ more
   capacity, longer tapes, or degraded late-stage resolution. Mulberry's framework gives the
   calculation *before* you start.
3. **The $q$ problem.** Strong induction ⇒ one symbol dominates ⇒ $\sum \xi_i^2$ rises ⇒ tree degrades
   exactly when the signal is most informative. **Keep a high-diversity constitutive background.**

## I.4 ⚑ What this means for the computational models

**1. Build the headline questions around ORDER, not TIME.** Pilarski: sequential editing
significantly improves topology, leaves timing unchanged. Mulberry: order is read directly from tape
position. ⇒ *"NF-κB preceded HPCS entry in this lineage"* needs only the **relative position of two
symbols in one tape** — no clock, no dating, no branch-length inference, robust to almost every
misspecification in the A-table. *"NF-κB peaked at day 12.3"* needs all of it.

> **Design the biology so the claim is an ordering claim, and the computational risk drops
> enormously.**

**2. The two-channel design makes the latent semi-supervised.** An ENGRAM channel driven by
*Slc4a11*, plus Slc4a11/mScarlet and full transcriptome measured at the tips ⇒ **the latent state is
observed at the leaves**, latent only internally. Best-identified version of the problem (§E.3) —
the reason to start there rather than with a purely latent signal.

**3. The tiers are good milestones because each has ground truth:** tier 1 validates $\xi(t)$ with no
latent field; tier 2 validates lineage assignment on a shallow tree with a known answer; tier 3
needs the full latent field, dropout model, and state-dependent rates.

## I.5 Open lead: Rudensky — regulatory T cells and TGF-β

**Status: raised by Jun Hong Choi; details unknown. Not yet worked out.** Placeholder with initial
thoughts to test against the actual proposal.

**Why it could be an unusually good ENGRAM fit — the *Foxp3* locus has textbook CREs.** ENGRAM's
design intent is recording *cis-regulatory element activity*, and Foxp3 has three well-mapped
conserved non-coding sequences with distinct, separable functions:

| Element | Function |
|---|---|
| **CNS1** | **TGF-β/Smad3-dependent**; required for *peripheral* Treg (pTreg) induction, dispensable for thymic |
| **CNS2 (TSDR)** | stability / heritable maintenance of Foxp3 expression |
| **CNS0** | Satb1-dependent, early priming of the locus |

⇒ Driving separate ENGRAM pegRNAs from CNS1 and CNS2 would record **induction** and **stabilization**
as *distinct, orderable symbols in the same tape.* That is a remarkably direct match to the
two-channel ordering design of §I.3–I.4.

**Candidate questions:**
- **Thymic vs peripheral origin.** A longstanding question with contested surrogate markers (Helios,
  Nrp1). A CNS1 record would settle it per-cell, retrospectively, without relying on markers.
- **Treg stability / "ex-Tregs."** Do Tregs lose Foxp3 in vivo, and how often? A permanent record of
  CNS2 activity would resolve a genuinely contested question.
- **Intratumoral Tregs: recruited or locally induced?** TGF-β is abundant in tumours; a TGF-β/CNS1
  record distinguishes cells that experienced tumour TGF-β from those that arrived already
  committed. Directly complementary to the §I.3 cancer application, and same institution.

**⚠ Serious technical caution.** The mouse embryo paper found **hematopoietic lineages edit
progressively more slowly as they differentiate**, attributed to PEmax silencing and/or
SAMHD1-limited dNTP supply. Lymphocytes may be a poor recording substrate. Check before committing —
this is row **A2** again, and it defeated recording in exactly this lineage. Possible mitigations:
transduced T-cell lines or adoptive-transfer systems rather than a germline GEMM; measure editing
rate in Tregs directly as a pilot.

**Also note:** T cells proliferate rapidly and are trivially isolatable by FACS, which is favourable
for both tree depth and cell recovery — *if* the editing rate holds up.

---

## Next up

**Session 2b — the BEAST 2 package specifics (p. 12).** BEAST 2 integration:
`SciPhySubstitutionModel`, `SciPhyTreeLikelihood` extending `GenericTreeLikelihood`,
caching, tape-level parallelism. What the MCMC actually samples (tree topology, branch
lengths, $\lambda$, $\xi$, birth–death parameters) and which operators move it.

**Session 3 — Validation and phylodynamics (pp. 13–15).** Well-calibrated simulation and
coverage, CCD summary trees, PI / wRF tree distances, the benchmark against UPGMA /
ordering-aware UPGMA / TiDeTree, and the birth–death-sampling priors in Tables 2–4
(including the piecewise-constant / Ornstein–Uhlenbeck skyline used for the gastruloid).

**Later — the extension work proper.** (i) Numerically verify the $H_m$ recursion in
§1c.7 against simulation. (ii) Identifiability of $\xi(t)$ from $N=6$ ordered positions.
(iii) Decide the SBI/likelihood boundary using the A1–A9 table.

---

# Session H.6 — Mulberry & Stadler, close read (Eqs. $f_k$, 1–4, Results 1)

Slow read of §2–3.1.1, adopting their notation. This is the foundational-theory task
flagged as thread #1. Kept in SciPhy register (§1a–1c). Numerics in this section were
checked against Monte Carlo (matches to 3 dp) before being written down.

## H.6.0 ⚠⚠ Notation — the house register, and the two source papers

**As of 2026-08-28 this repository uses one notation throughout: the HOUSE column.** The other two
columns exist only for reading the papers, whose symbols clash with each other and with ours.

| Concept | **HOUSE** | SciPhy paper | Mulberry & Stadler |
|---|---|---|---|
| sites per tape | $N$ | $N$ | $k$ |
| tapes per cell | $k$ | $k$ | $m$ |
| alphabet size (design) | $M$ | $M$ | $j$ |
| alphabet size (observed distinct) | $M_{\rm obs}$ | — | — |
| insert probabilities | $\xi=(\xi_1,\dots,\xi_M)$ | $\eta=(f_1,\dots,f_M)$ | $\xi$ |
| collision / homoplasy | $q=\sum_i\xi_i^2$ | $\sum_i f_i^2$ (buried in Eq. 9) | $q$ **(headline param)** |
| **editing rate** (per unit time) | $\lambda$ | $r$ | $\lambda$ |
| **birth–death pair** | $(b,\delta)$ | $(\lambda,\mu)$ | — |
| elapsed time | $\Delta t$ | $\Delta t$ | depth diff $d$; tips at depth 1 |
| expected edit count over an interval | $\lambda\,\Delta t$ | $r\,\Delta t$ | $r$ |
| edits realised | $\lvert c\rvert$ | $\lvert c\rvert$ | $x$ |

Three traps when reading the sources:

- ⚠ **Mulberry's $k$ and $m$ are the reverse of ours**: their $k$ = our $N$ (sites); their $m$ =
  our $k$ (tapes).
- ⚠ **Mulberry's $r$ is not a rate.** It is a *dimensionless expected count*: their $r$ = our
  $\lambda\Delta t$; their $\lambda$ = our $\lambda$. Their Table 1 flags this by calling $\lambda$
  the rate "**scaled by experimental duration**" — depths live in $[0,1]$, so $\lambda d$ =
  "expected edits by depth $d$". **Never let a bare $r$ stand in these notes.**
- ⚠ **$\delta$, not $d$, for the death rate**, because $d$ is already Mulberry's tree depth.

Mulberry's $q$ **is** SciPhy's $\sum_i f_i^2$ from §1c.2 / line 944 — the same object, promoted to
a first-class parameter. Everything about $\xi$ enters the theory through $q$ and nothing else.

> ### ⚑ Register decision (2026-08-28) — $\xi$ for insert probabilities, $f$ freed for functions
>
> Insert probabilities are $\xi$, and $q=\sum_i\xi_i^2$. **$f$ is reserved for genuine function
> symbols** — Mulberry's censored-Poisson pmf $f_k(x,\lambda d)$, and trajectory functions like
> $g_i$, $a_i$. The reason is not cosmetic: the entire ENGRAM extension is the claim that the
> insert probabilities *become functions of time and signal*, so it must be possible to write
> $\xi_i(t)$ and $\xi_i(\text{signal})$ without the symbol doing double duty.
>
> **The sweep through §0–§I was completed 2026-08-28**: $f_i\to\xi_i$, $\eta\to\xi$,
> $r\to\lambda$ (editing rate), $(\lambda,\mu)\to(b,\delta)$ (birth–death). Legacy $f_i$/$\eta$
> now appear only inside the SciPhy column of the table above.
>
> ⚠ **$\xi$ is a probability, not a rate.** It is *which symbol* is written given that an edit
> occurs; the rate at which edits occur is $\lambda$. The two are orthogonal, and §1a.3's
> jump-chain decomposition is what establishes that: the chain draws Categorical($\xi$)
> *independent* of the exponential holding time.
>
> Other live uses of these letters, all clearly scoped and deliberately left alone: $\lambda_a$ =
> CTMC exit rate from state $a$ (§1a.3, and it equals $\lambda$ for every non-absorbing state);
> $\mu$ = OU process mean (§G.2); Pagel's $\lambda$ (§G.3); and Park et al.'s own $\lambda$ in
> §1629, which is their organ-parsimony weight and their per-clone growth rate.

## H.6.1 $f_k(x,r)$ — the censored Poisson (their §2)

$$f(x,r)=\frac{r^x e^{-r}}{x!},\qquad
f_k(x,r)=\begin{cases} f(x,r) & 0\le x<k\\ 1-\sum_{x'=0}^{x-1}f(x',r) & x=k\\ 0 & x>k\end{cases}$$

**Identical to SciPhy Eqs. 3–5 / §1a.4 bounded Poisson**, with $r\mapsto\lambda d$. They cite
Seidel et al. for it. Let $N\sim\text{Poisson}(r)$; the tape records $X=\min(N,k)$; $f_k$ is
the pmf of $X$. Term by term:

- **Body** ($0\le x<k$): plain Poisson. Justified by constant exit rate = sequential
  architecture (one active site ever) ⇒ genuine homogeneous Poisson (§1a.4).
- **Boundary** ($x=k$): $f_k(k,r)=1-\sum_{x'=0}^{k-1}f(x',r)=\Pr(N\ge k)$. Tail piled onto
  the last state. This is §1a.6's absorbing wall as a probability. **No renormalisation.**
- **Support** ($x>k$): structural zero, same status as irreversibility.

**Terminology fix to propagate to §1a.4:** their word is **right-censored**, not "truncated,"
and theirs is correct. Censoring = observe $\min(N,k)$, keep all mass (pile at $k$).
Truncation = condition on $N\le k$ and renormalise. Line 2 divides by nothing ⇒ censoring.
Biologically right: a saturated tape didn't stop the editor, it ran out of writable sites;
the excess-attempt information is genuinely lost, not absent.

**Properties worth having in hand:**

1. Proper pmf on $\{0,\dots,k\}$, exactly (body + boundary telescope to 1).
2. **Saturation in closed form:** $f_k(k,r)=\Pr(N\ge k)=P(k,r)=\texttt{gammainc}(k,r)$, the
   regularised lower incomplete gamma (Poisson–Erlang duality: "$\ge k$ events by time $\lambda$"
   $\iff$ "sum of $k$ i.i.d. $\text{Exp}(1)\le r$"). Compute this way, never by summation.
3. **Monotonicity = the paper's whole tension:** $\partial_r f_k(k,r)=f(k-1,r)>0$ (saturation
   $\uparrow$), $f(0,r)=e^{-r}\downarrow$ (silence $\downarrow$). Every "optimal rate" result
   is a fight between these two monotone curves, and both live in $f_k$. Empirical shadow:
   §0.5, site 1 96.8% edited vs site 5 19.7%.
4. **⚠ Censored Poissons are NOT additive.** Ordinary Poissons compose (Chapman–Kolmogorov,
   §1a.3) — censored ones don't. You cannot get the depth-$d_2$ count by composing depth-$d_1$
   with an increment. **This is why every downstream expression is a disjunction over
   "saturated upstream, or not," instead of a clean product along a path.** (§1a.4 flagged the
   same: "truncation breaks additivity.")
5. **Marginal, not conditional.** SciPhy's $P_{a,b}(\Delta t)$ conditions on parent state $a$.
   Mulberry's $f_k(x,\lambda d)$ is marginal from an unedited origin — allowed because the
   questions are absolute ("saturated by depth $d$?"), not transitional. Buys simplicity,
   costs Markov composition (property 4).
6. **Purely the clock.** No $\xi$, no $j$. Symbols enter *only* later, *only* through $q$.
   This is the structural reason for §H.1(b): "$\xi$ affects accuracy only through $q$" is not
   empirical, it falls out of clock/jump independence (§1a.8). And it is why §H.1(c)'s tension
   is unavoidable: $q$ is the *sole* channel from $\xi$ to topology, and $q$ is minimised by
   uniform $\xi$ — ENGRAM peaks $\xi$ ⇒ no reparameterisation escapes it.

## H.6.2 Eq. (1) — the distance $D_{a,b}=\sum_{i=1}^m (k-k_i')$

$k_i'$ = largest index such that $S^i_{k''}(a)=S^i_{k''}(b)\ne 0$ **for all** $k''\le k_i'$
= **length of the longest common prefix** of the two tapes, $\lvert\mathrm{lcp}(a^{(i)},b^{(i)})\rvert$.

Three conditions fused in the clause, all load-bearing:
1. same character (not just both edited);
2. **$\ne 0$** — matching *zeros* is not sharing (every tape starts $0^k$; shared inactivity
   ≠ common ancestry);
3. **prefix closure** — a match behind an earlier mismatch does not count.

**Fig. 1 (right), worked:** A=`5 2 3 2 · ·`, B=`5 2 4 2 6 ·`. $k_1'=2$ (sites 1–2 agree; site
3 diverges 3 vs 4). **Site 4 agrees (2=2) but is NOT counted** — condition 3 — caption's
"non-sequential shared edits are not counted." That site-4 match is homoplasy the statistic
correctly refuses to credit. Contributes $k-k_1'=6-2=4$.

**Representation note:** SciPhy state = ordered tuple $\le N$ (§1a.7); Mulberry state =
fixed-length-$k$ vector over $\{0,\dots,j\}$, zero-padded. Equivalent *because* sequential
(zeros always a contiguous suffix). Theirs indexes sites cleanly; SciPhy's does prefix algebra.

**Four properties:**
- **(a) Not a metric:** $D_{a,a}=\sum_i(k-x_i)>0$ unless every tape saturated. $D$ measures
  *unshared recording capacity*, not divergence. Don't use metric intuitions.
- **(b) Depends on tree only via MRCA depth $d_v$:** with $N_i\sim\text{Poisson}(\lambda d_v)$,
  $D_{a,b}=\sum_i(k-\min(N_i,k))$. Tip depths (all 1) never appear. **This is where $f_k$
  enters:** each summand is $k-X$, $X\sim f_k(\cdot,\lambda d_v)$.
- **(c) Expected distance (not in paper, MC-verified):**
  $\mathbb{E}[D]=m\sum_{y=0}^{k-1}F_{\text{Pois}}(y;\lambda d_v)$, strictly decreasing $mk\to 0$.
  For $k=6$: $\approx k-\lambda d_v$ (one distance unit per expected edit) until $\lambda d_v\gtrsim k/2$,
  then collapses toward 0 as saturation destroys discrimination. The informative window is the
  middle; deep vs shallow MRCAs become indistinguishable once tapes fill.
- **(d) Homoplasy only *shrinks* $D$** (spurious prefix extension) ⇒ systematic bias making
  cells look **more** related, acting w.p. $q$ per site. Directional, not noise — can invert
  Eq. (2), not just blur it.

**⚑ Relation to SciPhy's lcp (third sighting, sharpened):** $\mathrm{lcp}(a^{(i)},b^{(i)})=
u^{(i)}_{\mathrm{MRCA}}$, the **maximal** element of $A_i=\mathrm{pre}(u_i)$ (§1c.2). So
$$D_{a,b}=mk-\big(\text{parsimony-reconstructed edit content at the MRCA}\big).$$
The paper = **SciPhy's ancestral-state machinery run in point-estimate (parsimony) mode**, so
accuracy is analytically tractable. Price: homoplasy must be modelled by hand (via $q$) instead
of integrated over. "Parsimony picks the top of the set; likelihood integrates over it" (line 951).

## H.6.3 Eq. (2) — resolution criterion $D_{a,b}<\min(D_{b,c},D_{a,c})$

$(a,b|c)$ = the rooted triple the true tree displays. With $v=\mathrm{MRCA}(a,b)$,
$u=\mathrm{MRCA}(a,b,c)$, $d_v>d_u$:

- **Green path** origin→$u$: edits shared by all three ⇒ cancel from Eq. (2).
- **Yellow path** $u\to v$, length $\ge\ell$: edits shared by $a,b$ only ⇒ **the sole source
  of signal.** ⇒ **This is why $\ell$ is *the* resolution parameter.**
- $(a,c)$ and $(b,c)$ both have MRCA $u$ ⇒ equal in distribution.

Holds in expectation automatically ($d_v>d_u$); the question is a *given realisation*.
- **$\min$**: must beat both alternatives (tie with one = unidentified triple).
- **strict $<$**: a tie is failure. With $q=0$ the only failure is a tie, occurring iff **zero
  edits land on the yellow path across all $m$ tapes** — the bridge to Eq. (3)/Result 1.

**It's the strict ultrametric three-point condition.** $\mathbb{E}[D]$ is genuinely ultrametric
(decreasing in MRCA depth) ⇒ **UPGMA is *correct here*, not just popular**: its molecular-clock
assumption is guaranteed by design (simultaneous harvest, constant $\lambda$). ⚠ **Fine print
for the extension:** ultrametric $\mathbb{E}[D]$ needs $\lambda$ shared across *lineages*.
Time-varying-but-lineage-common rate (§3.4; ZGA burst; **A1**) preserves it. State/lineage-varying
rate (**A2** hematopoietic slowdown, **A8** PE copy number) **destroys** it and breaks UPGMA
consistency. **A1 survivable, A2/A8 not** — not stated in the paper, directly relevant to us.

**"Every triplet identifiable" = the key logical move:**
- **Aho et al. 1981 (BUILD):** a consistent set of rooted triples ⇒ unique tree, poly-time.
  Converts global "P(whole tree right)" into local "P(this triple right)" = a one-edge calc.
- **Sufficient, NOT necessary:** real distance methods average over pairs and beat triplets ⇒
  bounds are **conservative lower bounds** (§H.1: NJ/UPGMA outperform triplet score). Note two
  *separate* loosenesses: (i) triplet reduction — conservative, structural; (ii) normal approx
  in Eq. (9) — approximate, either direction. Their "not rigorous bounds per se" = (ii).

**Two failure modes ⇒ two results:**
| Mode | Mechanism | Effect | Governed by | Bound |
|---|---|---|---|---|
| **Tie** | no edit on yellow path $(u,v)$, any tape | $D_{a,b}=D_{a,c}$ | $\lambda,\ell,d,m$ | $B_\infty$ §3.1.1 |
| **Inversion** | identical subseq arises independently on $(u,c)$ | $D_{a,c}\downarrow$ below $D_{a,b}$ | $q$ | $B_q$ §3.1.2 |

Setting $q=0$ deletes the inversion row ⇒ $B_\infty$ derivable in half a page.

## H.6.4 ⚠ Assumptions smuggled into Eqs. (1)–(2) — A-table candidates

1. **All $m$ tapes share $k$ and $\lambda$** (SciPhy has per-tape $r_z$, $N_z$). Heterogeneity
   turns every $(\cdot)^m$ into $\prod_i(\cdot)$ — mild, but §3.3 "how many tapes" answers
   become tape-set-specific.
2. **⚑⚑ Tape independence given tree (A5) — NOT benign.** The whole framework's power is the
   $m$-th power $p_0^m$. Per-cell frailty (§0.4 L3: shared PE dosage, cell-cycle, pseudo-
   processivity) makes tapes positively correlated ⇒ $\Pr(\text{all }m\text{ silent})>p_0^m$.
   **⇒ their bounds are anti-conservative in the one direction that matters** — every other
   looseness is conservative. Embryo #6 (no PEmax ⇒ effective $m=0$ despite 5 TAPEs) is the
   extreme. **Publishable-scale criticism; lands on the parameter §3.3's design advice is built
   from.**
3. **No dropout (A6):** $S^i_{k'}(v)$ always observed; ~50% of mouse matrix missing. Eq. (1)
   needs a missing-tape rule (drop → breaks cross-triple comparability of $D$; impute → not
   innocuous). Unaddressed.
4. **Ultrametric:** fine by design (§1c.1) but load-bearing — it's what legitimises UPGMA.
5. **$0$ = unedited conflated with unobserved** if real data fed in (same trap as line 920).

## H.6.5 Eq. (3) — $p_0$, the no-shared-edit-on-an-edge probability

$$\Pr(X^i_{u,v}=0)=\underbrace{f_k(k,\lambda d_u)}_{\text{saturated at }u}
+\underbrace{\big(1-f_k(k,\lambda d_u)\big)}_{\text{not saturated at }u}\underbrace{f(0,\lambda(d_v-d_u))}_{\text{no edit on }(u,v)}
\;\le\; f_k(k,\lambda d)+\big(1-f_k(k,\lambda d)\big)f(0,\lambda\ell)=:p_0(k,\ell,\lambda,d)$$

$X^i_{u,v}$ = shared edits on tape $i$ along edge $(u,v)$. It is **0 in exactly two disjoint
ways** — and the disjunction is *forced by non-additivity* (H.6.1 prop. 4), it is not a
modelling choice:
- **already saturated at $u$:** no capacity left ⇒ $f_k(k,\lambda d_u)=\texttt{gammainc}(k,\lambda d_u)$;
- **not saturated, but no event fires on the edge:** $(1-f_k(k,\lambda d_u))\cdot e^{-\lambda(d_v-d_u)}$
  (uses $f(0,\cdot)=e^{-\lambda\cdot}$).
Mutually exclusive and exhaustive ⇒ just add. (MC-verified: e.g. $k{=}6,\lambda{=}10,d_u{=}0.3,
d_v{=}0.5$ ⇒ analytic 0.2079, MC 0.2082.)

**The two bounding substitutions** (make it depend on $(k,\ell,\lambda,d)$ only, monotonically
worsening the estimate):
- $d_u\le d$ (**max depth of internal-split starts**) ⇒ raises the saturation term. ⚠ higher
  saturation ⇒ **higher** $p_0$ (a *full* tape can't record the split) — this is why $d$ is
  taken as a *max*, the pessimistic direction.
- $d_v-d_u\ge \ell$ (**min branch length**) ⇒ $e^{-\lambda(d_v-d_u)}\le e^{-\lambda\ell}$,
  raises the no-event term.
Both push $p_0$ **up** ⇒ $p_0$ is an **upper bound on the per-tape non-resolution prob**
(MC-verified: true 0.21 ≤ Eq3 0.61, etc.). Paper notes it "is not tight even" under synchronous
dense sampling, because it still uses max depth $d$.

**Then the $m$-th power (⚠ the A5-fragile step):** tapes independent ⇒
$\Pr(\text{no edit on }(u,v)\text{ on *any* of }m\text{ tapes})\le p_0^m$. With $q=0$,
irreversibility means *one* edit anywhere on the yellow path resolves *every* triple through $u$
(the edit is inherited by all of $u$'s descendants below $v$, absent in $c$'s subtree) ⇒
$p_0^m$ = **upper bound on P(split at $u$ not resolved)**.

**⇒ the informative-window story, now at edge level:** $p_0$ is U-shaped in $\lambda$. Small
$\lambda$: $e^{-\lambda\ell}\to1$, no edits, ties. Large $\lambda$: $f_k(k,\lambda d)\to1$, tapes
saturate before the split, also fails. Minimised at intermediate $\lambda$ — same window as
H.6.1(3)/§H.1(a). $m$ deepens the well ($p_0^m$) but doesn't move its floor in $\lambda$.

## H.6.6 Result 1 ($B_\infty$) — the first reconstruction bound

$$\boxed{\;\Pr(\text{exact reconstruction})\ \ge\ 1-(n-2)\,p_0(k,\ell,\lambda,d)^{\,m}\;}\qquad(q=0)$$

**Derivation = union bound over internal splits.**
- Tree on $n$ tips, ultrametric binary ⇒ $n-1$ internal nodes (§1c.1), of which the origin/root
  structure leaves **$n-2$ internal splits** whose resolution is at risk (traverse from the
  stem terminus, union-bound over the $n-2$ internal branches — they take a union bound because
  they can't assume independence between different splits).
- Each unresolved w.p. $\le p_0^m$ (H.6.5).
- $\Pr(\text{any split fails})\le (n-2)p_0^m$ (Boole) ⇒ complement is the bound.

**Reading it:**
- **Only $q=0$** — no homoplasy. Result 2 ($B_q$) restores it. $B_\infty$ is the ceiling: best
  case, infinite alphabet / uniform $\xi$.
- **⚠ Vacuous when $(n-2)p_0^m\ge1$** — goes negative (MC table: $m{=}10$ ⇒ bound $-6.4$;
  $\lambda{=}12$ oversaturates ⇒ $p_0=0.97$, bound $-50$). A union bound over $n-2$ terms is
  only useful once $p_0^m$ is driven below $1/(n-2)$, i.e. **$m$ must clear a threshold set by
  $n$** before *any* guarantee exists. This is the mathematical form of §H.1(d): the constraint
  is $m$ (tapes) and $\ell$ (resolution) — $n$ enters only weakly, as $\log n$, via the
  threshold $p_0^m<1/(n-2)$.
- **Design use (their §3.3):** invert it. Fix accuracy $c=1-\epsilon$, pick $\lambda^*$
  minimising $p_0$, solve $m^*=\lceil\log((1-c)/(n-2))/\log p_0(k,\ell,\lambda^*,d)\rceil$.
  This is where "$k{=}5$ ⇒ ~30 tapes for 512 tips at $\ell{\approx}0.1$" comes from (§H.1(d)).
- **Sensitivity ordering:** through $p_0^m$, sensitivity to $\lambda,\ell,m$ is exponential;
  to $n$, logarithmic. ⇒ **resolution $\ell$ and tape count $m$ dominate; sample size $n$ barely
  matters.** Reframes the scale question exactly as §H.1(d) already logs it.

**For us:**
- The $q=0$ ceiling is the *good* case; ENGRAM pushes $q$ up (§H.1(c)) so real performance sits
  below Result 1 and needs Result 2.
- The A5-fragility (H.6.4 #2) means the true exponent is $<m$, so even Result 1 is optimistic on
  real data with per-cell frailty. **Both corrections point the same way: Result 1 over-promises
  on ENGRAM data.** Quantifying that gap is a concrete contribution.

## H.6.7 Eq. (5) — $h_{k'}$, the homoplasy engine (their §3.1.2)

**Where we are.** §3.1.2 turns on $q>0$. The new failure mode is *inversion*: $c$ independently
acquires the same leading edits as $a,b$, shrinking $D_{a,c}$ below $D_{a,b}$ (H.6.2(d)). Eq. (5)
is the single primitive that quantifies it: **how long a shared prefix arises purely by chance
on two independent paths.** Everything $q$-dependent in the rest of the paper is built from $h$.

### The setup display

Two paths, $n_a$ and $n_c$ edits acquired independently. Given $n_a,n_c$, the probability of an
identical sequence of length $x$ from the first position:
$$q^x \ \text{(if the shorter path has no site } x{+}1)\qquad\text{vs}\qquad q^x(1-q)\ \text{(otherwise)}$$
$q^x$ = the first $x$ characters collide, i.i.d. w.p. $q$ each. The $(1-q)$ = **position $x+1$
must *differ***, or the shared prefix would be longer than $x$. This is $\Pr(\text{prefix}=x)$
exactly, not $\ge x$.

⚠ **Minor inconsistency in the paper.** The setup display gives the condition as "$q^x$ if
$n_a=k$"; the sentence after Eq. (5) gives it as "multiply by $1-q$ only if the length of the
**shorter path is strictly greater than $x$**", and Eq. (5)'s exponent is $\mathbf 1(i>x)$.
**The latter two are operative and mutually consistent; read the first display as the saturation
edge case ($x=n_a=k$).** The correct condition is $\min(n_a,n_c)>x$: you only need position
$x+1$ to differ if position $x+1$ *exists on both paths*. If the shorter path stopped at exactly
$x$ edits, the prefix ends by exhaustion, no $(1-q)$ needed.

### Equation (5) as printed

$$h_{k'}\big(x,\lambda(1-d),q\big)=\sum_{i=x}^{k'}\Big[f_{k'}(i,\cdot)^2\,q^x(1-q)^{\mathbf 1(i>x)}
\;+\;2f_{k'}(i,\cdot)\,q^x(1-q)^{\mathbf 1(i>x)}\!\!\sum_{j>i}^{k'}\!f_{k'}(j,\cdot)\Big]$$

**Arguments:** $k'$ = *remaining* capacity (sites still writable after what was consumed
upstream); $\lambda(1-d)$ = expected edits on a path from depth $d$ **down to the tips at depth 1**
— the two paths in question are $v\to a$ and $v\to b$, or $u\to c$, all running to the present.
(Upper limit printed as $k$ in the paper; harmless, since $f_{k'}(i,\cdot)=0$ for $i>k'$.)

**The two terms** are the paper's hand-decomposition of "min of two i.i.d. counts":
- first: both paths have **exactly $i$** edits ⇒ $f_{k'}(i)^2$;
- second: shorter has $i$, longer has $j>i$ ⇒ $f_{k'}(i)f_{k'}(j)$, **doubled** for the two
  orderings (which path is the short one).
The $(1-q)$ exponent uses the **shorter** length $i$ in both terms. ✓ consistent with the prose.

### ⚑ Eq. (5) collapses — much cleaner than printed (verified)

Let $S_i=\Pr(\text{len}\ge i)$, $M=\min(n_a,n_c)$ with $n_a,n_c\overset{iid}\sim f_{k'}(\cdot,r)$.
The bracket is $\xi_i^2+2\xi_i\sum_{j>i}\xi_j=\xi_i(2S_i-\xi_i)=S_i^2-S_{i+1}^2=\Pr(M=i)$. Hence

$$\boxed{\;h_{k'}(x,r,q)\;=\;q^x\Big[\Pr(M=x)\;+\;(1-q)\,\Pr(M>x)\Big]\;}$$

Equivalently, and this is the real content:
> **$X\mid M \sim$ Geometric$(1-q)$ censored at $M$**, where $M=\min$ of two i.i.d. censored
> Poisson path lengths. **Homoplasy is a geometric race — each site independently "fails to
> collide" w.p. $1-q$ — censored by whichever path runs out of capacity first.**
> Same censoring pattern as $f_k$ itself (H.6.1), one level up. The double sum is nothing but
> the textbook min-of-two-iid identity written out by hand.

**Verified numerically** (paper form vs closed form vs Monte Carlo of the generative process,
agreement to 4 dp across $k'\in\{4,6,9\}$, $r\in\{1,3,5,8\}$, $q\in\{1/64,1/16,1/4\}$):

| $k'$ | $\lambda$ | $q$ | $h(0)$ | $h(1)$ | $h(2)$ | $h(3)$ |
|---|---|---|---|---|---|---|
| 6 | 3 | 1/4 | 0.7743 | 0.1856 | 0.0349 | 0.0047 |
| 6 | 3 | 1/16 | 0.9436 | 0.0539 | 0.0024 | 0.0001 |
| 6 | 8 | 1/4 | 0.7502 | 0.1877 | 0.0469 | 0.0116 |
| 9 | 5 | 1/64 | 0.9846 | 0.0152 | 0.0002 | ~0 |

**$h$ is a proper pmf over $x\in\{0,\dots,k'\}$, summing to exactly 1** (conditional on $M$,
the censored geometric sums to $\sum_{x<M}q^x(1-q)+q^M=1$). Useful invariant for implementation.

### What it says quantitatively

Mean spurious prefix (uncensored limit) $=q/(1-q)$:

| $q$ | $1/64$ | $1/16$ | $1/8$ | $1/4$ |
|---|---|---|---|---|
| $\mathbb{E}[\text{spurious shared sites}]$ | 0.016 | 0.067 | **0.143** | **0.333** |

⇒ **The competition that decides $B_q$ vs $B_\infty$:** the yellow path contributes $\approx\lambda\ell$
*genuine* shared edits; homoplasy contributes $\approx q/(1-q)$ *spurious* ones to the out-group.
**When $\lambda\ell\lesssim q/(1-q)$, homoplasy wins and reconstruction fails.** This is the
mechanism behind §H.1's "$B_\infty$ significantly overestimates accuracy under low rates at
$q=1/4$" — at low $\lambda$ the real signal shrinks toward the noise floor, which is *fixed* by $q$.
Recall §1c.2 line 944: $q=\sum_i \xi_i^2\ge 1/j$ (Cauchy–Schwarz, equality iff uniform $\xi$), so
the mouse's $j=8$ floor is $q\ge0.125$ ⇒ **≥0.14 spurious sites per comparison, unavoidable.**

> ⚑ **For ENGRAM this is the sharp form of §H.1(c).** Fine time-resolution means small $\ell$;
> strong signal means peaked $\xi$ means large $q$. The two design goals push $\lambda\ell$ and
> $q/(1-q)$ toward each other from opposite sides. **The feasible region is where
> $\lambda\ell \gg q/(1-q)$** — that single inequality is the cleanest one-line summary of the
> experimental-design constraint, and it is not stated anywhere in the paper.

### Where $h$ gets used (Eqs. 6–7, next session)

- **Eq. (6), out-group** $\Pr(X^i_{\text{out}}=x)=\sum_{y'=0}^{k}f_k(y',\lambda d_u)\,
  h_{k-y'}(x-y',\lambda(1-d_u),q)$ — $y'$ genuine shared edits above $u$ (real common ancestry),
  then $x-y'$ **purely spurious** ones below, with capacity reduced to $k-y'$.
- **Eq. (7), in-group** adds the middle factor $f_{k-y'}(y,\lambda\ell)$ — the yellow path — before
  applying $h$ from $v$ down. **The in/out asymmetry is exactly that one extra term.** That is
  the entire signal, and Result 2 is the probability it survives the noise in $h$.

## H.6.8 Eqs. (6)–(8) — the two counts and the contest between them

### The one idea

**Resolution is a contest between two numbers.** Since $D_{a,b}=\sum_i(k-k_i')$ and $k$ is
constant, Eq. (2) $D_{a,b}<D_{a,c}$ is *identically* the statement
$$\sum_i k_i'(a,b)\;>\;\sum_i k_i'(a,c)\qquad\Longleftrightarrow\qquad \sum_i X^i_{\text{in}}>\sum_i X^i_{\text{out}}.$$
$X_{\text{in}}$ = shared prefix length between $a,b$; $X_{\text{out}}$ = between $a,c$. **Eq. (6)
is the pmf of one, Eq. (7) the pmf of the other, Eq. (8) says who wins.** Everything else in
§3.1.2 is bookkeeping. ⚠ The prose defining $X_{\text{in}}$ as edits "which arise *after* the
three lineages diverge at $u$" is misleading — the equations plainly include the pre-$u$ term
$y'$ in **both**. Trust the equations.

### The three zones and the budget

```
             origin ──── u ──── v ──── tips
  GREEN  origin→u   (len d_u)      : real, shared by a,b,c  → in BOTH counts, cancels
  YELLOW u→v        (len ≥ ℓ)      : real, shared by a,b    → in X_in ONLY  ⇐ THE SIGNAL
  BELOW  v→a,v→b / u→c (to depth 1): private; agreement only by chance → homoplasy, h()
```
$k$ sites is a **budget spent top-down**: green consumes $y'$, yellow consumes $y$, whatever
remains is available for spurious agreement below. Hence the capacity cascade
$k\to k-y'\to k-y-y'$.

### Eq. (6) — out-group: real + fake

$$\Pr(X^i_{\text{out}}=x)=\sum_{y'=0}^{k} \underbrace{f_k(y',\lambda d_u)}_{y'\text{ green edits}}\;
\underbrace{h_{k-y'}\big(x-y',\,\lambda(1-d_u),\,q\big)}_{x-y'\text{ spurious, remaining capacity}}$$

$a$ and $c$ share **two** things: the genuine green prefix, and a chance extension. Both their
post-$u$ paths run to depth 1 ⇒ both of length $1-d_u$ ⇒ $h$'s two-i.i.d.-paths setup is exactly
right. Marginalise over $y'$.

### Eq. (7) — in-group: real + **more real** + fake

$$\Pr(X^i_{\text{in}}=x)=\sum_{y=0}^{x}\sum_{y'=0}^{k}
f_k(y',\lambda d_u)\;\underbrace{f_{k-y'}(y,\lambda\ell)}_{\textbf{THE SIGNAL}}\;
h_{k-y-y'}\big(x-y-y',\,\lambda(1-d_v),\,q\big)$$

**Identical to Eq. (6) except for one inserted factor.** Side by side:

| Zone | Eq. (6) out-group | Eq. (7) in-group |
|---|---|---|
| green (real) | $f_k(y',\lambda d_u)$ | $f_k(y',\lambda d_u)$ — same |
| **yellow (real)** | **absent** | $f_{k-y'}(y,\lambda\ell)$ ⇐ **the entire signal** |
| below (fake) | $h_{k-y'}(x-y',\lambda(1-d_u),q)$ | $h_{k-y-y'}(x-y-y',\lambda(1-d_v),q)$ — later start, less capacity |

⇒ **The whole theory reduces to: does one extra $f(\cdot,\lambda\ell)$ factor outrun the noise
in $h(\cdot,q)$?** That is $\lambda\ell$ vs $q/(1-q)$ again (H.6.7), now exact instead of heuristic.
$\ell$ is used in place of $d_v-d_u$ — the same conservative substitution as Eq. (3).

### ⚑ The chaining trick — non-additivity resolved

H.6.1(4) said censored Poissons **don't add**. But they **chain exactly if you decrement the
capacity**:
$$\min(P_1,k)\;+\;\min\big(P_2,\,k-\min(P_1,k)\big)\;=\;\min(P_1+P_2,\,k)\quad\text{(exact, both cases)}$$
So $f_k(y',\lambda d_u)\,f_{k-y'}(y,\lambda\ell)$ is not an approximation — it is the honest joint
law. **This is the same conditioning move as Eq. (3)'s disjunction, generalised.** Eq. (3)
conditioned on "saturated or not"; Eqs. (6)–(7) condition on *how much budget is left*. Same idea,
finer grain. Good mental model: **the tape is a budget, and every equation in this section is
"spend from the top, pass the remainder down."**

### Eq. (8) — the contest

$$\Pr\Big(\sum_{i=1}^m X^i_{\text{in}}-\sum_{i=1}^m X^i_{\text{out}}\le 0\Big)$$
$=\Pr(\text{triplet not resolvable})$. Note **$\le$**, not $<$: a tie is a failure (Eq. 2 is strict).
Exact via $m$-fold convolution when $mk$ small; hence the normal approximation in Eq. (9).

### Verification (Eqs. 6–7 vs direct triplet simulation)

$k{=}6,\lambda{=}8,d_u{=}0.35,d_v{=}0.55,q{=}1/4$, 400k sims. Both proper pmfs (sum $=1.000000$):

| $x$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Eq. 6 out | .0458 | .1397 | .2144 | .2210 | .1721 | .1081 | .0991 |
| MC out | .0453 | .1400 | .2154 | .2200 | .1716 | .1082 | .0996 |
| Eq. 7 in | .0094 | .0436 | .1014 | .1578 | .1844 | .1729 | .3306 |
| MC in | .0092 | .0437 | .1018 | .1573 | .1831 | .1730 | .3320 |

### ⚑⚑ NEW FINDING — $X_{\text{in}}$ and $X_{\text{out}}$ are strongly positively correlated, and Eq. (9) ignores it

Eq. (9) combines dispersions as $\sigma_{\text{in}}+\sigma_{\text{out}}$ with **no covariance term**
— i.e. it treats the two counts as independent. **They are not**: on the same tape they share the
same green-path realisation $y'$ (and are further coupled through $N_y$, $N_a$ and the shared
capacity budget). Measured at the parameters above:

$$\text{corr}(X_{\text{in}},X_{\text{out}})=\mathbf{0.68},\qquad
\mathrm{Var}(X_{\text{in}}-X_{\text{out}})=1.65\ \text{true}\ \text{vs}\ 5.19\ \text{assuming independence}\ (\mathbf{+215\%})$$

**Direction:** $\mu_{\text{in}}-\mu_{\text{out}}>0$ always (the yellow term), and for $Z\sim N(\mu,\sigma^2)$
with $\mu>0$, $\Pr(Z\le0)=\Phi(-\mu/\sigma)$ **increases** with $\sigma$. ⇒ **inflated variance ⇒
inflated failure probability ⇒ the bound is conservative.** Magnitude is large:

| $m$ | 1 | 3 | 5 | 10 | 20 | 30 |
|---|---|---|---|---|---|---|
| TRUE $\Pr(\text{fail})$ | .3056 | .0630 | .0178 | .00105 | .00001 | ~0 |
| independence assumption | .3678 | .2060 | .1314 | .04996 | .00900 | .00181 |
| over-estimate factor | 1.2× | 3.3× | 7.4× | **48×** | **900×** | — |

> ⚑ **This explains an observation the paper reports but does not account for:** *"$B_q$ is much
> more robust, and we do not find any parameter regimes in which this bound is violated. It is,
> however, not particularly tight."* **The un-modelled positive covariance is the reason for
> both halves of that sentence** — never violated (conservative direction) and never tight
> (1–3 orders of magnitude at realistic $m$).

**Bookkeeping of the three correlation issues so far — they do NOT all point the same way:**

| Correlation | Modelled? | Direction of error |
|---|---|---|
| $X_{\text{in}}$ vs $X_{\text{out}}$ within a tape | ✗ ignored (Eq. 9) | **conservative** — bound too pessimistic, large effect |
| Between different splits/triplets | ✗ avoided by Boole | **conservative** — union bound |
| **Between tapes (row A5, per-cell frailty)** | ✗ assumed independent | **⚠ ANTI-conservative** — bound too optimistic |

⇒ For us: the paper's headroom is larger than it looks *and* its one unsafe assumption is the one
we care about. **Both are quantifiable.** A corrected $B_q$ that (i) keeps the covariance and
(ii) puts a frailty term on the tapes is a concrete, self-contained methods contribution.

## H.6.9 Eqs. (9)–(10) — the normal approximation and the worst-case depth

### Why a normal appears at all — nothing biological

**$m$ (tape count) is the sample size in a CLT.** Define the per-tape **margin**
$Z^i=X^i_{\text{in}}-X^i_{\text{out}}$ — one tape's "vote" on whether the triplet resolves.
Tapes are i.i.d. (A5), so Eq. (8) is
$$\Pr\Big(\sum_{i=1}^m Z^i\le 0\Big),\qquad Z^i\ \text{i.i.d., bounded on }[-k,k].$$
A sum of $m$ i.i.d. bounded integer rvs ⇒ CLT ⇒ approximately normal. **There is no Gaussian
anywhere in the editing process.** The normal is purely a computational stand-in for the $m$-fold
convolution of $Z$'s pmf. Nothing is being *modelled*; something is being *approximated*.

Measured skewness of $Z$ at $k{=}6,\lambda{=}8,q{=}1/4$: **$-0.25$** — mild, so Berry–Esseen error
$\sim\text{skew}/\sqrt m$ is small. The CLT is well-behaved here; that is not the problem.

### Notation: $\Phi(\alpha,\mu,\sigma)$

Normal CDF with mean $\mu$, dispersion $\sigma$, **evaluated at $\alpha$** — i.e.
$\Phi\big((\alpha-\mu)/\text{sd}\big)$. Unusual way to write it; $\alpha$ is the *evaluation point*,
not a significance level.

### The continuity correction $\alpha$

$S=\sum Z^i$ is **integer-valued**, so $\Pr(S\le0)=\Pr(S<1)$ and a continuous approximation
should be evaluated at the midpoint $\mathbf{0.5}$ (textbook continuity correction). They use
$\boldsymbol{\alpha=1}$ instead — deliberately overshooting, because $\mu_S>0$ so evaluating
further right *raises* the failure probability. Their stated reason ("$\alpha=0.5$ occasionally
leads an overestimate of the reconstruction probability for small $m$") is confirmed below.

### ⚠⚠ Two notation problems in Eq. (9) as printed

$$\Pr\Big(\sum X^i_{\text{in}}-\sum X^i_{\text{out}}\le0\Big)\approx\Phi\big(\alpha,\ \mu_{\text{in}}-\mu_{\text{out}},\ \sigma_{\text{in}}+\sigma_{\text{out}}\big)$$

1. **$m$ is missing.** $\mu_{\text{in}}$ is defined as the **per-tape** mean $E[X^i_{\text{in}}]$,
   but $S$ has mean $m(\mu_{\text{in}}-\mu_{\text{out}})$ and variance $m(\sigma^2_{\text{in}}+\sigma^2_{\text{out}})$.
   As literally printed Eq. (9) **has no $m$ in it** — impossible, since Eq. (10) declares
   $p^*_{\text{trip}}(k,\ell,\lambda,q,m)$. Verified: the printed form is constant in $m$
   (0.468959 at every $m$). **Typo; read in the $m$.**
2. **Variance or SD?** They define $Var[X^i]=\sigma^2$, so "$\sigma_{\text{in}}+\sigma_{\text{out}}$"
   literally adds **standard deviations** — wrong for independent sums (variances add). Adding SDs
   is *also* conservative ($\sigma_A+\sigma_B\ge\sqrt{\sigma_A^2+\sigma_B^2}$) but wildly loose.
   Numerics below say they must mean **variances**; the add-SD reading is off by ~11× at $m=30$.

### Numerical check — exact vs. the variants

$k{=}6,\lambda{=}8,q{=}1/4,\ell{=}0.2,d_u{=}0.35$. "EXACT" = $m$-fold convolution under the
paper's *own* independence assumption (so this isolates the CLT error only).

| $m$ | EXACT conv | $\alpha{=}1$, add **var** | $\alpha{=}0.5$ | $\alpha{=}1$, add **SD** | as printed (no $m$) |
|---|---|---|---|---|---|
| 1 | 0.36676 | 0.45615 | 0.37082 | 0.46896 | 0.46896 |
| 3 | 0.20533 | 0.24265 | 0.20481 | 0.31083 | 0.46896 |
| 5 | 0.13143 | 0.15108 | 0.12924 | 0.23279 | 0.46896 |
| 10 | 0.04996 | 0.05502 | 0.04771 | 0.12921 | 0.46896 |
| 20 | 0.00887 | 0.00918 | 0.00804 | 0.04769 | 0.46896 |
| 30 | 0.00174 | 0.00170 | 0.00150 | 0.01919 | 0.46896 |
| 50 | 0.000075 | 0.000066 | 0.000058 | 0.00344 | 0.46896 |

- **$\alpha{=}0.5$ dips below exact from $m\ge3$** ⇒ anti-conservative ⇒ exactly why they chose
  $\alpha=1$. Their stated justification checks out.
- **$\alpha{=}1$ + add-variances tracks exact well** and is conservative up to $m\approx20$.
  ⚠ **But it goes slightly anti-conservative at $m\gtrsim30$** (0.00170 vs 0.00174; 0.000066 vs
  0.000075 — a 12% *under*estimate of failure at $m=50$). Small in absolute terms, but the
  "conservative" claim is not uniform in $m$. Worth knowing before quoting the bound as safe.
- **Add-SDs** is ~11× too pessimistic at $m=30$ ⇒ almost certainly not the intended reading.

### ⚑ The approximation is a convenience, not a necessity

They write "*for small values of $mk$ we can compute this directly as an $m$-fold convolution;
however, in general an approximation is necessary.*" **Measured: the exact 30-fold convolution
takes 0.13 ms.** Support of $S$ is $2mk+1$ points; repeated `np.convolve` (or one FFT) is
$O(mk\log mk)$. For every parameter set in the paper this is instant.

> ⇒ **Replacing Eq. (9) with the exact convolution is a trivial code change that removes the
> CLT error entirely** — and it removes the $m$-scaling ambiguity, the $\sigma$-vs-$\sigma^2$
> ambiguity, and the continuity-correction fudge in one stroke. Cheap improvement; take it.

### ⚑ The irony that matters for us

The normal approximation **needs $m$ large**. §3.3's entire purpose is **finding the minimum
$m$**. So the approximation is least reliable exactly where the paper leans on it hardest — and
at $m=1$–3 the error is 10–25%. They concede this ("*a significant limitation … is that $m$ needs
to be large*"). ⇒ **Any $m^*$ we quote from $B_q$ at small $m$ should be recomputed exactly.**

### Eq. (10) — $p^*_{\text{trip}}$, worst case over depth

$$p^*_{\text{trip}}(k,\ell,\lambda,q,m)=\max_{d'\in[0,d]}\Phi\big(1,\ \mu_{\text{in}}-\mu_{\text{out}},\ \sigma_{\text{in}}+\sigma_{\text{out}}\big)$$

**Why a max.** $\mu_{\text{in}},\mu_{\text{out}},\sigma$ all depend on **where in the tree the
triplet sits** — i.e. on $d_u=d'$ — via Eqs. (6)–(7). Result 2 is a Boole union bound over *all*
$\binom n3$ triplets, and a union bound needs **one** number that upper-bounds every triplet's
failure probability. So: take the worst depth. Same logical move as Eq. (3)'s $d_u\le d$, now
done by explicit optimisation instead of a crude substitution.

⚠ **Wording bug:** the paper says "this **minimum** occurs at the maximal depth" / "this
**minimum** can occur higher in the tree", but Eq. (10) is a **max**. They mean the minimum
*reconstruction* probability = maximum *failure* probability. Read accordingly.

**Why the worst case is usually deepest.** Two competing effects as $d'$ increases:
- ↑ failure: more budget spent above $u$ ⇒ **less capacity left for the yellow path to write the
  distinguishing edits**. In the limit the tape saturates before $u$ and the signal is
  structurally unwritable. (This is the paper's "irreversible nature of edits and limited number
  of target sites.")
- ↓ failure: paths below are shorter ⇒ less runway for homoplasy.
Usually the first dominates ⇒ max at $d'=d$.

**Verified, both regimes** ($\ell=0.15$, $m=10$, scan $d'$ on a 0.05 grid):

| Regime | argmax $d'$ | max allowed | $p_{\text{trip}}$ at $d'=0$ → $0.75$ |
|---|---|---|---|
| A: $k{=}6,\lambda{=}8,q{=}1/16$ (fast, diverse) | **0.85 (boundary)** | 0.85 | 0.0013 → 0.276, monotone ↑ |
| B: $k{=}13,\lambda{=}2,q{=}1/4$ (slow, low-div) | **0.75 (interior)** | 0.85 | 0.257 → 0.385, turns over |

✓ Reproduces their claim exactly: **high $k$, high $q$, low $\lambda$ moves the worst case off the
boundary.** Mechanism: with $k{=}13$ and $\lambda{=}2$ (≈2 expected edits vs 13 sites) the budget
**never** exhausts, so the dominant ↑ effect is switched off and the homoplasy-runway term wins at
the deepest depths. ⇒ **You cannot just evaluate at $d'=d$; the optimisation is real.** Note this
also means Result 1's crude "take $d_u\le d$" substitution is not merely loose but can point at
the wrong depth entirely.

## H.6.10 Result 2 / Eq. (11) — $B_q$, and why it is so much weaker than $B_\infty$

$$\boxed{\;\Pr(\text{exactly reconstruct }\mathcal T\text{ on }n\text{ tips})\;\ge\;1-\binom{n}{3}\,p^*_{\text{trip}}\;}\qquad(0<q<1)$$

Boole over **all $\binom n3$ triplets**, each failing w.p. $\le p^*_{\text{trip}}$ (Eq. 10).
Structurally identical to Result 1 — only the *unit of failure* changed.

### ⚑ Why $\binom n3$ here but $(n-2)$ in Result 1

| | Result 1 ($q=0$) | Result 2 ($q>0$) |
|---|---|---|
| unit | $n-2$ internal **branches** | $\binom n3$ **triplets** |
| why | one edit on the yellow path resolves **every** triplet through that branch *at once* — the edit is deterministically inherited by all of $v$'s descendants and absent from $c$'s side | homoplasy is **triplet-specific**: whether $c$ accidentally matches depends on *which $c$*. Two triplets sharing the same $(u,v)$ can have different outcomes ⇒ **no collapse** |

**The collapse at $q=0$ is the whole reason Result 1 is tight.** Turning $q$ on destroys it.

$$\binom n3\Big/(n-2)=\frac{n(n-1)}{6}:\quad n{=}128\to \mathbf{2{,}709\times},\quad n{=}512\to \mathbf{43{,}605\times},\quad n{=}1000\to \mathbf{166{,}500\times}$$

**Vacuity thresholds** (bound is negative until the per-unit failure prob clears these):

| $n$ | Result 1 needs $p_0^m<$ | Result 2 needs $p^*_{\text{trip}}<$ |
|---|---|---|
| 128 | $7.9\times10^{-3}$ | $2.9\times10^{-6}$ |
| 512 | $2.0\times10^{-3}$ | $4.5\times10^{-8}$ |
| 1000 | $1.0\times10^{-3}$ | $6.0\times10^{-9}$ |

### ⚑⚑ NEW FINDING — the paper's explanation for non-tightness is the *minor* cause

Paper: *"$B_q$ … does not reduce to $B_\infty$ when $q=0$, **due to the fact that we must bound
over all triplets**."* That is true but accounts for only a constant factor $n(n-1)/6$.
**The dominant cause is that Eq. (9) estimates an extreme-tail probability with a CLT.**

At $q\to0$ the truth is exactly $p_0^m$ (triplet fails iff no yellow edit on any tape).
Decomposed at $\lambda^*{=}4.5$, $k{=}6$, $\ell{=}0.15$, worst depth $d'{=}0.85$:

| $m$ | TRUTH $p_0^m$ | indep. + **exact** convolution | indep. + **normal** (the paper) |
|---|---|---|---|
| 5 | $7.9\times10^{-2}$ | $3.4\times10^{-1}$ | $3.8\times10^{-1}$ |
| 10 | $6.2\times10^{-3}$ | $2.6\times10^{-1}$ | $2.9\times10^{-1}$ |
| 20 | $3.8\times10^{-5}$ | $1.7\times10^{-1}$ | $1.8\times10^{-1}$ |
| 30 | $2.4\times10^{-7}$ | $1.2\times10^{-1}$ | $1.3\times10^{-1}$ |
| 40 | $1.5\times10^{-9}$ | $8.6\times10^{-2}$ | $9.1\times10^{-2}$ |

Per-triplet inflation: **46× at $m{=}10$, $5.3\times10^{5}$ at $m{=}30$, $8\times10^{11}$ at
$m{=}60$** — and it **grows without bound in $m$**, whereas the $\binom n3$ prefactor is constant
in $m$. At $m=30,n=128$ the tail error beats the prefactor by ~200×.

**Root cause, and it is simple:** at $q=0$, $Z^i=(\text{yellow edits})\ \ge 0$ **always**. So
$\sum_i Z^i\le0$ iff *every* $Z^i=0$ — a pure extreme-tail event. **Both the independence
assumption and the normal put probability mass on negative values that can never occur.**
⚠ Correcting H.6.9's claim: the decomposition above shows **the independence assumption is the
dominant error** (0.12 vs $2.4\times10^{-7}$), and the normal adds only a little on top
(0.128 vs 0.120). **Fixing the covariance matters far more than replacing the CLT.**

> **The general principle worth carrying:** Eq. (8) asks for a **large-deviation** probability;
> the CLT is a **central**-limit statement. A CLT controls the centre with *absolute* error
> $O(1/\sqrt m)$ — but the target here is $\sim10^{-7}$. An absolute error of $0.1$ swamps a
> $10^{-7}$ target by six orders of magnitude. **You cannot estimate a rare-event probability
> with a CLT**; you need exact convolution (0.13 ms, H.6.9) or a Chernoff/large-deviation bound.

### Practical cost, at the optimal rate

$k{=}6,\ \ell{=}0.15,\ d{=}0.85$; $\lambda^*=4.50$ (minimises Eq. 3), $p_0=0.6015$:

| $n$ | $m^*$ for 90% from $B_\infty$ | $m^*$ for 90% from $B_q$ | ratio |
|---|---|---|---|
| 128 | **15** | **534** | 36× |
| 512 | **17** | **706** | 42× |

⇒ **$B_q$'s $m^*$ numbers are ~1.5 orders of magnitude too large and should not be used as design
targets.** Note also how weakly both scale in $n$ (15→17 for a 4× larger tree) — §H.1(d) again.

### How to actually use the two bounds

The paper's own guidance is right, for the right reason:
- *"We include $B_\infty$ in our results since it is a tight bound in the case of $q=0$, whereas
  $B_q$ is not. Therefore, estimates derived from $B_\infty$ may be useful in high diversity
  settings."* ⇒ **$B_\infty$ is the usable one.**
- $B_q$ is a **robustness check**: it confirms the *direction* (homoplasy hurts, and hurts more at
  low $\lambda$ / high $q$) without giving usable *magnitudes*.
- ⇒ For our design work: **quote $m^*$ from $B_\infty$, use $B_q$ only qualitatively**, and if we
  need real numbers at moderate $q$, compute them ourselves (see below).

### ⚑ Two concrete improvements, both cheap

1. **Keep $\mathrm{Cov}(X_{\text{in}},X_{\text{out}})$ and use exact convolution.** Removes the
   dominant error; would make $B_q\to B_\infty$ as $q\to0$ up to the prefactor. This is the fix
   that matters.
2. **Shrink the prefactor.** Requiring *all* $\binom n3$ triplets is sufficient but nowhere near
   necessary — a classical result (Semple & Steel; Steel 1992) is that a rooted binary tree on $n$
   leaves is determined by as few as **$n-2$** well-chosen rooted triples. ⚠ Caveat: the guarantee
   would then only apply to a reconstruction procedure restricted to that defining set (BUILD on
   those triples), not to UPGMA. Still — **the $\binom n3$ is an artifact of the proof strategy,
   not a fact about the problem.** Pleasing symmetry if it works: $n-2$ branches ↔ $n-2$ triples.

⇒ Combined with the **A5 frailty** correction (H.6.8), (1)+(2) constitute a well-defined,
self-contained methods contribution: *a corrected reconstruction bound for sequential recorders
that is tight at $q\to0$ and honest about tape correlation.*

### Running tally of every approximation in §3.1, and its direction

| # | Approximation | Direction | Magnitude |
|---|---|---|---|
| 1 | triplet reduction (sufficient, not necessary) | conservative | large (real methods beat triplets) |
| 2 | Boole over triplets/branches | conservative | $\binom n3$ vs $n-2$ = $n^2/6$ |
| 3 | $d_u\le d$, $d_v-d_u\ge\ell$ (Eq. 3) | conservative | moderate |
| 4 | $\mathrm{Cov}(X_{\text{in}},X_{\text{out}})=0$ (Eq. 9) | conservative | **dominant** — grows in $m$ |
| 5 | normal for the tail (Eq. 9) | conservative | small on top of #4 |
| 6 | $\alpha=1$ continuity correction | conservative to $m\!\approx\!20$, then mildly anti | small |
| 7 | **tape independence (A5)** | ⚠ **ANTI-conservative** | unknown, potentially large |

**Six conservative, one anti-conservative — and #7 is the one nobody has quantified.**

## H.6.11 §3.2 — Validation & optimal rates

### What §3.2 establishes (four claims)

1. **Accuracy is single-peaked in $\lambda$.** Too slow ⇒ no information; too fast ⇒ tapes
   saturate before the deep splits. Optimum at intermediate $\lambda$. (This is the $p_0$ U-shape
   of H.6.5 showing up in real reconstructions.)
2. **The width of the usable window grows sharply with $k$, and also with $m$.** Not just the peak
   height — the *range of workable rates*. Practically: more sites per tape buys **tolerance to
   getting $\lambda$ wrong**, which matters because $\lambda$ is not precisely controllable.
3. **Both bounds reproduce the shape**, and $B_q$ is never violated in any regime they test;
   $B_\infty$ is violated (overestimates accuracy) at low $\lambda$ when $q$ is high.
4. **Both bounds skew toward *lower* editing rates than the simulations.** Their diagnosis:
   maximising over triplet depth (Eq. 10) **overestimates the saturation effect at high $\lambda$**.
   More pronounced for $B_q$ than $B_\infty$.

### ⚑ Reproduced independently (`fig32.html` / `sim32.py`)

Balanced tree, $n{=}64$, synchronous ($n_g{=}6$, $\ell{=}1/7$, $d{=}6/7$), $m{=}10$, $q{=}1/16$,
100 replicates/point, UPGMA, scored as exact-topology (RF $=0$):

| $k$ | peak accuracy | $\lambda$ with sim $\ge90\%$ | $\lambda$ with $B_\infty\ge90\%$ |
|---|---|---|---|
| 5 | **0.58** (never reaches 90%) | none | none |
| 9 | 0.99 | 6, 8, 10 | 5, 6, 8 |
| 13 | 1.00 | 6, 8, 10, 13, 16 | 5, 6, 8, 10, 13 |

✓ Claims 1, 2 and 4 all reproduce cleanly. **The $B_\infty$ window is shifted left by exactly one
grid step in both $k$** — the paper's "skew towards lower editing rates," confirmed quantitatively.
Also note **$k=5$ cannot resolve a 64-tip tree with 10 tapes at *any* rate** — a hard structural
limit, not a tuning problem. (Typewriter as published is $k=5$; cf. §H.1 "not yet in a regime
where exact reconstruction is possible for ~1000 cells.")

### Fig. 3 / Fig. 4 details worth keeping

- **$q=1/64$ (high diversity):** $B_\infty$ only *slightly* overestimates accuracy at low rates.
- **$q=1/4$ (low diversity):** $B_\infty$ **significantly** overestimates. ⇒ the $q\ge1/j$ floor
  (mouse $j{=}8\Rightarrow q\ge0.125$) puts real data uncomfortably close to the bad regime.
- **$B_q$ never violated, but not tight — *"and this is a problem which worsens with respect to
  the tree size."*** ⇒ they know about the $\binom n3$ growth (H.6.10).
- **Sub-sampled trees (Fig. 4): both bounds get looser.** They deliberately keep $\ell$ from the
  *full* process even though sampling stretches branches, arguing researchers know the minimum
  *cell-division time*, not the minimum branch length of a subsample. Conservative and honest —
  but it means Fig. 4's looseness is partly self-inflicted. **For our designs, $\ell$ should be
  the resolution we demand, not a biological constant** (their own §3.3 framing).

### §3.2.1 — scoring, and the concession that matters

- **Triplet score** = proportion of simulated distance matrices where *all* triplets satisfy
  Eq. (2). Corresponds exactly to the theory.
- **RF $=0$** = proportion where the reconstructed topology is exactly right.
- *"Theoretically, a reconstruction method based on triplets … would provide the same results as
  the triplet score. In practice, we expect any appropriate tree-reconstruction method to perform
  better."* Fig. 5: **UPGMA and NJ are both more optimistic than the triplet score**, and they
  expect *"a likelihood-based approach (such as Seidel et al.) to perform even better."*

> ⚑ **The concession:** the theory is calibrated to a method nobody uses. Bounds are conservative
> by an amount that (a) is not quantified and (b) varies with parameters. **So the bounds order
> designs correctly but do not measure them.** Use for ranking/screening, not for absolute $m^*$.

## H.6.12 ⚑ Methodological: bounds vs simulation for experimental design

**Question: for real design work, why not just simulate?** Largely right — but the division of
labour is not "simulation is more accurate," it is *what kind of question you are asking*.

**Where simulation clearly wins**
- **Accuracy.** $B_q$'s $m^*$ is ~40× off (H.6.10); $B_\infty$'s is conservative by an unquantified
  amount (§3.2.1).
- **Model richness — decisive for us.** Every row of the A-table that the bounds structurally
  cannot express is trivial to simulate: dropout (**A6**), per-cell frailty (**A5**),
  state-dependent rates (**A2**), time-varying $f$ (**A3**), heterogeneous $k_z,\lambda_z$.
- **It scores the method you will actually run** (SciPhy / LAML), not a triplet oracle.
- The paper's own §3.2 *is* a simulation study.

**Where the theory genuinely earns its keep (four things simulation can't do)**
1. **Inversion, not evaluation.** §3.3 solves *for* $m^*$ in closed form. Simulation is a noisy
   *forward* map; inverting it means grid-searching a 5–6-dimensional space
   $(k,m,\lambda,\ell,n,q)$, with a noisy inner optimisation for $\lambda^*$ at every point.
   $\lambda^*=\arg\max p_0$ is instant analytically.
2. **Scaling laws readable by inspection.** "$\ell$ and $m$ exponential, $n$ logarithmic" falls
   straight out of $1-(n-2)p_0^m$. A simulation grid would let you *fit* that; the formula lets
   you *know* it.
3. **Structural results — the real value here.** ⚑ *"Accuracy depends on $\xi$ only through
   $q$"* is a **theorem** (H.6.7: clock/jump independence ⇒ symbols reach topology solely via
   collisions). Simulation would show two $\xi$ vectors giving similar accuracy; it would never
   tell you the reduction is *exact and universal*. Ditto the $\lambda\ell\gg q/(1-q)$ feasibility
   inequality, and the existence (not location) of the U-shape.
4. **Extrapolation.** Simulation interpolates inside the grid you ran. Park's $n\approx10^5$ needs
   an $O(n^2)$ distance matrix ($10^{10}$ entries) per replicate; the formula is free.

**The trap to avoid:** *a simulator is also a model.* If it omits A5/A6 its answer is exactly as
wrong as the bound's — just more confidently wrong, because it comes with error bars that describe
Monte Carlo noise rather than model error. The algebra at least announces its assumptions.

> **Working rule for our design work:**
> **Theory for structure and screening** (which knob matters, roughly where the $\lambda$ window
> is, order-of-magnitude $m$, and the $q$ / $\ell$ trade-offs) → **simulation for the number**,
> with our actual inference method and our actual nuisance processes.
> **Do not quote $m^*$ from $B_q$.** Quote it from $B_\infty$ as a screen, then simulate.

**What this paper is actually worth to us — none of it is a number:**
(a) $q$ as the sole channel from $\xi$ to topology ⇒ **strong signal degrades the phylogeny**
(§H.1c); (b) the $\lambda\ell$ vs $q/(1-q)$ feasibility inequality; (c) §3.4's multi-rate result
reframing **A1** as a design *feature*; (d) sensitivity to $\ell$ not $n$, reframing the scale
question. All four are structural, and none would have come out of a simulation sweep.

*(Secondary but real: a closed-form design calculation is far easier to defend in a grant or
methods section than "I ran a simulation." Reviewers can check algebra.)*

## H.6.13 ⚑⚑ §3.4 — Multiple editing rates (the highest-value section)

### The problem in one sentence

**One clock cannot time two timescales.** A single rate $\lambda$ resolves a branch of length
$\ell$ only if $\lambda\ell\gtrsim1$; it avoids saturation only if $\lambda\cdot 1\lesssim k$.
Both at once requires $\ell\gtrsim 1/k$. Early development has short cell cycles ⇒ short deep
branches ⇒ exactly the case where no single $\lambda$ works.

### ⚑ The dynamic-range framing (mine, not the paper's)

> **A single tape is an exposure setting.** With $k$ sites at rate $\lambda$ it can time events
> over a window of roughly $[1/\lambda,\;k/\lambda]$ — under-exposed (no edits) below,
> over-exposed (saturated) above. **The dynamic range of one tape is therefore $\approx k$-fold.**
> If the tree's branch lengths span more than $k$-fold, **one rate provably cannot cover it** and
> you must bracket. Typewriter has $k=5$ ⇒ **a 5-fold dynamic range.** That is tiny.

Confirmed numerically below: with $k=5$ the two-rate benefit switches on right around
$\alpha\approx0.2$–$0.35$, i.e. a 3–5× branch-length spread. ✓ The heuristic predicts the
transition.

### ⚠ What "multiple editing rates" actually means

**Two populations of tapes with different constant rates in the same cell** — $\lambda_1$ on
$m_1$ tapes, $\lambda_2$ on $m_2$, with $m_1+m_2\le m$. **Not** a rate that changes over time.
(Correction applied to §H.1(a).) The ensemble behaves like a watch with an hour hand and a second
hand: fast tapes saturate early but resolve the early burst; slow tapes are near-blank early but
still have capacity when the late divisions happen.

⇒ **Experimentally implementable**: different promoters, spacer-mismatch levels, or pegRNA
strengths give different per-tape rates on the same construct. Cheap.

⇒ ⚑ **And SciPhy is already the right model for it.** SciPhy has a **per-tape clock rate $r_z$**
(§0.4 Level 2). We logged that as machinery for *cis* nuisance variation — but it is precisely the
parameter that makes a multi-rate recorder **inferable**. A method assuming one shared rate could
not exploit this design at all. **Design and model already match; nobody has said so.**

### Their procedure (§3.4) — greedy, and explicitly not guaranteed

Given $k,m$ fixed and $\ell_1,\ell_2,n_1,n_2$ known ($\ell_1<\ell_2$: fast early, slow late):
1. If all $m$ tapes at one rate tuned to $\ell_1$ already hit accuracy $c$ — stop.
2. $\lambda_1=\arg\max_\lambda p_0(k,\ell_1,\lambda,d_1)$, $d_1=(n_1-1)\ell_1$ *(the arg**max** is
   a typo for arg**min** — you minimise the failure probability)*.
3. $m_1=\max(m,\lceil m^*\rceil)$ with $m^*$ from the §3.3 formula *(also a typo — must be
   $\min$, or $m_1$ could exceed the budget)*.
4. $m_2=m-m_1$; $\lambda_2=\arg\min_\lambda p_0(k,\ell_2,\lambda,1-\ell_2)$.
Yields $\lambda_1>\lambda_2$, prioritising early divisions. They state plainly: *"this procedure
does not necessarily have a guarantee on the accuracy of the entire tree reconstruction."*

### Reproduction (`sim34.py`) — $n=256$, $k=5$, $m=20$, $q=1/64$, 8 generations, shift after 4, 30 sims

$\ell_1=\alpha\ell_2$; RF = mean Robinson–Foulds to truth (max possible 508):

| $\alpha$ | $\ell_1$ | $\ell_2$ | $\lambda_1$ | $\lambda_2$ | **RF one rate** | **RF two rates** | $p_0(\text{late})^m$ under one rate |
|---|---|---|---|---|---|---|---|
| 0.10 | 0.022 | 0.222 | 28.6 | 3.9 | **186.4** | **37.1** | **1.00** |
| 0.20 | 0.040 | 0.200 | 15.9 | 3.9 | 58.1 | 20.6 | 0.92 |
| 0.35 | 0.061 | 0.174 | 10.5 | 3.8 | 3.4 | 2.1 | 0.31 |
| 0.50 | 0.077 | 0.154 | 8.3 | 3.8 | 1.1 | 0.3 | 0.068 |
| 1.00 | 0.111 | 0.111 | 5.7 | 3.8 | 0.7 | 0.5 | 0.011 |

✓ Reproduces the claim: benefit grows sharply as $\alpha$ falls, negligible at $\alpha=1$.
**Mechanism confirmed by the last column:** with one rate tuned to $\ell_1$, the *late* splits fail
with probability **1.00** at $\alpha=0.1$ — the tapes are fully saturated before the late divisions.
Exactly the paper's *"they quickly become saturated and are unable to resolve the later divisions."*

### ⚑⚑ NEW FINDING — their allocation rule is badly suboptimal

Scanning $m_1$ directly at fixed $\lambda_1,\lambda_2$ (30 sims each):

| $m_1$ fast / $m_2$ slow | 20/0 | **18/2** | 16/4 | 14/6 | 12/8 | **10/10** | 8/12 | 4/16 |
|---|---|---|---|---|---|---|---|---|
| RF, $\alpha=0.10$ | 186.4 | **37.1** | 9.2 | 2.3 | 0.9 | **0.1** | 0.3 | 1.6 |
| RF, $\alpha=0.20$ | 59.0 | **19.2** | 6.3 | 1.7 | 0.4 | 0.4 | **0.1** | 0.5 |

**Their rule picks $m_1=18$ (RF 37.1); a balanced 10/10 split gives RF 0.1 — ~370× better.**

**Why it fails:** step 3 allocates to the early phase *first*, using $m^*$ derived from
$B_\infty$ — which is conservative (§3.2.1), so $m^*$ is inflated (18 of 20 tapes) and the late
phase is starved. **A greedy allocation driven by a conservative bound over-allocates to whichever
phase it processes first.** Fix: optimise $(m_1,\lambda_1,\lambda_2)$ jointly against total
expected failure, or just scan $m_1$ — it is one dimension and costs nothing.

⇒ Concrete, checkable contribution: **the multi-rate design is much stronger than the paper's own
numbers show**, because their allocation heuristic is throwing most of the benefit away.

### Takeaways

1. **Rate heterogeneity across tapes is a cheap design lever that buys what $k$ cannot.** $k$ is
   hard to increase (they cite Liao et al.); rate diversity is easy.
2. **Budget the design by dynamic range:** estimate the spread of cell-cycle lengths across the
   experiment; if it exceeds $\sim k$-fold, a single rate cannot work, full stop.
3. **A1 splits into two distinct things** — temporal/uncontrolled (model it) vs across-tape/
   engineered (exploit it). SciPhy's per-tape $r_z$ already supports the second.
4. **For ENGRAM (§I.3–I.4):** we need a constitutive reference channel anyway. **Make it
   multi-rate at no extra cost.** Pairs perfectly with the $q$ problem: the signal channel has
   peaked $\xi$ (high $q$, poor topology), so the constitutive channels must carry the lineage —
   and they should span the timescales, not just be numerous.
5. ⚠ **Circularity caveat:** the whole procedure assumes $\ell_1,\ell_2,n_1,n_2$ are **known in
   advance**. For a system where the growth phases are what we want to *discover* (HPCS kinetics,
   metastatic seeding), that is exactly what we don't have. **A rate-diverse design is robust
   precisely because it doesn't require knowing them** — a stronger argument for rate diversity
   than the paper's own optimisation framing.

### Next
§3.5 — insertion probabilities: the $\xi$-only-through-$q$ result validated (Fig. 8: skewed over 21
characters vs uniform over 12, both $q=1/12$, indistinguishable) and the estimator
$\hat q=\sum_{i,j}c_{ij}/\sum_{i,j}n_{ij}$ (Fig. 9). Then §4 Discussion (their stated limitations
and future work — the open invitation logged in §H.1).

## H.6.14 The corrected bound — consolidated spec

Four changes, in decreasing order of effect. Each is independently checkable.
1. **Keep $\mathrm{Cov}(X_\text{in},X_\text{out})$** (H.6.8, H.6.10). Dominant error; at $q\to0$
   it is the difference between $10^{-7}$ and $10^{-1}$.
2. **Exact $m$-fold convolution instead of Eq. (9)** (0.13 ms). Removes the CLT-for-a-tail error,
   the missing-$m$ ambiguity, the $\sigma$-vs-$\sigma^2$ ambiguity and the $\alpha$ fudge.
3. **Add a per-cell frailty term on the tapes** (row **A5**) — the only *anti*-conservative
   assumption in the whole paper (H.6.8 ledger).
4. **Optimise $(m_1,\lambda_1,\lambda_2)$ jointly** in the multi-rate design instead of greedily
   (H.6.13; ~370× in my test).
Optional 5: shrink the $\binom n3$ prefactor toward $n-2$ via a defining triple set.

---

# Session I.6 — ⚑⚑ Implementing rate heterogeneity in Typewriter/ENGRAM

*Prompted by §H.6.13. The design question: Mulberry's multi-rate result needs tapes with different
$\lambda$ in the same cell. Can Typewriter do that, and does it force pegRNA↔tape pairing?*

## I.6.1 The governing constraint (restating §0.4's identifiability boundary)

- **Symbol is *trans*.** All Typewriter pegRNAs share **one 20-nt spacer** and differ only in RTT.
  Any pegRNA writes to any tape in the cell.
- **Rate is *cis*.** pegRNA/PE abundance is cell-wide, so it sets a *global* rate and **cannot**
  differentiate tapes within a cell. PE dosage is likewise trans — **tuning PE does not work.**

⇒ **Per-tape rate heterogeneity must come from a *cis* change.** And in this architecture almost
every cis change that affects rate sits *inside the protospacer*, which forces a matched pegRNA.

## I.6.2 ⚑ Why attenuation-by-mismatch is architecturally blocked

Protospacer $=\underbrace{\text{key}}_{1..3}+\underbrace{\text{monomer}_n}_{4..17}+\underbrace{\text{monomer}_{n+1}[1..3]}_{18..20}$, PAM $=\text{monomer}_{n+1}[4..6]$.
**Each monomer base is multi-purposed:**

| Monomer position | Roles it plays |
|---|---|
| 1–3 (`TGA`) | protospacer 18–20 (**PAM-proximal seed**) as $m_{n+1}$; protospacer 4–6 as $m_n$ |
| 4–6 (`TGG`) | **the PAM** as $m_{n+1}$; protospacer 7–9 as $m_n$ |
| 7–14 | protospacer 10–17 as $m_n$ — **inside the seed** |

⇒ **Every monomer position is either double-purposed or deep in the seed.** There is no
"spare" position to place a graded PAM-distal mismatch. The only genuinely PAM-distal handle is
the **key**, and the key is *rewritten by every edit* (insert $=[\text{symbol}][\text{key}]$), so
changing it changes the pegRNA's RTT — i.e. a different pegRNA anyway.

> **This is why Choi et al. had to *screen* 48 tapes rather than design an attenuation series.**
> The constraints are tangled; efficiency is not a tunable dial within one family.

## I.6.3 The options, and whether each forces pegRNA↔tape pairing

| Route | Mechanism | Pairing required? | Verdict |
|---|---|---|---|
| **PE dosage / PE variants** | trans, cell-wide | n/a | ✗ **Does not work** — scales all tapes equally (this is row **A8**) |
| **Mismatch attenuation, shared spacer** | graded R-loop / PBS capture | no | ✗ **Blocked** by I.6.2 |
| **Orthogonal TAPE families** (48-TAPE menu; TAPE-27 >50% more efficient than TAPE-1) | different basal spacer/key/monomer ⇒ different intrinsic efficiency | **YES** — one pegRNA spacer family per tape family | ✓ Works; menu already measured. Rates are *fixed by the menu*, not tunable |
| **Chromatin / locus control** (Bxb1 landing pads at loci of known accessibility) | cis, sequence-neutral | **NO** | ✓ Works, no pairing. But rates are empirical, coarse, and cell-type-dependent (Li et al. 2024) |
| ⚑ **Orthogonal prime editors** (SpCas9-PE + FnCas9-PE / engineered orthogonal PE) | separate editor per family; **rate set by each editor's expression level** | **YES** (by construction) | ✓✓ **Best.** Converts a hard *cis* problem into an easy *trans* one: $\lambda_1,\lambda_2$ continuously tunable by promoter strength |

**Answer to "does it require pegRNAs that only target certain tapes?" — effectively yes**, for any
route with *tunable* rates. The exception (landing-pad chromatin control) needs no pairing but
gives you whatever rates the loci happen to have.

> ⚑ **The orthogonal-PE route is the one to pursue.** Rate becomes a promoter-strength knob rather
> than a sequence-screening exercise, and $\lambda_1/\lambda_2$ is continuously tunable — which is
> exactly what §3.4's joint optimisation over $(m_1,\lambda_1,\lambda_2)$ assumes. Prime editors on
> non-Sp scaffolds exist (FnCas9-PE, Genome Biol 2022) as do engineered orthogonal PE systems
> (Nat Commun 2024). **Verify current performance before committing** — orthologue PEs have
> historically been less efficient.

## I.6.4 ⚑⚑ The same architectural move solves the $q$ problem — and ENGRAM does *not* do it

**Confirmed from the ENGRAM protocol paper:** ENGRAM records *"many different events as unique
insertions to a **shared DNA Tape**"* — 4-bp insertions, up to 256 CRE channels, **all writing to
one tape type.**

⇒ **In the published architecture, a strongly induced channel dominates $\xi$ on *every* tape**, so
$q$ rises everywhere and the lineage degrades everywhere. This is §H.1(c) with no escape hatch.

**And the alphabet does not save you.** With signal share $p$ and the rest uniform over $j$:

| signal share $p$ | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.70 |
|---|---|---|---|---|---|---|---|
| $q$ ($j=256$) | 0.006 | 0.013 | 0.043 | 0.092 | 0.161 | 0.251 | 0.490 |
| $q/(1-q)$ = spurious sites | 0.006 | 0.013 | 0.044 | 0.101 | 0.193 | 0.335 | 0.962 |

At $p=0.5$, $q\approx0.25$ **whether $j$ is 16 or 4096** (0.267 / 0.250). **$q$ is governed by the
skew, not the alphabet size** — so "add more symbols" is not a fix. (Direct corollary of
$q=\sum\xi_i^2\ge\max_i\xi_i^2$.)

**Design rule, from $\lambda\ell\gg q/(1-q)$ (H.6.7).** With $\ell=0.1$ and $\lambda=5$ edits per
experiment, $\lambda\ell=0.5$:
> **On a shared tape, the induced signal symbol must stay below ~30% of all insertions**
> ($q<0.09$, 5× margin) — **40–45% is the hard ceiling.** Above that the lineage signal is
> swamped by homoplasy.

**⇒ The fix is the same orthogonality that buys rate heterogeneity.** Put the signal pegRNAs on
spacer family **S** and the constitutive/lineage pegRNAs on family **L**. Then family L's $\xi$
stays uniform and low-$q$ *no matter how hard the signal fires*, while family S is free to be
maximally skewed — which is exactly what you want, because a peaked $\xi$ is *informative* about
signal and only harmful to topology.

> ⚑⚑ **One architectural decision buys three things at once:**
> **(1)** immunity of the lineage channel to the $q$ problem; **(2)** independently tunable
> $\lambda_L$ vs $\lambda_S$ (Mulberry §3.4 multi-rate); **(3)** a clean constitutive reference
> channel for normalising out **A2/A8** rate confounding (§I.3 risk 1, the HPCS quiescence
> problem). **This is the single most important design conclusion from the Mulberry read, and it
> is a change to ENGRAM's published architecture, not an application of it.**

**⚠ Cost to check in the ENGRAM papers:** two spacer families means two pegRNA scaffolds and two
tape arrays per cell, splitting the integration budget; and the signal channel loses the passive
lineage information it currently gets for free. Whether the $m$ split is affordable is the first
thing to work out from Chen et al.

## I.6.5 ⚠ THREE REVISIONS after reading Chen et al. in full (see §I.7)

1. **Shared tape is their central design *claim*, not an oversight.** *"1,024–4,096 unique
   biological signals could theoretically be recorded within the same cell, all competing to write
   to a shared DNA Tape. The advantage of a shared recording medium … is particularly manifest in
   the combination of ENGRAM and DNA Typewriter."* My orthogonal-tape proposal cuts against their
   thesis — keep it, but frame it as a cost-benefit trade, not a fix to an error.
2. **The two-channel design is already Shendure's stated roadmap** (§I.7.5). Not novel. What *is*
   missing is the quantitative constraint on it.
3. ⚠ **I was wrong that PE dosage is useless.** It is useless for differentiating tapes *within a
   cell* (still true, row A8). But **Dox-inducible PEmax gives a controllable, time-varying
   *global* rate** — and Chen et al. *use it exactly that way* to gate 24-h recording windows in
   gastruloids. That is a demonstrated third route to Mulberry's multi-rate result, and the
   easiest one. Revised route table:

| Route to multi-rate | Requires pegRNA↔tape pairing? | Status |
|---|---|---|
| **Inducible PE, rate varying in *time*** | **NO** | ✓✓ **Demonstrated in Chen et al.** Easiest. Needs exogenous control (Dox) + knowing when the phases are |
| Orthogonal PE/tape families, rate varying *across tapes* | YES | Hardest; but works in vivo with no exogenous control, and no phase knowledge needed |
| Chromatin/landing-pad locus control | NO | Coarse, uncontrolled |
| PE dosage across tapes | — | ✗ Impossible (trans) |

## I.6.6 Questions we set for ENGRAM — all now answered in §I.7

1. **Are channels spacer-orthogonal, or all on one tape?** (Protocol says shared — confirm, and
   check whether orthogonal tapes were tried.)
2. **How peaked does $\xi$ get at full induction?** ⇒ $q$ directly, ⇒ feasibility via I.6.4.
3. **Transfer function** promoter activity → pegRNA abundance → insertion frequency. This *is* the
   $\xi_i(\text{signal})$ map the whole extension needs. Is it linear? Saturating?
4. **Temporal integration window** — does a 2-hour burst register, or only sustained signal? Sets
   the achievable $\ell$, hence everything in H.6.
5. **Is there already a constitutive reference channel?** (Tet-On, Wnt, NF-κB were the three demo
   channels — was a fourth always-on channel included?)
6. **Crosstalk** — the tweet claims "minimal crosstalk"; quantify it, because crosstalk is exactly
   what would break the family-L isolation argument above.

---

# Session I.7 — ⚑⚑ Chen et al. 2024, ENGRAM (Nature 632:1073–1081) — close read

*Read against the six questions set in §I.6.5. All six answered. Three findings change the project.*

## I.7.1 Architecture — confirmed, and it is deliberate

- **Csy4-based Pol2 pegRNA release.** CRE→minP drives a Pol2 transcript containing
  `csy4–pegRNA–csy4`; Csy4 (Cas6f) excises the pegRNA. This is the trick that lets an *enhancer*
  drive a guide RNA. Three designs; **5′ ENGRAM** used throughout (3′ FT had the best S/N —
  13.9% vs 0.58% editing ±TNF, 23.8-fold — but 5′ is easier to clone).
- ⚑ **Shared spacer, shared tape, by design:** *"many ENGRAM recorders … may **share a common
  spacer** while encoding different insertions, such that signal-specific symbols will be written
  to a **shared location**."* And in the Discussion: *"1,024–4,096 unique biological signals …
  all competing to write to a shared DNA Tape. The advantage of a shared recording medium … is
  particularly manifest in the combination of ENGRAM and DNA Typewriter."*
  ⇒ **Confirms §I.6 exactly. Shared tape is the thesis.**
- **DNA Tape**: endogenous HEK3 (2–3 copies) or synthetic HEK3 via piggyBac (**~20 copies/cell**).
  ⇒ **$m\approx20$**, right in Mulberry's simulated range (10/30/50/70).
- **DTT** in the Typewriter-coupled experiment: **five units ⇒ $k=5$.**
- Coupling is trivial by construction: *"combining ENGRAM and DNA Typewriter requires only that
  **some or all** of the symbols of a DNA Typewriter be driven by ENGRAM recorders."*

## I.7.2 ⚑⚑ FINDING 1 — ENGRAM's editing rate is ~7× below Mulberry-optimal

Measured in gastruloids: **8–15% of DNA Tapes edited per 24-h PEmax window**, *"trending downwards
with time."* Over a 5-day course that is $\lambda_{\text{eff}}\approx0.4$–$0.75$ edits per tape per
experiment. (Corroborating: only 2.8% of endogenous HEK3 edited at 5 days in K562.)

| Setting | Mulberry $\lambda^*$ | ENGRAM $\lambda\approx0.58$ | gap | per-split failure, $m=20$ |
|---|---|---|---|---|
| $k{=}5$, resolve 24-h windows ($\ell{=}0.2$) | **3.85** | 0.58 | **6.6×** | $1.3\times10^{-5}$ → **0.098** |
| $k{=}5$, resolve 12-h windows ($\ell{=}0.1$) | 3.75 | 0.58 | 6.5× | $4.8\times10^{-3}$ → **0.314** |
| $k{=}6$, 24-h windows | 4.55 | 0.58 | 7.8× | $9.2\times10^{-7}$ → 0.098 |

> **ENGRAM as published sits deep in Mulberry's "λ too low" failure regime — the left limb of the
> Fig. 2 curve.** Per-split failure degrades by ~7,600×. **This is the single most actionable
> number from the whole reading programme**, and it is *good news*: the fix is more editing
> (PEmax over PE2 already helped; longer induction; better pegRNAs), not a redesign. It also means
> **any simulation we run should sweep λ upward from the published operating point**, not around it.

⚠ *"Trending downwards with time"* = rows **A1/A2** confirmed inside the ENGRAM system itself,
not just in the mouse embryo.

## I.7.3 ⚑⚑ FINDING 2 — the published proof-of-concept runs at $q=0.5$

The ENGRAM+Typewriter experiment used exactly **two symbols** (Tet-On, WNT) on a 5-unit DTT.
Uniform over $j=2$ ⇒ $\boxed{q=0.5}$ — **double Mulberry's "low diversity" regime ($q=1/4$), where
$B_\infty$ already fails badly.** The 26 PCA variables were "two symbols × five positions + four
bigrams × four position-pairs", confirming $j=2$.

⇒ **Direct confirmation of §H.1(c)/§I.6.4:** a signal-only tape cannot carry lineage. Not a
hypothetical — it is the actual parameter regime of the published demonstration. They weren't
doing lineage there, so no error; but it means **nobody has yet run ENGRAM+Typewriter at a $q$
where the lineage is recoverable.**

## I.7.4 The $\xi_i(\text{signal})$ transfer function — measured, and it factorises

$$\xi_i(\text{signal})\;\propto\;\underbrace{\beta_i}_{\text{barcode-intrinsic bias}}\times\underbrace{\sigma(\text{signal})}_{\text{sigmoid dose–response}}$$

- **Dose–response is sigmoid, not linear.** Fitted by nonlinear regression; EC50s: Dox
  **0.17 µg/ml**, TNF **2.5 ng/ml**, CHIR **2.2 µM**. WNT is *"almost switch-like across a
  fourfold range"* of CHIR.
- **Recording is a function of both intensity and duration** (Fig. 3f,g) — i.e. it integrates.
- **Relative activity is faithful**: 15.1-fold ENGRAM vs 15-fold MPRA for a high/low CRE pair;
  strong MPRA-vs-ENGRAM correlation over 300 CREs.
- ⚑ **Barcode-intrinsic bias is large but predictable.** 948 degenerate 5-mers; editing scores span
  several log2 units; **minimum free energy of the pegRNA secondary structure is the top predictive
  feature**; lasso model reaches Pearson **0.907** on held-out barcodes.
  ⇒ **$\beta_i$ is measurable *and* designable.** We can pick barcodes to flatten $\beta$, which
  directly lowers $q$ at no biological cost. **Free win.**

## I.7.5 ⚑⚑⚑ FINDING 3 — the two-channel design is Shendure's stated roadmap, and its constraint is unstated

Discussion, verbatim:
> *"we envision that hundreds to thousands of biological signals could be coupled to the ordered
> writing of signal-specific insertions ('symbols') to DTT(s) … **A further set of non-specific
> symbols, stochastically written to the same Tape(s), would facilitate the capture of cell
> lineage.**"*

⇒ **Exactly the two-channel design (§I.3–I.4, §I.6.4) — same tape, symbol-partitioned.** So it is
*their* target, not ours. **What is entirely absent is any quantitative constraint on the mix.**
That is precisely what Mulberry supplies and what §I.6.4 derives:

> **The signal symbols' share of insertions must stay below ~30% (hard ceiling 40–45%), or the
> "non-specific symbols" they are relying on for lineage are swamped by homoplasy.**

**And there is a second, simpler consequence they don't state:** because $q$ depends only on the
*share*, and the non-specific set can be made arbitrarily diverse and abundant, **the mix is a
tunable design parameter** — dial constitutive pegRNA abundance up until $q$ is safe. Cost: signal
dilution, hence more tapes. **That trade-off curve is a clean, self-contained contribution**, and
it does *not* require abandoning the shared tape.

⇒ **Revised recommendation (supersedes §I.6.4's orthogonal-tape emphasis):** the shared-tape,
symbol-partitioned design is fine *provided the share is controlled*. Reserve orthogonal tapes for
the case where the signal must be allowed to saturate.

## I.7.6 Other things worth keeping

- **Crosstalk is low.** Three recorders co-deployed, 8 agonist combinations: barcode abundances
  *"highly dependent on the combination of stimuli applied … consistent with the notion that these
  signalling pathways are orthogonal."* 15–20 integrations/cell. ✓ Q6 answered.
- **Temporal resolution: 6 h of stimulation suffices to exceed background.** NF-κB kinetics faster
  than WNT. ⇒ **$\ell_{\min}\approx6$ h**; in a 5-day run that is $\ell\approx0.05$. **A real number
  for our design calculations** — and note §H.1(d): halving $\ell$ costs a lot of tapes.
- **Background editing is real but bounded** and *plateaus* after several days (0.58% for 3′ FT).
  ⚑ With ~300 recorders each carrying its own barcode, background writes a *diverse* symbol set at
  low rate — i.e. **a nascent high-diversity constitutive channel already exists for free.**
  Worth quantifying: is background alone enough to keep $q$ low?
- **Bigrams work, but are not sufficient.** *"Within each programme class we could distinguish not
  only different orders but also different timings, solely based on the ratios of bigrams"* — but
  **not across classes** (serial vs layered vs pulse). Needed all 26 variables + PCA/random forest.
  ⇒ **Confirms §1a.9 *and* bounds it**: the adjacent-position dependence carries the temporal
  signal, but a bigram summary loses information a full likelihood would keep. **This is the
  argument for our model, stated by their own data.**
- **PEmax was Dox-inducible** to gate 24-h windows ⇒ see §I.6.5 revision 3: a demonstrated
  time-varying global rate.
- **Their stated open problems** (two more open invitations):
  1. *"the deconvolution of ENGRAM signals, when coupled to DNA Typewriter, will pose new
     algorithmic challenges."* ⇐ **this is our project, named by them.**
  2. Not suited to fast timescales, or to signals not couplable to CREs.
  3. Thousands of recorders hard by random integration → they propose consolidating to a single
     *"recorder locus."*
  4. Recording endogenous gene expression by integrating minimal cassettes into gene bodies
     **failed** (Supplementary Note 2) — barcode proportions did not track expression. Three
     candidate explanations offered. ⇒ **Don't build on that route.**
- **Readout is bulk gDNA.** They flag single-cell readout via a *"transcriptionally expressed DNA
  Tape"* as the next step, as done for DNA Typewriter. ⇒ Our per-cell latent-state model needs
  that; check whether the mouse-embryo/Park pipelines already have it.

## I.7.7 Where this leaves the project

**The gap is now sharply defined and small:**
- Shendure has stated the destination (signal symbols + non-specific lineage symbols on a shared DTT).
- Chen et al. have built every component and demonstrated order-recovery from bigrams.
- Mulberry has supplied the reconstruction theory but only for a *homogeneous, single-$\xi$* tape.
- **Nobody has connected them.** Specifically, nobody has answered: *given a signal channel with
  share $p$ and a constitutive channel with diversity $j$, at rate $\lambda$, on $m$ tapes of $k$
  sites — can you recover both the tree and the signal history, and at what resolution $\ell$?*

**That question is our thesis, and all three inputs to it now exist.**

Immediate next steps, revised:
1. **Simulate the joint design space** $(p,\lambda,m,k,j,\ell)$ under the ENGRAM-measured
   parameters, scoring tree accuracy *and* signal-history recovery together. The rate finding
   (§I.7.2) says sweep $\lambda$ up 5–10× from published.
2. **Derive the share/dilution trade-off curve** (§I.7.5) — self-contained, publishable, and it is
   the missing quantitative constraint on Shendure's own roadmap.
3. Park et al. compatibility measurement (§D.4) — still cheap, still decisive.
4. The corrected bound (§H.6.14).

---

# Reading brief — Schiffman et al. 2024, PATH (Nat Genet 56:2174–2184)

*Pre-read orientation only. Fill in with a proper close read; this is the map, not the territory.*

## What PATH is, in four lines

1. **Metric.** Adapts **Moran's $I$** (phylogenetic autocorrelation) to somatic evolution ⇒
   heritability vs plasticity of a cell state. Plus **cross-correlation** between distinct states.
2. **⚑ The actual novelty.** They *formally link the correlation metric to a Markov transition
   model*, so correlations become **transition-rate inferences**. Their words: *"phylogenetic
   autocorrelation has not been linked with a model of phenotypic transition dynamics, as we
   accomplish here."*
3. **PATHpro** extends this to state-specific **proliferation** rates (multitype branching process).
4. **Cheap.** Accuracy *"comparable to MLE"* but far faster; scales to many states × many cells.

**Input = a phylogeny + leaf annotations.** It does not build trees. That is the whole reason it
sits downstream of our problem.

## ⇒ Verdict on priority

**Worth reading — but as a survey, not an equation-by-equation session like Mulberry.** It is not
blocking any of the §I.7.7 next steps. **Its priority is coupled to the Park thread (§D.4):
Schiffman is a Park et al. coauthor**, so Park's heritability analysis is very likely PATH. Read it
*immediately before* touching Park's data, not before the design simulations.

## The five things to attend to

1. **⚑⚑ Branch-length imputation.** PATH can *impute* branch lengths from a somatic-evolution model
   + the sampling rate — **barcode-independent**, unlike reconstructed lengths which depend on scar
   count. Reported: accuracy rises with barcode length, **but for *short* barcodes imputation
   outperforms reconstructed branch lengths** (Fig. 2c). **Typewriter is $k=5$–6 — short.**
   ⇒ Directly bears on our architecture (§F.3/§H.2): possibly we do not need good branch lengths
   at all for state-dynamics questions. **Check the Supplementary Note for the imputation model.**
2. **⚑ The correlation → transition-rate derivation** (Methods + Supplementary Note). This is the
   part to actually work through. It is the moment-matching analogue of what we would do with a
   likelihood, and it is the thing to either reuse or position against.
3. **⚑ The lineage-tracing confounder benchmarks** (Fig. 3): accuracy vs *number of cut sites*, vs
   *dropout*, vs *sampling rate*. This is the closest published analogue to Mulberry's design
   analysis but for **state inference rather than topology**. Compare their conclusions with
   §H.6 — do the two agree about what limits you?
4. **⚑⚑ The equilibrium assumption.** *"both approaches assume cell state frequencies are near
   equilibrium."* ⚠ **Our applications are explicitly non-equilibrium** — gastruloid
   differentiation, HPCS emergence, metastatic seeding. **This is the single biggest reason PATH
   cannot simply be adopted for the signal-history project**, and the cleanest statement of what
   our method would add. Note exactly how they phrase and defend it.
5. **What they claim about Park-type data**, so we know what is already asserted before we
   re-analyse it.

## Read past quickly

- The three biology applications (pancreatic EMT, GBM, B-ALL) — read the *conclusions*, skip the
  mechanism. Come back only if the GBM astrocyte-bridge result becomes a template for an
  HPCS-intermediate claim.
- Moran's $I$ background / Pagel's $\lambda$ / Blomberg's $K$ — standard, skim.
- The CNV–phenotype cross-correlation section, unless we end up wanting genotype–phenotype links.

## ⚑ The conceptual question to hold while reading

PATH is a **moment-based** estimator: compress the tree+states into a correlation statistic, then
invert a model to get rates. We have been framing our options as **likelihood vs SBI** (§1a.9,
A5–A9). PATH is a third option we have not costed: **moment matching** — cheap, scalable, robust to
branch-length error, but statistically inefficient and blind to anything the summary discards.

> **And we already have direct evidence of that trade-off in our own domain.** ENGRAM's bigram
> analysis (§I.7.6) recovered order *and* timing within a programme class — but **could not
> distinguish serial vs layered vs pulse patterns**; that needed all 26 variables plus PCA.
> **Bigrams are to signal history what Moran's $I$ is to state dynamics: a fast summary that
> captures a great deal and provably not everything.**

⇒ Read PATH with this question: *what does the summary throw away, and is it the part we need?*
If the answer is "the non-equilibrium transient" — which points 4 and the ENGRAM bigram result both
suggest — then **that is the argument for our method, stated in one sentence.**

## Where it would slot in if we adopt it

Plausible hybrid: **our likelihood for the tape/signal layer** (which PATH cannot do) →
**PATH/PATHpro for downstream state-dynamics questions** on the resulting tree (which is cheap and
already validated). Worth deciding early whether we are competing with PATH or composing with it.
Current read: **composing.**

---

# D.4b — The Park compatibility check, as an executable protocol

*Concretising §D.4's "one-afternoon computation." Includes a reframing: it is a **residual**
analysis, not just a measurement, because we already have a prediction.*

## The object being measured

**Character** = a pair (tape $z$, prefix $p$), $|p|\ge1$. Its **clade** $S_{z,p}$ = the set of cells
whose tape $z$ carries $p$ as a prefix.

**Two characters are compatible** iff their clades are **nested** ($S_1\subseteq S_2$ or
$S_2\subseteq S_1$) or **disjoint** ($S_1\cap S_2=\varnothing$). Partial overlap = incompatible:
no single tree can carry both.

**Only cross-tape pairs need testing.** Within a tape the characters form a trie and are nested by
construction (§D.4 fact 1) — a single tape can never be internally incompatible. *This is the
sequential-recorder-specific fact that makes the whole thing cheap.*

$$\text{compatibility fraction}=\frac{\#\{\text{compatible cross-tape pairs}\}}{\#\{\text{cross-tape pairs tested}\}}$$

## ⚑ Reframing: this is a residual analysis

We **already have the prediction**, from our own notes: Park uses NNNNGGA, $M\approx256$,
$\sum_i \xi_i^2\approx\mathbf{0.4\%}$ — better than Mulberry's "high diversity" regime ($q=1/64$).
Spurious sharing of a length-$L$ prefix costs $\approx q^L$, so **homoplasy alone predicts
incompatibility in the low tenths of a percent.**

> **⇒ The informative quantity is not the fraction. It is the GAP between predicted (~0.4%) and
> measured.** Everything in that gap is dropout, sequencing/PCR error, doublets, tape
> misassignment, copy-number artefacts (§0.2b `tape/copynumber.py`) or genuine model violation —
> i.e. **rows A6/A7/A9 measured directly, on real data, for the first time.**
> That is a more interesting result than "the perfect-phylogeny route is live", and it comes free.

> ### ⚠⚠ CORRECTION (2026-08-28) — the measured $q$ is **4.3× larger** than this prediction
>
> Measured on the delivered edit tables (`analyses/2026-08_park-compatibility`, Step 0):
>
> | table | cells | $M_{\text{obs}}$ | edits | $q=\sum_i \xi_i^2$ | $1/q$ | missing |
> |---|---|---|---|---|---|---|
> | Initial  | 37,810 | 244 | 13.9 M | 0.01691 | 59.2 | 62.2% |
> | Mouse1   | 12,232 | 197 | 5.81 M | 0.01740 | 57.5 | 51.2% |
> | Mouse2   |  6,899 | 184 | 3.34 M | 0.01717 | 58.2 | 50.0% |
> | Mouse3   |  2,904 | 167 | 1.32 M | 0.01726 | 57.9 | 53.2% |
> | Subclone | 39,606 | 218 | 22.1 M | 0.01689 | 59.2 | 42.7% |
>
> **$q\approx0.0170$, not $0.004$.** The 0.4% figure assumed a *flat* distribution over all
> $M=256$ symbols ($1/256=0.39\%$). Neither half holds: only 167–244 of the 256 are ever observed,
> and the realised distribution is far from flat — the **effective alphabet is $1/q\approx58$,
> not 256.** The single commonest symbol `ATATGGA` carries 4.4% on its own.
>
> Two consequences:
>
> 1. **The comparison to Mulberry inverts.** Their "high diversity" regime is $q=1/64=0.0156$.
>    Park sits at $0.0170$ — *slightly worse* than that benchmark, not "better than" it, and not
>    "near homoplasy-free" (the same error appears at §1629's comparison table, line ~1681).
> 2. **The §D.4b null moves up by $\approx4.3\times$ at $L=1$.** Homoplasy alone no longer predicts
>    "low tenths of a percent". Recompute the expected incompatibility at $q=0.0170$ *before*
>    reading any gap as dropout/error. The residual framing survives intact — only the number
>    being subtracted changes.
>
> Consistency across five independently sequenced tables (0.0169–0.0174) makes this a property of
> the pegRNA pool, not of any one sample. Reproduce with
> `analyses/2026-08_park-compatibility/src/01_symbol_composition.py`.

> ### ⚠⚠ CORRECTION 2 (2026-08-28) — "$q^L$" is the wrong mechanism, and it understates homoplasy
>
> The paragraph above says spurious sharing of a length-$L$ prefix costs $\approx q^L$, so that
> homoplasy vanishes for $L\ge2$. **That describes the wrong event.**
>
> $q^L$ is the probability that two lineages which diverged *before any edit on the tape*
> independently fill it with the same $L$ symbols — all $L$ positions must collide. Real, correctly
> computed, and rare.
>
> But a character $(z,p)$ is acquired on the edge where its **last** symbol $p_L$ is written, in a
> lineage that already carries $p_1\dots p_{L-1}$. It gains a second origin when $p_L$ is written
> again, independently, elsewhere within the clade that already shares $p_1\dots p_{L-1}$ — those
> first $L-1$ symbols being shared *by descent*, and costing nothing. **One collision, not $L$.**
>
> ⇒ There is no $q^L$ collapse. Deep prefixes are not automatically clean, and "low tenths of a
> percent" was never the right null.
>
> **The correct scaling.** Let $m$ be the number of independent write events at site $L$ within the
> parent prefix's clade. Each draws from $\xi$, so the expected number of colliding pairs is
>
> $$\binom{m}{2}\,q$$
>
> — linear in $q$, quadratic in $m$. Setting it to 1 gives a **birthday threshold**
>
> $$m^{*}=\sqrt{2/q}\approx\sqrt{2/0.0170}\approx 11$$
>
> so roughly eleven independent writes at one site in one prefix clade already produce an expected
> recurrence. The governing variable is $m$ — i.e. **clade size** — not prefix length.
>
> ⚑ **Measured (`analyses/2026-08_park-compatibility`, scripts 06–09).** $m$ is not observable, but
> it is estimable from the distinct-symbol count $s$ by inverting $\mathbb{E}[s\mid m]=\sum_i(1-(1-\xi_i)^m)$
> against the measured $\xi$, and better still by the Poissonised MLE
> $\log L(m)=\sum_{i\in A}\log(1-e^{-m\xi_i})-m(1-W_A)$. Over **1,567,321 prefix nodes**:
>
> - **96.0% sit below $m^{*}$.** Median $s$ is 1–3 almost everywhere — the single-origin behaviour a
>   perfect phylogeny predicts.
> - Mean recurrences rise with clade size: 0.01 (clade 3) → 0.53 (21–50) → 5.75 (101–500) → 50.5 (>500).
> - **Homoplasy is concentrated: the top 1% of nodes carry 65% of it, the top 10% carry 87%** — and
>   this holds *within a single clone*, not only across the pooled arms.
>
> ⇒ By this section's own step-5 rule that is the **favourable** world: conflict concentrated on a
> few characters is removable by deleting them, rather than spread thin enough to make the
> maximal-compatible-set problem genuinely hard. The two corrections above both *raised* the
> predicted homoplasy, but Park's clades turn out small enough that $m^{*}$ is rarely reached.
> The metastasis mice are near homoplasy-free at their actual clade sizes; the subclone arm
> (median clade 933 at the root of the trie) carries almost all of it, by design.


## Procedure

1. **Work within a clone.** Park has ~75 clones; cells in different clones share no ancestry, so
   cross-clone pairs are trivially disjoint and would inflate the fraction. Report per clone.
2. **Filter characters.** Keep those with clade size $\ge3$ cells (singletons are uninformative and
   noise-dominated). Record how many characters survive — that is itself a data-quality readout.
3. **Bitset the clades.** `np.packbits` per character; nestedness/disjointness are then
   `popcount(A & ~B)` / `popcount(A & B)`. Millions of pairs are seconds.
4. **⚠ Three-valued missing data.** A cell whose tape $z$ is undetermined is neither in nor out of
   $S_{z,p}$. **Compute the fraction both ways** — missing-as-absent vs missing-excluded — and
   report both. The spread between them is the dropout sensitivity, and with ~50% missingness in
   this technology it will not be small. *Do not report a single number.*
5. **Characterise the conflict graph, not just its size.** Degree distribution over characters:
   - conflicts **concentrated** on a few characters ⇒ drop them, route is live and easy;
   - conflicts **spread thin** ⇒ maximal-compatible-set is genuinely hard, route needs real
     machinery.
   §D.4 step 2 is NP-hard in general but easy when conflicts are sparse — **this is the step that
   tells you which world you are in.**

## ⚑⚑ The real output is $C$, not the fraction

The fraction is a means. The payoff metric is:

> $C$ = **how many of the $n-1$ internal nodes the compatible skeleton fixes outright.**

That is what determines the search-space reduction (§D.4: "search only $n-1-C$"). Build the laminar
family from the maximal compatible set, count its internal nodes, report $C/(n-1)$ per clone.
**A high compatibility fraction with a low $C$ would be a null result** — lots of agreeing
characters that all describe the same few clades. Measure $C$ or the check has not answered the
question.

## Decision rule

| measured compatibility | reading |
|---|---|
| $\gtrsim95\%$ **and** $C/(n-1)$ high | route is live; proceed to the skeleton + local likelihood design |
| 80–95% | live but needs conflict resolution; check whether conflicts are concentrated |
| $<80\%$ | something other than homoplasy dominates — needs an error model first, different project |
| any, with a large missing-as-absent vs missing-excluded spread | **dropout is the binding constraint**, which is row **A6** and argues for the SBI route |

## Sequencing with the second experiment

§D.4's follow-on (line ~2092) is the same afternoon's work once the skeleton exists, and is
arguably the bigger prize: **assign each edit to the branch where it first appears, bin by time,
plot symbol composition vs time. Park's pegRNAs are not CRE-driven ⇒ composition SHOULD be flat.**
A clean negative control for the entire $\xi(t)$ estimator, on real data, before ever pointing it at
ENGRAM. **Do these together.**

## ⚠ Step 0, which precedes both

Neither computation can start without a **per-cell × per-tape ordered-genotype matrix** in memory,
with the missingness encoding understood. Getting that, and reproducing one published summary
number from it (e.g. ~254 edits/cell by day 7, or ~150,000 unique insertion patterns), is the
actual first task — and a better first contact with the data than either analysis.
