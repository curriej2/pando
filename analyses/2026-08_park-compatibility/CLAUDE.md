# Analysis: park-compatibility

**Question.** Does the cross-tape character-compatibility structure of the Park data support the
perfect-phylogeny route (§D.4), and how large is the gap between homoplasy-predicted and measured
incompatibility — i.e. how much dropout/error/model violation is really there?

**Protocol.** `notes/sciphy_notes.md` §D.4b, executed in order. §D.4 for the underlying argument.

**Inputs.** `/data1/choij10/justin/pando/data/cancer_metastasis/` (symlink target of `data/`;
gitignored, never copy out). Justin's own copy of the Park, Chang et al. 2026 tables from Jihye Park.
- `{Initial,Mouse1,Mouse2,Mouse3,Subclone}_EditTable_filtered.csv` — cell × (166 tapes × 6 sites),
  wide, 997 columns. Values are NNNNGGA 7-mers; missing is the literal string `None`.
- `clonalbc_percell_hamming1_corrected.csv` — 167,736 rows, `CellID,Sample,ClonalBC_raw,ClonalBC`.
  **This is the clone assignment**, needed for §D.4b Procedure step 1 (work within a clone).

**⭐⭐ STATE (2026-09-14, session 12) — READ THIS FIRST. ⭐ THE EVIDENCE FOR HERITABLE TAPE
SILENCING IS CONSIDERED SUFFICIENT (Justin's call, recorded). The figure programme for it is CLOSED;
the next thread is the SIMULATOR.**

**The four figures that carry it**, all in `figures/`, all five arms:
| fig | script | what it shows |
|---|---|---|
| `fig6a_cost_benefit` | `77`,`80` | the ladder — held-out gain against parameters spent, per rung |
| `fig6b_deviance` | `77`,`80` | $Q$, the lineage share of what the technical model left |
| `fig6c_unanimity_{arm}_k{3,5}` | `83`,`84` | per-tape unanimity, missing vs present, against the diagonal |
| `fig6d_setdiag_{arm}_k5` | `83`,`85` | the companion — capture matched, relatedness not |

**⭐ THE HEADLINE NUMBERS.** Ladder (floor 100, $k$=20, nats/cell obs−null): Subclone **+8.31±0.09**
($Q$=**15.9%**) · Pre-TX **+4.59±0.12** (7.1%) · Mouse2 +2.47 · Mouse3 +2.19 · Mouse1 +1.48,
**positive 5/5 at 16–90 se**, $w$=+0.78…+1.55. Per-tape unanimity ($k$=5, common tapes, median fold
relatives ÷ capture-matched): Pre-TX **9.77× missing vs 1.22× present (ratio 8.0)** · Subclone
**3.84× vs 1.15× (3.3)** · mice 1.16–1.30× vs ~1.03×.

**⚠⚠ THE FOUR CORRECTIONS THIS THREAD PAID FOR — do not reintroduce any of them.**
1. **The null must EXCLUDE SELF.** $M_3$ forces $\sum_{c\in C}r_{cz}=0$, so a self-excluding
   neighbour set carries $-1/(n_C-1)$ of the cell's own residual. `74`'s row-permutation null does
   not, so obs−null cannot cancel it. Fatal at $n_C=31$ (Pre-TX went NEGATIVE), negligible at 3,387.
   ⇒ **null = $k$ random clone-mates, self excluded.** ⚠ `74`'s small-clone numbers remain in doubt.
2. **A cell could be its own nearest relative** (`74`, pre-fix): self and cross-clone pairs both
   $-\infty$, tied under `argsort`. 71.5% of Pre-TX cells at $k$=50. Fixed by blocking per clone.
3. **Stratifying on $\tilde p$ does not hold the truth fixed**, because $\tilde p$ is an ESTIMATE.
   Simulated with NO lineage structure, the "arbitrary clone-mates" curve still rose 14.9%→28.2%.
   ⇒ the model-stratified Fig C overclaimed; the per-tape version needs no model.
4. **Capture matching must be cell-by-cell on RANK, not by decile.** Deciles left Mouse2 at rel 37.8
   vs mat 40.9 and **inflated the folds by up to 20%**. Related sets really are worse captured
   (Mouse2 37.4 vs 54.5 unmatched), so this control is load-bearing, not a formality.

**⚠ Three framing rules learned here.** (i) **Per tape, never pooled** — the pooled unanimity number
was a mixture dominated by near-dead tapes ($0.22^{20}=5\times10^{-14}$ expected vs 4.3% observed).
(ii) **Common-tape convention**: all-$k$-missing needs a high rate, all-$k$-present a low one, so the
two panels otherwise sit on opposite ends of the range. (iii) **$k$=5 maximises the discriminator**;
$k$=2 is the variogram in disguise, $k$=20 leaves Pre-TX only 9 common tapes.

**⇒ NEXT: THE SIMULATOR.** Calibration target per arm, through the identical scripts: the ladder's
4-vector $(\Delta_1..\Delta_4)$ — including the **sign** of $M_2$, which is negative exactly where
the cell filter is strict — the variogram's slope and convexity, the nats-by-$k$ inversion, and the
per-tape unanimity fold. ⚠ These are **joint** constraints on (tree, editing, dropout), not clean
per-parameter ones: the variogram's $x$-axis IS the editing process. Calibrate tree and editing first
against statistics that do not involve dropout, then bring these in. Their value is falsification,
not fitting — tune on the ladder, validate on the variogram.

⚠ **PARKED with a known flaw: the $\beta_z$ skew test** (`78`). Dispersion is hugely inflated
($z$=11–89 in 5/5 arms) but moment skew is not, because centring $\gamma$ on the clone-weighted mean
turns a one-directional effect two-sided when there are few clones (4–11 here). Rebuild against a low
quantile or an out-of-arm baseline before drawing any conclusion about directionality.

---

**STATE (2026-09-11, session 11). THE LADDER IS BUILT (`77`) AND IT FOUND A
BUG IN THE NULL THAT ALSO AFFECTS `74`.** Full spec and results in README "The nested ladder".

**⚠⚠ THE NULL MUST EXCLUDE SELF.** $M_3$ forces $\sum_{c\in C} r_{cz}=0$, so a neighbour set
excluding $c$ carries $-1/(n_C-1)$ of $c$'s OWN residual. `74`'s null permutes residual ROWS across
fixed neighbour SLOTS, giving a uniform subset of ALL $n_C$ cells INCLUDING $c$ — expectation zero,
so it carries NO leak and obs−null cannot cancel it. Verified by simulation (400 reps): observed
−0.0500 at $n_C=21$ vs theory $-1/20$, old null +0.0005, self-excluding null −0.0500. **Negligible at
$n_C=3387$ ($3\times10^{-4}$), fatal at $n_C=31$** — Pre-TX returned a NEGATIVE lineage rung until
fixed. ⇒ **the null is now "$k$ random clone-mates, self excluded"**, which is also the question we
mean: are the NEAREST relatives more informative than ARBITRARY ones?
⚠ **`74`'s small-clone numbers are in question** (Pre-TX pooled +2.86); its large-clone ones are not.

**⭐ THE LADDER, floor 100, $k=20$ (nats/cell, obs−null, 5 splits × 2 masks).** Fitted WITHIN CLONE,
held-out ENTRIES, two-stage so $\Delta_4$ is a lower bound.
Subclone **+8.31±0.09 ($Q$=15.9%)** · Pre-TX **+4.59±0.12 (7.1%)** · Mouse2 +2.47±0.09 (4.4%) ·
Mouse3 +2.19±0.09 (3.8%) · Mouse1 +1.48±0.05 (2.4%). **Positive 5/5 at 16–90 se.** $w$ = +0.78 to
+1.55, consistent across arms. $Q$ = share of REMAINING deviance removed — **the talk number**.
**Rungs:** M1 clone +0.07…+0.95 (clone identity buys almost nothing — the clonal signal is not a rate
shift), M2 cell **+5.3…+9.8 on the mice but NEGATIVE on Subclone/Pre-TX** (they carry a ≥100-tape
filter vs ≥20 — where QC is strict a per-cell parameter is noise; fig 3c as a held-out loss),
M3 tape-in-clone **+26…+48**, held-out $R^2$ 0.42–0.49. M4 adds **one** parameter against 3,232–64,975.

⭐⭐ **Subclone is FLAT across the floor sweep** (+8.29/+8.31/+8.31; its clones are all large, so
the floor removes 21 of 38,636 cells) — **the control showing the rise elsewhere is about which
clones are retained, not about the floor procedure**. ⭐ **$k$ has a genuine optimum near 20 even at
$n_C$ in the thousands**: Subclone +6.82/+8.31/+7.57 at $k$=5/20/50, so the fall at $k=50$ is not the
$k/n_C$ effect but the neighbour set reaching past the genuinely related cells.
**⭐ Clone floor confirms Justin's intuition, monotone 4/4**: Pre-TX +0.75→+3.05→**+4.59** at floors
20/50/100, Mouse3 +1.11→+1.50→+2.19, Mouse2 +2.13→+2.32→+2.47, Mouse1 +1.14→+1.35→+1.48.
⚠ Part is a cleaner test, part is escaping $k/n_C$ — at floor 20 with $k=20$ many clones have
$k_{\rm eff}=n_C-1$. **Rule: $k\lesssim0.2\,n_C$**, and always report $k$ with $n_C$.
⚑ **Floor 100 chosen by arithmetic**: train/param at $M_3$ = $(n_C|B|/2)/(n_C+|B|-1)$ runs 8.1 → 22.8
→ 41.5 (limit) at $n_C$ = 20/100/∞, so the steep gain is 20→200. Floor 200 deletes Pre-TX (max clone
127) and Mouse3. ⚠⚠ Clone size makes the MEASUREMENT cleaner, NOT the effect bigger — Pre-TX c7 (127
cells) gave +10.64 vs Mouse1 c36 (1,607) at +1.73.
⚠ Report per arm, never pooled: at floor 100 Subclone is 81% of retained cells.

---

**STATE (2026-09-11, end of session 10). The owed list is CLOSED and the
simulator is unblocked.** All five arms are measured on both statistics, at 100% pair coverage, and
two defects in `73`/`74` were found and fixed before the runs (details in README "Session 10").

**⚠⚠ DEFECT FIXED — a cell could be its own nearest relative.** In `74`, self and cross-clone pairs
were both marked $-\infty$, so they TIED at $+\infty$ under `argsort(-rel)` and were ordered by
index: the cell itself entered its own top-$k$ whenever its clone held fewer than $k+1$ cells. It
does **not** cancel in obs − null (the null swaps the self term for a random clone-mate). Exposure at
$k=50$: Pre-TX **71.5%** of cells, Mouse3 29.7%, Mouse1 19.1%, Mouse2 9.3%, Subclone 0.1%.
⚑ **The published Mouse2 c76 number was unaffected** (0.0% exposure) and reproduces to $10^{-15}$.
Both scripts now work **one clone at a time**, which fixes it structurally and is also 443× less
arithmetic on Pre-TX (68 min 53 s → **33 s**). `--nsub 0` = no cap; validated bit-for-bit against the
pre-fix outputs, kept in `results/preblock/`.

**⭐ Variogram, 100% coverage, obs − null lowest→top bin.** Pre-TX −0.0529→**+0.2618** ($t=95.9$) ·
Subclone −0.0102→**+0.1406** ($t=31.8$, 127.0 M pairs) · Mouse2 −0.0106→+0.0419 · Mouse1
−0.0250→+0.0369 · Mouse3 −0.0311→+0.0190. **Null flat in every bin of every arm.** Interior monotone
5/5 (dips only at the extreme bins, and only Mouse3's exceed 1 se).
⚠⚠ **CORRECTION: "Subclone is a cliff, not a gradient" was a 5%-subsample artefact.** At 100%
coverage it is **strictly monotone, 13/13 bins**; the flat middle is gone and the top bin moved
+0.1523→+0.1406. `--nsub` does not just add noise, it **reshapes** the curve — leave the cap off.
⚑ Only the SLOPE is a finding: the pair-weighted mean of obs − null equals the all-pairs reference to
4 dp in every arm, because $\sum_{c\in C} r_{cz}=0$ pins the total. The negative low bins are the
arithmetic complement of the positive high ones, so do not quote "−0.05 → +0.26" as a range.

**⭐ Prediction task, all five arms. Largest clone, nats/cell (obs − null), best $k$:** Subclone c2
(10,996 cells) **+11.29** · Pre-TX c7 (127) **+10.64** · Mouse2 c76 (3,387) +4.04 · Mouse3 c110 (210)
+3.33 · Mouse1 c36 (1,607) +1.73. **Pooled:** Subclone **+9.26** · Pre-TX **+2.86** · Mouse2 +2.65 ·
Mouse1 +0.88 · Mouse3 +0.62.
**⭐⭐ THE BEST $k$ INVERTS WITH CLONE SIZE — quote $k$ relative to $n_C$, never absolutely.** Large
clones rise with $k$ (Mouse2 c76 +1.48→+3.41→+4.04); small ones collapse (Pre-TX pooled
+2.86→+1.57→**+0.18**). Same $\gamma$ constraint: as $k\to n_C$ the neighbour mean tends to
$-r_{cz}/(n_C-1)$, an anti-signal. Pre-TX pooled belongs at $k=5$; its $k=20$ reading is an artefact
of $k/n_C$, not biology.
**⭐ The conditional table now rests on fat counts** (Subclone pooled, 38,636 cells): *where the model
says 7%, the tape is missing in **3%** of cells whose relatives all have it (n=268,916) and **90%** of
those whose relatives all lack it (n=5,715).* That largely retires the "extreme cells hold 35–345
entries" caveat. ⚠ Still quote the table, not the OR ($e^{w_f}=19.1$ for a 3%→90% span).

**✅ SD floor swept, and it does not matter** — a 20× change moves the top bin <0.5% (only the 5×
coarsening to 0.05 does anything, and it *shrinks* the effect). ⚑ Why: the floor binds on 3.8–23.4%
of entries, but those are entries where $\tilde p\approx0$ and the model is *right*, so
$r^{*}\approx-0.001$ whatever the floor. It only rescales the rare confident-and-wrong entries, and
there are too few to move a mean over millions of pairs. `--sdfloor` is now a CLI flag.

**✅ Fig 4e built** (`75_fig_variogram.py`, three standalone PNGs: variogram, conditional table,
nats-by-$k$). ⚠ Four label collisions found by rendering and inspecting — see README.

**⭐ NEXT: the simulator.** Target, per arm, through the identical scripts: the variogram's slope and
convexity (not its level), the top-bin value (+0.019 to +0.262, a **14× spread across arms** any
adequate simulator must generate), and the nats-by-$k$ curve including its inversion.

---

**STATE (2026-09-10, end of session 9). The thread has moved OFF event
counting and onto two statistics that need no events at all.** `73_dropout_variogram.py` and
`74_profile_prediction.py`. Justin's judgement, recorded: events are hard to justify and probably
useless for calibrating simulations; prevalence (1.3–25%) is a weak headline; the collapse is a
greedy heuristic that fragments 4.2× on Subclone. **Stop leading with event counts.**

**⚠⚠ THE CONFOUND both scripts exist to defeat.** Lineage relatedness is read from the edit data and
dropout decides which edits are readable — two cells that both lack tapes 1–50 look related BECAUSE
their dropout matches. ⇒ **DISJOINT TAPE SPLIT: relatedness from half A only, dropout from half B
only**, over several random splits. Never remove this.

**⚠ Pearson residuals** $r^{*}=(X-\tilde p)/\sqrt{\tilde p(1-\tilde p)}$, because a Bernoulli's
variance depends on its mean (a missing tape at $\tilde p=0.05$ is a 4.4 SD surprise; at 0.5 it is
1.0 SD). SD floored so $|r^{*}|\le10.1$ (`--sdfloor`). ✅ **Swept 2026-09-11: it does not matter** — a 20×
change moves the curve <0.5%. The floor binds on 3.78% (Subclone) to 23.35% (Pre-TX) of entries, not
the single "16%" (which was Mouse2's), but almost all of those are entries the model gets right.

**Variogram (`73`) — ⚠ SUPERSEDED by the 100%-coverage run at the top of this file.** These were
the 2026-09-10 numbers at 5–97% coverage: Subclone −0.0135→+0.1523, Mouse1 −0.0298→+0.0335, Mouse3
−0.0311→+0.0338, Mouse2 −0.0160→+0.0381, Pre-TX still running. Subclone's shape and top bin both
moved at full coverage; the others moved slightly.

**⭐ HOW TO RUN THE PREDICTION TASK (`74`)** — the digestible result, 8 s on Mouse2 c76:
`74_profile_prediction.py <arm> [--clone C | --pooled] [--nsub 0] [--nsplit 5] [--kk 5,20,50]
[--sdfloor 0.01]` — **use `--nsub 0`** (no cap; the blocked loop makes full coverage cheap), and
`--part cpu` on `submit.sh` for these, never `lesliec`.
1. split the 166 tapes into A and B; 2. for each cell take its $k$ nearest relatives IN THE SAME
CLONE from half A, **excluding itself**; 3. neighbour signals over half-B tapes, $u$ (mean Pearson
residual, for the fit) and $f$ (fraction of relatives missing the tape, for the table); 4. fit ONE
scalar $w$ with the existing predictor as a fixed offset,
$\operatorname{logit}\Pr(X_{cz}=1)=\eta_{cz}+w\,u_{cz}$, $\eta=\alpha_c+\beta_z+\gamma_{Cz}$;
5. evaluate on HELD-OUT cells → nats per cell; 6. null = permute cell labels on the B residuals
within clone.
⚠⚠ **Quote observed − null only.** The $\gamma$ fit forces $\sum_c r_{cz}=0$, so a random clone-mate
is negatively correlated with $c$ by $\approx-1/(n_C-1)$ — ~−5% in a 20-cell clone, comparable to
the signal and opposite in sign.
**Mouse2 c76: +4.04 nats/cell at $k=50$** on the full 3,387 cells (the +3.04 on the record was a
1,500-cell subsample at $k=20$). Conditional table, model prediction
held fixed: *where the model says 14%, the tape is missing in 9% of cells whose relatives all have it
and **99%** of those whose relatives all lack it.*
⚠ Quote the TABLE, not the odds ratio ($e^{w_f}=4.7$ understates a 9%→99% contrast — it is linear,
the relationship is a cliff). ⚠ On Mouse2 c76 the extreme cells hold 35–345 entries — but Subclone pooled gives the same contrast
on 1,968–398,226, so prefer it for the headline. Always print counts.

**⚠⚠ THE ORGANISING FACT — the effect is CONCENTRATED, not diffuse.** Most entries carry no lineage
signal; a minority carry an overwhelming one. That reconciles the variogram's 0.024, the event
route's $p=10^{-46}$, the 1.3–25% prevalence and the 9%→99% table. ⇒ **stratify, never average**,
when conveying magnitude. ⚠ Aggregating over tapes does NOT inflate a correlation (ceiling
$\sqrt{0.024}=0.155$); what aggregates is the likelihood gain.

**⚠⚠ COMPLETENESS IS PARTLY THE ATTRIBUTION RULE.** $\hat\pi$ median 1.000 is inflated because the
collapse is ordered by $\Lambda_{\rm hard}$, maximised at the largest COMPLETE clade. Pre-collapse,
above-threshold combos have median $\hat\pi$ 0.88–0.98. Controlling for power via
slack $=(k-k_{\min})/(m-k_{\min})$: strong in Pre-TX 0.87, Mouse2 0.97, Subclone 0.93; **equivocal in
Mouse1 0.57 and Mouse3 0.56**. Honest headline: *where there was room not to be complete, losses are
complete in three of five arms.*

**Direction 2 PARKED, not discarded** — evidential not modelling; map is many-to-one
($\bar j\approx1.6$), depletion partial (6–47%), knock-on to $q$ a fraction of a percent. The
split-half result (87.3%, 33.9× null) stands if the biology ever needs defending. Owed if resumed:
multiplicity check ($\theta$ vs $(j-1)/j$) and the out-of-sample prediction into the mouse clades.

**OWED NEXT — ✅ (1)–(4) ALL DONE 2026-09-11, see the session-10 state block at the top.**
~~(1) Pre-TX pooled variogram; (2) `74` on the other four arms and pooled; (3) sweep the Pearson SD
floor; (4) raise `--nsub`~~ — (5) **the simulator** is now the only one left, and is next.

**State (2026-09-10, end of session 8). ⭐⭐ THE TEST IS NOW ANALYTIC — `67_exact_merge.py`.**
No calibration draws, no clade-size strata, no permutations needed at all (they are kept only as a
validation set). Read README "Exact permutation moments and analytic p-values" before touching event
calling. Per (clade $S$, tape $z$) with $r_c=X_{cz}-\tilde p_{cz}$ and score $T=k-E$:
the $\gamma_{C,z}$ fit forces $\sum_{c\in C}r_c=0$, so a clade is a simple random sample **without
replacement** from its clone's residuals, giving the closed form
$\mathbb{E}_\pi[T]=m\bar r$ and $\mathrm{Var}_\pi[T]=m\frac{n-m}{n-1}\sigma^2$
(population = the clone's cells **in that block**, not the whole clone). Then
$p=\max(p_{\rm normal},p_{\rm hypergeometric})$ and a single Bonferroni cut
$p\le\mathcal{E}/N_{\rm test}$.
⚠⚠ **CORRECTION: the model variance $V$ is wrong by ~2×** — the permutation SD of $z=(k-E)/\sqrt V$
is 0.51–0.67, not 1. My Monte Carlo "verification" simulated the model null with $\tilde p$ FIXED,
i.e. assumed what it tested. Calling was unaffected (the margin cancels a constant factor) but every
reading of $z$ as "standard deviations" was wrong.
⚠⚠ **Edgeworth was dropped before building** — at $z=6$ its correction is 37× the leading term, so
the expansion diverges where decisions are made. The hypergeometric replaces it (exact in the
constant-$\alpha$ limit; Fisher's exact test on clade × missingness).
⚑ **The strata are gone and are not needed**: $\mathrm{Var}_\pi$ depends on $m$ and $n$ explicitly,
so $p$ already conditions on clade size. ⚑ **Testability is structural**: $\sigma^2>0$ and $m<n$
(53.3% on Mouse3) — the rest are clades that ARE their whole block population, where $T\equiv0$.
**Mouse 3: 92 events / 2.97% at $\mathcal{E}=1$; 128 / 3.79% at $\mathcal{E}=12$**, and at matched
stringency the three routes converge (128 / 127 / 132, 85 clone×tape pairs common to all three).
Verified: predicted/empirical perm SD ratio 1.046, $p_{\rm norm}/p_{\rm emp}=2.31$ and
$p_{\rm hyper}/p_{\rm emp}=13.0$ (both conservative). Runs in **9 s / 1.76 GB** on stored parts.
⚠ Owed: the normal tail is validated only to $p\approx10^{-3}$ while events sit near $10^{-9}$;
Bonferroni over dependent tests over-corrects; the hypergeometric's 13× conservatism costs real power.

**State (2026-09-09, end of session 7). ⭐ THE DETECTOR WAS REBUILT — read README
"The detector, rebuilt" before touching event calling.** An event is now defined on **two axes**:
detect with the one-sided **score test** $z=(k-E)/\sqrt V$ (nested family
$X_c\sim\mathrm{Bern}(\sigma(\eta_c+\delta))$, $\delta=0$ **is** $H_0$), attribute with
$\Lambda_{\rm hard}$ (maximised at the largest clade that is still complete = the Dollo rule), and
**report $\hat\pi=k/m$** rather than assume it. Mouse 3: **127 events / 3.19% of missing**,
$\hat\pi$ median **1.000** vs $\bar e=0.376$ predicted — the first NON-CIRCULAR completeness claim,
and it converges on the committed routes (128/3.09%, 132/3.15%).
⚠⚠ **$\Lambda_{\rm soft}$ is NOT a completeness-agnostic $\Lambda_{\rm hard}$** — it is non-nested
against $H_0$, decomposes as elevation + *dispersion-of-$\tilde p$*, and 72% of its Mouse3 calls were
large clades with $z<0.5$ (no elevation at all). Do not resurrect it as a detector.
⚠⚠ **Two older numbers are superseded**: the partial fraction is **29.9%**, not 8.9% (the hard route
cannot see partials, so its figure was circular); and the **$\hat\pi$-by-depth gradient is gone**, so
the "graded losses are largely clade coarseness" reading from 09-03/09-07 was a detector artefact.
⚠ Calibration draws are **appended** (`62 --calib N`), not held out, so their number is a post-hoc
decision costing one scan each. NCAL is the only source of precision on the null count.
**Fig 5 built (`66`, panels a/b/c standalone).** Panel a is **faceted by clade-size stratum**, one
threshold line per facet — that line *is* the rule, and every called event sits above it. ⚠ Do not
re-pool it: the earlier single "null ceiling" mixed clade sizes and put 58 of 127 events below a line
they were never tested against.
**Only Mouse 3 has been run on the new detector.** The other four arms are the same code
(`64_submit_percombo_soft.sh <arm>`; Pre-TX 24 G, Subclone 16 G).

**State (2026-09-03, end of session 4).** **Fig 3 done (a–e)**; **row A9 measured** — tape loss is
heritable, heritable *below* the clone (monotone subclade-depth gradient, 5/5 arms), and discrete
rather than a graded rate. Fig 4 panels a/b built. Event catalogue, clone-wide layer, soft variant
and the `MAX_D=6` test all done. **$B=1{,}000$ permutations launched across all 15 configurations**;
clone-wide layer already back at $p=1/1001$ in 5/5 arms, the rest pending the `perm_collect` job.
See README "Row A9", "Panels c/d groundwork", "B = 1,000 permutations — launched" and the full
$\rho$ derivation under "Methods, in full".

**⭐ NEXT, and it is the crux.** The effect is claimed to be a **tree** property, not a clonal one,
so the sharpest test is that it appears **WITHIN a single large clone** at subclade resolution —
where a clone-level or batch effect cannot reach. Start with **Mouse2's 3,387-cell clone**: it has
near-zero *clone-level* excess yet holds 98.9% of that arm's pooled weight, so if tree structure is
real it must be visible inside it. Then Subclone's 10,996-cell clone. Panels **c and d may carry the
figure alone**; a and b become supporting.

**⚠ Estimator rules, learned the hard way.**
0. **Every layer needs its own margin.** $\Lambda$ scored against $\alpha_c+\beta_z$ alone gave
   FDR 64–73% at *every* threshold, because clone-wide losses look spectacular in permuted data too.
   Adding a per-(clone,tape) $\gamma$ dropped it to ≈0%. Whatever layer you are testing, fit the
   layer above it.
1. **Equal-clone weighting, not pooled.** $\hat\rho$ is a ratio of sums, so clones enter weighted by
   $n_C(n_C-1)$ — one clone held 98.9% of Mouse2 and 84.0% of Mouse1, and pooling *hid* the signal
   (Mouse2 +0.009 pooled vs +0.246 equal-weighted).
2. **Significance saturates; report effect size.** With millions of within-clone pairs, 77–100% of
   tapes reach FDR 5% and median $z$ runs +8 to +422. $z$ is an ordering statistic only — the
   per-tape null is right-skewed. Dropout characterised on both margins — see
README "Fig 3 redesign". Compatibility measured on Mouse3 and Initial;
Mouse1/Mouse2 never finished (5 OOMs, all on their single largest clone) and are now **moot** —
see the strategic reassessment in README.md. $C$ measured on Mouse3 (0.42–0.57 depending on
heuristic). Subclone ground-truth test passed.

**⚑ Direction has changed.** The skeleton is a step sideways: its output is arbitrary (45% spread
across heuristics), only 9 of 2,547 clones exceed 1,000 cells so the scale problem it solves barely
exists here, and Felsenstein pruning handles dropout natively. **Current work is a figure programme
motivating the likelihood route**, not further skeleton development.

**Analysis unit (settled §D.4c).** Per `ClonalBC` clone, cells = clonal-barcode table ∩ group edit
table, **with the paper's tape filter applied (≥100 recovered tapes for Initial/Subclone, ≥20 for
Mouse1–3) — it is NOT applied in the delivered CSVs.** Clone-size floor is ours to choose, for
statistical power, not inherited from the paper.

## Figure programme (for PI presentations)

| fig | status | script |
|---|---|---|
| 1 recorder / $q$ | ✅ done, 2 panels | `20_fig1_recorder.py` |
| 2 homoplasy | ✅ done, **two versions** — `simple` (present) and `mle` (reserve) | `21_fig2_homoplasy.py <simple\|mle>` |
| 3 dropout | ✅ **done, 5 panels a–e**, each a standalone PNG | `27` a · `28` b · `30` c · `31` d+e |
| 4 dropout & lineage (row **A9**) | ⚑ a/b built; catalogue + `MAX_D=6` test done; **c/d now built — see the next row** | `32`–`44` |
| 4d **excess curve** — observed vs a closed-form null, 5 arms | ✅ built (`fig4d_excess.png`); 11,119–499,263× at each arm's Bonferroni cut | `72` |
| 4e **variogram + prediction** — structured dropout with no events | ✅ **built**, three standalone PNGs (`fig4e_a/b/c_*`), all five arms at 100% coverage; the simulator calibration target | `73`,`74`,`75` |
| 4c **the worked example** — what a loss looks like, how the test works, and a near-miss | ✅ **built**, three standalone PNGs (`fig4c_a/b/c_*`), all in Mouse2 clone 76 | `69` |
| 4c/d **the detection plane** — what a silencing event IS; fulfils the planned c/d | ✅ **built on Mouse 3**, three standalone PNGs (`fig5a/b/c_*`); panel a faceted by clade size, one threshold line per facet | `62`,`63`,`66` |
| 6a **the ladder** — held-out gain against parameters spent, per rung | ✅ **built**, five arms | `77`,`80` |
| 6b **$Q$** — the lineage share of what the technical model left | ✅ **built**, five arms | `77`,`80` |
| 6c **per-tape unanimity** — missing vs present against the diagonal | ✅ **built**, five arms × $k$=3,5 | `83`,`84` |
| 6d **the companion** — capture matched, relatedness not | ✅ **built**, five arms | `83`,`85` |
| 5 compatibility spread + homoplasy null | needs the simulator | — |
| 6 method comparison under simulation | planned; needs simulator | — |
| 7 calibration / honest uncertainty | planned; needs simulator | — |
| 8 runtime–accuracy frontier | planned | — |

**No composed grid.** Justin does not want a–e assembled into one multi-panel figure; the panels are
presented individually. So `27/28/30/31` each own their PNG and there is no assembly script —
do not build one.

**Fig 3 as built (2026-09-01).** The old figure described *our handling* of `None`; the redesign
measures the assay. The original decomposition panel was dropped outright — nobody disputes that
trailing `None` is biology, and if challenged the answer is the flat per-cell depth profile (panel e),
not the internal-gap rate, which only excludes *random* site loss and not terminal truncation.

| panel | claim | headline |
|---|---|---|
| a | dropout is not a coin flip on entries | VIF 22.6, $\rho_{\rm cell}=0.131$; $R_c$ bimodal |
| b | which tape it is matters more than which cell | $\rho_{\rm tape}=0.250$; rates 0.006–0.962 |
| c | the shelf is a QC choice, not an arm difference | $\rho_{\rm cell}$ 0.13–0.18 → 0.012–0.024 at a common ≥100 cut |
| d | a tape we cannot see recorded less | $\rho=+0.34$; deciles 2.92 → 4.95 sites |
| e | …but a cell we see badly recorded as much | $\rho=+0.05$; 4.88 → 4.98 sites |

d and e **must keep their shared $y$-axis** — the contrast is the argument and it dies if they are
scaled independently.

**Figure conventions.** Light theme, palette from the `dataviz` skill (surface `#fcfcfb`, ink
`#0b0b0b`/`#52514e`, grid `#e1e0d9`). **Orange `#eb6834` means "the null" figure-wide** and is not
available as a series colour. Five-arm series use categorical slots 1,3,4,5,7 —
Mouse1 `#2a78d6` · Mouse2 `#1baf7a` · Mouse3 `#eda100` · Pre-TX `#e87ba4` · Subclone `#4a3aa7`
(validated: CVD $\Delta E$ 9.1, normal-vision 19.6 on the adjacent pairlist). Where a panel needs
both five arms and a null, the **null is drawn neutral** (`#cfcec6`/`#8f8e86`), because orange fails
the normal-vision floor against magenta and yellow.
⚠ The skill's `validate_palette.js` will not run on the cluster (node v10.24: no ESM, no `??=`).
Port it to Python — same thresholds, same Machado–Oliveira–Fernandes matrices — rather than eyeballing.
Minimal annotating text — Justin's explicit preference. Panels should **teach the quantity**
(show the formula / the read-off), not just report it. Always render and inspect the PNG for label
collisions before calling a figure done.

## Conventions here
- Scripts in `src/`, figures in `figures/`, small tables in `results/`.
- Record findings in `README.md` as you go, not at the end.
- Scripts that only stream the CSVs use stdlib `csv` and run under any python.
- Numerics/figures need the project env: `/data1/choij10/justin/envs/pando/bin/python`
  (python 3.12 + numpy/scipy/matplotlib, built from the system `miniforge3` module).

---

# B = 1,000 permutations → proper p and q values — **COMPLETE 2026-09-04**

**Result.** All 15 configurations pooled to exactly 1,000 permutations; every merge passed its
completeness assertion. **$p=1/1001$ in all twenty runs.** ⭐ Lead with the **null maximum**: the
most extreme of 1,000 random labellings gave **111** candidates in Subclone against **279,973**
real (118–2,522× across arms). Null is right-skewed 6–100×. Every hard-catalogue event clears
$q\le1.9\times10^{-4}$. Cell × tape entries inside called blocks: 2,916–266,688 per arm.
⚠ Disclose the one marginal result — clone-wide Mouse 2, FDR 4.36%, max $q$ 0.044.
⚠⚠ **The thresholds did NOT move** (Mouse 3 soft 16.3→16.2; all others already on the 10-nat floor):
the $B=3$ null was unbiased in the mean and only missed the tail. Full tables in README
"B = 1,000 — results".

**✅ THE CATALOGUE TO QUOTE (2026-09-07): `*_lam2_b2` / `*_lam4_b2`.** Floor 2 (soft 4),
clade-size-stratified, threshold per stratum at expected false $\le2$. Subclone 6,597 events /
20.07% of all missing · Pre-TX 2,019 / 2.46% · Mouse2 452 / 7.02% · Mouse1 424 / 4.34% · Mouse3
132 / 3.15%. Only 1.1–1.8× floor 10 (0.8× Subclone): the 3–14× at a 5% *rate* was false positives,
because FDR is computed on candidates while events are what we report. Use `--budget`, not
`--target`. Two corrections to my own reasoning are recorded in README — the cliff lives at 5.9–12.4
nats (so floor 10 was accidentally near-right), and the 4–5 cell stratum is live, not dead.

**⚠⚠ Before any "how widespread" claim, read README "How events are actually called".** The
10-nat cutoff is a hard-coded **scan floor**, and the adaptive FDR$\le5\%$ step is **degenerate in
13 of 15 configurations** (FDR at the floor is 262–2,511× below target for the hard catalogue, so
`LAM` snaps to the floor; it binds only for Mouse 3 soft). ⇒ counts and the dropout share are set by
an arbitrary constant, and the **arm ranking flips** with it (Pre-TX 1,783 events at $\ge10$, 110 at
$\ge20$). Quote a curve, never a number. ⚑ The floor is conservatism, not rigour — the numbers are
loose lower bounds. ⚠ Do not inherit the 5% target; prefer an absolute one.

**Launch state (2026-09-03).** Implemented and submitted across everything: 5 arms × {`42` clone-wide, `40` hard
sub-clone, `43` soft at depth 4, `43` soft at depth 6} = 15 permutation configurations + 5 direct
clone-wide runs. 300 array tasks, all `RUNNING` from the moment of submission; a dependent
`perm_collect` job pools the parts and re-runs each observed scan once with `--nullfile`.
**Driver: `src/47_submit_B1000.sh`** (re-runnable) → `src/48_collect_B1000.sh` → per-config
`46_perm_merge.py`. Design, the count-vector trick, and the ⚠ correction to the monotonisation
direction are written up in README "B = 1,000 permutations — launched".

**First results in (the five `42` clone-wide runs):** $p = 0.000999 = 1/1001$, the floor, in every
arm at both the chosen threshold and the 10-nat scan floor — i.e. **not one of 1,000 permutations
came close**. Loss counts are unchanged from $B=200$ (Mouse 1 725, Mouse 2 410, Mouse 3 398), so
$B$ bought resolution on $p$, not a different answer.

Motivation, for the record: a presentation audience expects a significance test, and our
permutation p-values were floored at $1/(B+1)$ with $B=3$–200. The counts were overwhelming
(28,367 candidates vs 1.4 null) but the *formal* claim was only $p<0.17$ for the catalogue.

---

## ⭐ NEXT, once the B=1,000 run lands — two directions, written up in full in README
## "Two directions from the B = 1,000 run"

**A. How widespread is heritable silencing?** Lead with **cell × tape entries inside called
blocks**, not event counts — events are an artefact of how the clade search carves the tree.
⚠⚠ **$B$ does not lower the detection floor.** Three distinct knobs: $B$ (p resolution), the FDR
*threshold* (will fall where a 3-permutation null had pushed it up — Mouse 3 soft 16.3 nats should
drop toward 10), and the scan **floor** `--lam`, hard-fixed at 10 nats, below which nothing is ever
collected and which this run's grid cannot see beneath. ⇒ **the follow-up is `40`/`43` with
`--lam 4`**, which subsumes the floor-10 run; script `45` already shows real signal there
(4.6× enrichment at $\ge4$ nats on Mouse 3). Deliberately *not* done by restarting the current
jobs — read this run's FDR curve first, then size it. Expect a chunk of new low-$\Lambda$ events to
be clade coarseness; the $\hat\pi$-by-depth stratification in `43` is the control.

**A2. ⚠⚠ Per-combo permutation p-values — asked 2026-09-03, answered: not at this floor.**
Well-defined (clade slots persist under the permutation) and cheap (one exceedance *counter* per
combo, ~22 MB, never $B$ values per combo). But at the 10-nat floor it is strictly weaker than the
pooled count: both spend $B=1{,}000$, yet pooling buys $B\times N_{\rm combos}\approx10^9$ null
draws (tail resolution $10^{-7}$) against 1,000 (resolution $10^{-3}$). Measured: 8 permutations of
Mouse 3's hard scan gave **zero** null candidates above 10 nats, so every real event censors at
$1/1001$ — the statistic saturates where the signal is strongest. ⚠ And a second floor no $B$ can
lift: an $m$-cell clade in an $n_C$-cell clone has only $\binom{n_C}{m}$ permuted compositions, so
$p\ge1/\binom{n_C}{m}$ — a 4-cell clade in a 6-cell clone **cannot reach $p<0.05$**, which is a
structural mechanism for "no events in clones <20 cells". ⇒ Fold into the `--lam 4` run instead:
lower floor + **null stratified by clade size** (the real repair for "weaker but real") + per-combo
counters, which stop being censored at $\Lambda\ge4$. Full argument in README.

**B. ⭐ Does the co-integrated symbol vanish when a tape is silenced?** The first *orthogonal* test
of the mechanism — everything so far infers silencing from missingness, which is what technical
dropout also looks like. pegRNA and tape share one cassette and pegRNAs act in **trans**, so
silencing integration $z$ should remove symbol $s(z)$ from **every other tape** in those cells.
⚠ The $z\mapsto s(z)$ map is unknown (`TargetBC` 10-nt vs `NNNN` 4-nt, never linked), so **do not
assume it — recover it**, and let its structure be the test: five-way cross-arm concordance (all
**166 TargetBCs verified identical across arms**, 2026-09-03), near-injectivity with the collision
rate predicted for 166 draws from 256 ($E=122$ distinct vs 100–106 observed),
$\mathrm{corr}(\beta_z,\xi_{s(z)})>0$.
**⚠ Step 0 before anything: measure the fraction of (tape, site) slots polymorphic within a clone
and within a clade.** Tapes are ~4.5–5/6 saturated, so most content may predate the loss and only
post-loss insertions can show depletion — that fraction is the power calculation for the whole idea.
Then 3–5 hand-inspected examples, then screen the **clone-wide** layer first (strongest: 4,763
losses over 1,188 Pre-TX clones), then the sub-clone version.

## ⚠ First, the permutation itself — state it correctly

We do **not** shuffle clone labels for the sub-clone test. The rule is:

> **Permute the label you are testing, blocked by the label above it.** Whole cell ROWS move
> (a cell's entire 166-tape vector travels together); the labels stay fixed at their positions.

| test | shuffle | blocked within | preserves | destroys |
|---|---|---|---|---|
| **clone-level** (`34`, `37`, `38`, `42`) | clone membership | **harvest sample** | cell profiles, tape marginals, clone sizes, clone→sample composition | which cells are in which clone |
| **sub-clone** (`35`, `36`, `40`, `43`, `45`) | subclade membership | **clone** | all of the above **plus each clone's own rate for every tape** | which cells are in which subclade |

That last preservation is why clone-wide losses are invisible to the event catalogue — they survive
the within-clone permutation untouched — and why they needed script `42` with the within-sample
permutation instead.

Implementation (already in the scripts, do not redesign): sort cells by the blocking label, permute
indices inside each block, invert the map, and index the data arrays with it. Labels and prefix
codes are never touched.

**And yes — the entire scan is redone per permutation** (all anchors × depths × clades × tapes).
That is what makes $B=1{,}000$ expensive.

## Measured per-scan cost (elapsed ÷ (nperm+1), from `sacct` 2026-09-03)

| script | Mouse3 | Mouse1 | Mouse2 | Initial | Subclone | **B=1,000 single job** |
|---|---|---|---|---|---|---|
| `42_clonewide` | ~0.02 s | ~0.03 s | ~0.03 s | ~0.12 s | ~0.12 s | **~2 min — just run it** |
| `40_event_catalogue` | 2 s | 4 s | 5 s | 40 s | 39 s | mice ~1 h; Initial/Subclone **~11 h** |
| `43_soft_events` | 9 s | 27 s | 23 s | 3.2 min | 3.1 min | mice ~7 h; Subclone **~52 h** |

⇒ `42` needs nothing. `40` and `43` need **parallelisation**.

## The design to implement

1. **Add `--permpart i/N --seed s`** to `40_event_catalogue.py` and `43_soft_events.py`. Each task
   runs $B/N$ permutations with an independent seed and writes **only the threshold-count vector**
   (counts of candidates ≥ each value on the existing $\Lambda$ grid) to
   `results/permcounts_{script}_{arm}_p{i}.json`. Do not store candidate lists — the counts are all
   that the FDR curve needs, and they are a few hundred integers per permutation.
2. **Submit N = 20 parts** per arm. Subclone soft becomes ~2.6 h per task. Right-size from the table
   above; do **not** use the 128 G / 12 h defaults (every oversized ask this project made sat in `PD`
   behind jobs that would have run immediately).
3. **A merge script** (`46_perm_merge.py`) pools the parts and emits:
   - **global permutation p-value** — $p = \frac{1+\#\{b:\,C_b \ge C_{\rm obs}\}}{B+1}$ where $C$ is
     the total candidate count above the chosen threshold. With $B=1{,}000$ and zero exceedances
     this licenses **$p < 0.001$**.
   - **per-event q-values** — $q(\Lambda) = \frac{\text{mean null count}(\ge\Lambda)}{\text{observed count}(\ge\Lambda)}$,
     made monotone by a running minimum from the top, then attached as a column to every row of
     `events_{arm}.tsv.gz` and `clonewide_{arm}.tsv.gz`.
4. **Re-run the observed scan once** (cheap) to attach $q$ to the event rows.

## ⚠ What the p-value does and does not establish — say this in the talk

It rejects exactly one null: *dropout is exchangeable among cells within a clone*. It does **not**
by itself say "silencing happened" — any lineage-correlated effect, technical or biological, breaks
exchangeability too. The technical explanation is excluded separately by **capture-independence**:
event cells have median $R_c$ = 116 against 114 for all cells, so they are if anything *better*
captured than average. Lead with the comparison ("1.4 null events vs 28,367 real"), not the p-value —
it makes the null explicit, which a p-value hides.

## Process gotchas already paid for

- **`cd` to the repo root before touching `notes/` or the root `CLAUDE.md`.** Two commits this
  session silently dropped those edits because the shell was left in this directory.
- **Assert on every `str.replace`** into a notes file — a no-op replace prints success otherwise.
- **Right-size from a *completed* job's `sacct`**, never a mid-run `sstat`.
