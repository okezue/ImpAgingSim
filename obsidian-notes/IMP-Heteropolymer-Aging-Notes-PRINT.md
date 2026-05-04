---
title: "IMP Heteropolymer Aging — Print Version"
tags: [glass-physics, heteropolymer, aging, IMP-model, print]
date: 2026-04-20
cssclasses:
  - print
---

<style>
@media print {
  body { font-size: 11pt; line-height: 1.45; }
  h1 { font-size: 22pt; margin-top: 0; }
  h2 { font-size: 16pt; margin-top: 1.2em; border-bottom: 1px solid #000; padding-bottom: 4px; }
  h3 { font-size: 13pt; margin-top: 1em; }
  table { font-size: 10pt; }
  pre, code { font-size: 9.5pt; }
  blockquote { border-left: 3px solid #000; padding-left: 10px; }
}
.pgbreak { page-break-after: always; break-after: page; }
.term { font-weight: bold; }
.def { margin-bottom: 0.6em; }
</style>

# IMP Heteropolymer Aging

**Thesis.** In a frustrated heteropolymer, **how the random couplings are correlated along the backbone matters as much as how strong they are**. This opens a new axis of glassy-aging behavior even when the disorder amplitude $\epsilon$ is held fixed.

**Reading order.** Study summary → Glossary → Big picture → Model → Diagram → Ensembles → Protocol → Observables → Findings → Interpretation → Implications → Next steps.

---

## Study summary (plain English)

We built two versions of the same polymer that were statistically indistinguishable one bond at a time — every single random pairwise coupling had the same average and the same spread — but differed in whether the couplings had hidden patterns linking them to a backbone sequence. One version was the standard "every pair is its own independent coin flip." The other was built so that beads sharing the same "flavor" along the chain got coordinated pulls toward each other, and beads of opposite flavor got coordinated pulls apart, while keeping the overall strength of the randomness identical. We then heated both polymers until they were loose and wandering, suddenly dropped the temperature (the **quench**), waited different amounts of time, and watched how the contact patterns forgot themselves, how much the different simulated copies disagreed with each other, and how unevenly the beads moved. Using the zero-disorder case as a built-in sanity check (where both versions *must* agree, and did), we found that at every nonzero disorder strength the two versions gave systematically different dynamics: the correlated version actually **forgot its contacts faster, had up to 44% less trajectory-to-trajectory disagreement at late ages, showed its characteristic "heterogeneity time" nearly three times earlier, had a broader spread of relaxation timescales, less bead-level hopping intermittency, and slightly tighter local caging** — and sweeping the correlation strength continuously showed these effects tune smoothly and can even flip direction depending on age.

==**Why this is not an obvious conclusion.**== The naive expectations point the wrong way on almost every axis. If you'd asked a physicist beforehand, most would have said structured, correlated couplings should make the polymer *more* rigid, *more* glassy, and *slower* to decorrelate — because coherent pulls sound like they'd lock the chain into a preferred shape. The opposite happened for contact decorrelation. Second, the entire heteropolymer-aging literature for 30+ years has tacitly assumed that once you fix the *strength* of the disorder, you've essentially fixed the physics — the variance of the random couplings is usually treated as the only knob. Finding that two ensembles with identical variance give quantitatively and qualitatively different glassy behavior says that variance is not enough, which quietly undermines a default assumption in the field. Third, the result doesn't fit a "correlated disorder is just effective-weaker (or effective-stronger) disorder" story: at intermediate ages the ordering between the two ensembles actually reverses, which no single-parameter rescaling can explain. Fourth, the way different observables split is genuinely counterintuitive — correlated disorder simultaneously makes the *contact network* more flexible while making *individual beads* more locally trapped, a tension you wouldn't predict from any one-sentence intuition. And fifth, the suppression of dynamical heterogeneity doesn't mean "less glassy" — the spread of relaxation timescales actually *broadens* under correlations, so the system is in some senses *more* glassy while in others *less* so. Put together, the result says something subtle: correlations don't uniformly soften or harden the landscape — they **reshape its topology**, reducing the number of mutually-incompatible deep basins (which is what drives trajectory disagreement in uncorrelated disorder) while sharpening the basins that remain. That's a structural statement about energy landscapes, not a quantitative tweak, and it's the kind of thing that's easy to miss if you only ever simulate the i.i.d. default.

<div class="pgbreak"></div>

## 0. Glossary

**Quench.** A sudden drop in temperature. In simulations, instantly change the inverse temperature $\beta = 1/(k_B T)$ from a small value (hot) to a large value (cold). After a quench the system is out of equilibrium.

**Aging.** The dynamics depend on how long it has been since the quench, not just on the current time. The "age" is the waiting time $t_w$. Aging breaks time-translation invariance.

**Heteropolymer.** A chain of beads where the beads are not all the same — different chemistries, different interactions. Examples: proteins, intrinsically disordered proteins.

**Frustration.** No single configuration can simultaneously satisfy every local preference. Three Ising spins on a triangle each wanting to anti-align is the canonical example.

**Quenched disorder.** The random coupling $\eta_{ij}$ is drawn once when the polymer is built and never changes. Contrast: annealed disorder, which fluctuates.

**Glassy dynamics.** Slow, history-dependent, multi-timescale relaxation. Hallmarks: aging, subdiffusion, dynamical heterogeneity, broken time-translation invariance, stretched-exponential correlations.

**Dynamical heterogeneity.** Different parts of the system relax at very different rates. Trajectories from the same starting configuration disagree about how the system evolves. Quantified by $\chi_4$.

**Subdiffusion.** Motion slower than normal diffusion. Normal: $\langle \Delta r^2 \rangle \sim t$. Subdiffusion: $\langle \Delta r^2 \rangle \sim t^\alpha$ with $\alpha < 1$.

**Ensemble.** A recipe for drawing the random couplings $\eta_{ij}$. We compare two: i.i.d. Gaussian, and a sequence-structured correlated ensemble.

**Marginal variance.** The variance of a single random variable. Both ensembles here have $\mathbb{E}[\eta_{ij}^2] = 1$, so any difference must come from correlations *between* different $\eta_{ij}$'s.

**Lennard–Jones potential.** $V(r) = R/r^{12} - A/r^6$. Repulsion at short range, attraction at long range. Minimum at $r_{\min} = (2R/A)^{1/6}$.

**Metropolis Monte Carlo.** Sampling algorithm. Propose a small random change, compute $\Delta H$, accept with probability $\min(1, e^{-\beta \Delta H})$.

**KWW function.** Stretched exponential $f(t) = \exp[-(t/\tau)^\beta]$. $\beta = 1$ is normal exponential. $\beta < 1$ signals a broad distribution of relaxation timescales.

**Four-point susceptibility $\chi_4$.** Measures how much trajectories disagree with each other. Its peak height $\chi_4^\star$ is the standard glass-physics diagnostic for dynamical heterogeneity.

**Inverse temperature $\beta$.** $\beta = 1/(k_B T)$. Large $\beta$ = cold; small $\beta$ = hot.

**MC sweep.** One unit of simulation time, defined as $N$ attempted single-bead moves.

<div class="pgbreak"></div>

## 1. Why this matters

Glasses, gels, and frustrated polymers all show **aging**: their dynamics slow down the longer you wait after preparing them. The standard minimal model in this space is the IMP heteropolymer, which puts a random attractive coupling between every pair of monomers. Everyone uses the simplest random ensemble — independent Gaussians, no correlations between pairs. Real biological sequences have correlations along the backbone. Does that matter for the dynamics?

We find: **yes, dramatically**. The correlation structure of the disorder is an independent control knob, comparable in importance to the disorder amplitude.

Three things to keep in mind throughout:

1. Glasses age. A freshly-quenched glass relaxes fast at first; the older it gets, the slower it relaxes.
2. Heteropolymers are minimal models for glasses — a chain of beads with mismatched, frustrated interactions captures the essential ingredients.
3. The random couplings $\eta_{ij}$ have a structure we choose. Standard practice has been i.i.d. We push back on that.

**Two meanings of "correlation" — keep them separate.**

- **Input correlation** — built into the random couplings $\eta_{ij}$. The i.i.d. ensemble has none; the correlated ensemble injects some via the backbone field $\sigma_i \sigma_j$. **This is the new control knob we tune.**
- **Output correlation** — what $\chi_4$ diagnoses. It correlates the overlap $Q$ across trajectories to measure how much they *disagree* — i.e., dynamical heterogeneity. This is what we *observe*.

The whole study is: turn the input knob and watch what happens to the output. Side motivation — real protein sequences have input-correlation, so the i.i.d. default that 30+ years of literature has used is a convenience nobody had tested.

---

## 2. The IMP heteropolymer

A 3D off-lattice bead-spring chain of $N=30$ monomers at positions $\{\mathbf{r}_i\}_{i=1}^{N}$, after Iori, Marinari, and Parisi (1992). The energy of any configuration $\mathbf{r}$ given a frozen disorder field $\eta$ is

$$H(\mathbf{r}\mid\eta) \;=\; \sum_{1 \le i < j \le N} \left[\, h\,\delta_{j,i+1}\,r_{ij}^{2} \;+\; \frac{R}{r_{ij}^{12}} \;-\; \frac{A}{r_{ij}^{6}} \;+\; \frac{\sqrt{\epsilon}\,\eta_{ij}}{r_{ij}^{6}} \,\right]$$

with $r_{ij} = \|\mathbf{r}_i - \mathbf{r}_j\|$ and $\delta_{j,i+1}$ the Kronecker delta (only fires for backbone neighbors).

| Term | Symbol | Role |
|------|--------|------|
| Harmonic spring | $h\,\delta_{j,i+1}\,r_{ij}^2$ | Holds adjacent beads together. Only on when $j = i+1$. |
| Hard repulsion | $R/r_{ij}^{12}$ | Beads can't overlap. Blows up as $r \to 0$. |
| Soft attraction | $-A/r_{ij}^{6}$ | Standard van-der-Waals tail. |
| Quenched random coupling | $\sqrt{\epsilon}\,\eta_{ij}/r_{ij}^{6}$ | The hero. Each pair has its own random bias. $\epsilon$ controls strength. |

The first three terms together are a clean Lennard-Jones polymer. The fourth term **adds the disorder**: some pairs become extra-attractive ($\eta_{ij} > 0$), others extra-repulsive ($\eta_{ij} < 0$). Because the chain is a connected object in 3D, **it can't satisfy all of these preferences at once**. That's the **frustration**.

**Parameters.** $A = 3.8$, $R = 2.0$, $h = 1.0$, $N = 30$. LJ minimum: $r_{\min} = (2R/A)^{1/6} \approx 1.007$. Contact cutoff: $r_c = 1.25\,r_{\min} \approx 1.259$. Two non-bonded beads ($|i-j|>1$) count as "in contact" when $r_{ij} < r_c$.

<div class="pgbreak"></div>

## 3. Diagram

```
       σ = (+ − − + + − +)    ← backbone flavors

           ╭─── η₁,₄ ───╮
           │            │
           ▼            ▼
       (+)─(−)─(−)─(+)─(+)─(−)─(+)
        1   2   3   4   5   6   7
                ▲           ▲
                ╰── η₃,₆ ───╯

   ───  Backbone bond:  h·r²ᵢ,ᵢ₊₁
   ─ ─  Quenched LJ + ηᵢⱼ on every non-bonded pair
        (only two ηᵢⱼ drawn; rest omitted)
   (+)/(−) = σᵢ ∈ {+1, −1}
```

Solid lines are the backbone — every adjacent pair is glued together by a harmonic spring. Dashed arrows show example non-bonded couplings $\eta_{ij}$. In reality, every non-adjacent pair has one. For the $N = 7$ chain drawn above that's $\binom{7}{2} - 6 = 15$ random couplings; for $N = 30$, it's $\binom{30}{2} - 29 = 406$.

The $(+)$ and $(-)$ are the **backbone flavors** $\sigma_i$. In the i.i.d. ensemble these are irrelevant — every $\eta_{ij}$ is its own independent Gaussian. In the **correlated ensemble**, like-flavor pairs ($\sigma_i \sigma_j = +1$) get a positive bias and unlike-flavor pairs get a negative bias.

After the quench, the chain crumples into a globule. The random $\eta_{ij}$'s decide which non-bonded contacts get rewarded and which get punished. **Frustration → rugged energy landscape → glassy aging.**

<div class="pgbreak"></div>

## 4. The two disorder ensembles

Both ensembles have **identical mean (0) and identical variance (1)**. The only difference is in the correlation structure between *different* $\eta_{ij}$'s.

**(A) i.i.d. ensemble (standard).**

$$\eta_{ij} \sim \mathcal{N}(0,1) \quad \text{independently for each pair } (i,j)$$

Every pair gets its own independent draw from a standard Gaussian.

**(B) Correlated ensemble (the new construction).**

First, generate a backbone "flavor" sequence $\sigma_i \in \{+1, -1\}$ from a two-state Markov chain with persistence $\pi$:

$$P(\sigma_{i+1} = \sigma_i) = \pi$$

The autocorrelation along the backbone decays geometrically:

$$\langle \sigma_i \sigma_{i+k}\rangle = (2\pi - 1)^{|k|}$$

Then construct the pair couplings as a weighted mix of structured and noise:

$$\eta_{ij} \;=\; \kappa\,\sigma_i\sigma_j \;+\; \sqrt{1 - \kappa^2}\,\xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0, 1)$$

where $\kappa \in [0, 1]$ controls how much of $\eta_{ij}$ comes from the structured part vs. noise.

The $\sqrt{1-\kappa^2}$ weighting preserves the marginal variance:

$$\mathbb{E}[\eta_{ij}^2] \;=\; \kappa^2 \underbrace{\mathbb{E}[(\sigma_i\sigma_j)^2]}_{=\,1} + (1-\kappa^2)\underbrace{\mathbb{E}[\xi_{ij}^2]}_{=\,1} + 2\kappa\sqrt{1-\kappa^2}\underbrace{\mathbb{E}[\sigma_i\sigma_j\,\xi_{ij}]}_{=\,0} \;=\; 1$$

So both ensembles look **identical "one bond at a time"**, but the correlated one has hidden multi-bond structure.

**Defaults.** $\kappa = 0.7$, $\pi = 0.9$. Backbone correlation length $\approx 1/\ln(2\pi - 1) \approx 5.8$ monomers.

**Internal control.** When $\epsilon = 0$, the disorder term vanishes from the Hamiltonian, so the dynamics **must be identical** regardless of which ensemble was used to generate $\eta$. Any difference at $\epsilon > 0$ is therefore real physics, not numerics. This is our null test.

### Visual comparison — the two ensembles

```
   (A) i.i.d. ENSEMBLE  — each ηᵢⱼ independent

         ●──●──●──●──●──●──●
         1  2  3  4  5  6  7

         η₁,₃ = rand    η₁,₄ = rand
         η₂,₅ = rand    η₃,₆ = rand
         η₄,₇ = rand    ...   (no link to backbone)


   (B) CORRELATED ENSEMBLE  — σ-sequence biases ηᵢⱼ

         σ:  +  −  −  +  +  −  +
             │  │  │  │  │  │  │
             ●──●──●──●──●──●──●
             1  2  3  4  5  6  7

         η₁,₄ ∝ σ₁σ₄ = (+)(+) = +1  →  attract
         η₂,₅ ∝ σ₂σ₅ = (−)(+) = −1  →  repel
         η₃,₆ ∝ σ₃σ₆ = (−)(−) = +1  →  attract
         + independent noise of variance (1 − κ²)

   Both ensembles: identical mean, identical variance.
   Only the correlation pattern differs.
```

<div class="pgbreak"></div>

## 5. Aging protocol

1. **Initialize** the chain as a 3D random walk with bond length $r_{\min}$.
2. **Pre-equilibrate** for 2,000 MC sweeps at $\beta_0 = 0.05$ ($T = 20$). The chain is hot and explores freely.
3. **Quench** at $t = 0$: instantly drop the inverse temperature to $\beta = 1.0$ ($T = 1$). Chain is now far out of equilibrium.
4. **Wait** for $t_w \in \{0, 100, 300, 1000, 3000, 10000\}$ MC sweeps.
5. **Measure** how the chain evolves over a "lag time" $t$, comparing configurations at times $t_w$ and $t_w + t$.

An aging system depends on **two times**, not one:

- $t_w$ — the **age** (time since the quench)
- $t$ — the **lag** (time elapsed during the measurement)

In a system at equilibrium, only $t$ matters. In an aging system, **both matter independently** — a 10-second measurement looks different on a 1-minute-old sample than on a 1-hour-old sample.

**Dynamics: single-bead Metropolis Monte Carlo.** Each MC sweep consists of $N$ attempted moves. In each move:

1. Pick a bead $i$ uniformly at random.
2. Propose a displacement $\delta \in [-0.25, 0.25]^3$.
3. Compute the energy change $\Delta H$.
4. Accept with probability $\min(1, e^{-\beta\,\Delta H})$.

After each sweep we reset the center of mass to the origin so the polymer doesn't drift away.

---

## 6. Observables

### 6a. Contact overlap $Q(t_w, t)$ — structural memory

$$Q(t_w, t) \;=\; \frac{1}{N_{\text{pairs}}} \sum_{|i-j|>1} \mathbf{1}[r_{ij}(t_w) < r_c]\;\mathbf{1}[r_{ij}(t_w + t) < r_c]$$

The fraction of contact pairs that exist at *both* times $t_w$ and $t_w + t$. The normalized version $\widetilde{Q}(t_w, t) = Q(t_w, t)/Q(t_w, 0)$ starts at 1 and decays toward 0.

Aging signature: **older $t_w$ → slower decay of $\widetilde{Q}$**.

### 6b. Four-point susceptibility $\chi_4(t_w, t)$ — dynamical heterogeneity

$$\chi_4(t_w, t) \;=\; N_{\text{pairs}}\left(\langle Q^2\rangle_{\text{traj}} \;-\; \langle Q\rangle_{\text{traj}}^2\right)$$

Run many independent trajectories from the **same** aged disorder configuration. If they all relax similarly, the variance is small; if they take wildly different paths, the variance is large. **$\chi_4$ measures between-trajectory disagreement.**

Peak height: $\chi_4^\star(t_w) = \max_t \chi_4(t_w, t)$. Peak location: $t^\star(t_w) = \arg\max_t \chi_4(t_w, t)$.

<div class="pgbreak"></div>

### 6c. Mean-square displacement (MSD)

$$\langle\Delta r^2\rangle \;=\; \frac{1}{N}\sum_i \|\mathbf{r}_i(t_w + t) - \mathbf{r}_i(t_w)\|^2$$

Normal liquid: $\langle\Delta r^2\rangle \sim t$. Glass: $\langle\Delta r^2\rangle \sim t^\alpha$ with $\alpha < 1$. We find $\alpha \approx 0.64$–$0.68$.

### 6d. Non-Gaussianity parameter $\alpha_2$

$$\alpha_2(t_w, t) \;=\; \frac{3}{5}\frac{\langle r^4\rangle}{\langle r^2\rangle^2} \;-\; 1$$

Gaussian distribution: $\alpha_2 = 0$. Positive values mean **heavy tails** — most beads barely move, but a few jump a lot. Intermittent hopping is a hallmark of glassy dynamics.

### 6e. Structural displacement $D_4$

$$D_4(t_w, t) \;=\; \left\langle\bigl(d_{ij}^2(t_w + t) - d_{ij}^2(t_w)\bigr)^2\right\rangle_{i<j}$$

Mean-squared change in pairwise squared distances. Sensitive to *all* rearrangements, not just contact-breaking ones.

### 6f. KWW stretched-exponential fit

$$\widetilde{Q}(t_w, t) \;\approx\; \exp\!\left[-\left(\frac{t}{\tau_{\text{KWW}}}\right)^{\beta_{\text{KWW}}}\right]$$

| Stretching exponent | Meaning |
|---|---|
| $\beta_{\text{KWW}} = 1$ | Simple exponential — one timescale |
| $\beta_{\text{KWW}} < 1$ | **Stretched exponential** — broad spectrum of timescales |

We find $\beta_{\text{KWW}} \approx 0.30$–$0.41$.

### 6g. Partial relaxation time $\tau_c(t_w)$

The first lag $t$ at which $\widetilde{Q}(t_w, t)$ crosses a threshold $c$. For $c = e^{-1} \approx 0.368$, this is the standard "$\alpha$-relaxation time." We also track $c \in \{0.6, 0.8\}$. Right-censored if never reached within $t_{\max} = 10^4$.

### 6h. Aging exponent $\mu$

Fit $\tau_c(t_w) \sim t_w^\mu$ in log–log space. $\mu = 1$ is simple aging; $\mu < 1$ is subaging; $\mu > 1$ is superaging.

<div class="pgbreak"></div>

## 7. Findings

**1. At $\epsilon = 0$, both ensembles agree exactly.** Internal control passes. Disorder term vanishes, dynamics are identical.

**2. At $\epsilon > 0$, the ensembles diverge.** Correlated $\widetilde{Q}$ decays **faster** than i.i.d. $\widetilde{Q}$ at the same $(\epsilon, t_w)$. Counterintuitive — structured couplings might be expected to create more rigidity.

**3. i.i.d. disorder enhances late-age dynamical heterogeneity.** At $t_w = 10^4$:
- $\epsilon = 3$: $\chi_4^\star \approx 0.283$ (i.i.d.) vs. $\approx 0.206$ (correlated) → **38% enhancement**
- $\epsilon = 6$: $\chi_4^\star \approx 0.295$ (i.i.d.) vs. $\approx 0.205$ (correlated) → **44% enhancement**

**4. $\chi_4^\star(t_w)$ is non-monotonic in age.** Peak heterogeneity rises from early ages, peaks around $t_w \sim 10^3$, then *decreases* by $t_w = 10^4$. At intermediate ages the system is exploring; at late ages it gets locked into a narrower set of pathways and trajectories agree more.

**5. The heterogeneity timescale $t^\star$ is delayed under i.i.d. disorder.** At $\epsilon = 6$, $t_w = 10^4$: $t^\star \approx 2555$ (i.i.d.) vs. $t^\star \approx 918$ (correlated) — **nearly threefold**.

### Visual comparison — the headline result

```
   PEAK HETEROGENEITY  χ₄*  AT LATE AGE  (t_w = 10⁴)

                    0.00    0.10    0.20    0.30
                    ├───────┼───────┼───────┤

   ε = 3  i.i.d.    ▉▉▉▉▉▉▉▉▉▉▉▉▉▉ 0.283
          correl.   ▉▉▉▉▉▉▉▉▉▉     0.206   → −38%

   ε = 6  i.i.d.    ▉▉▉▉▉▉▉▉▉▉▉▉▉▉▉ 0.295
          correl.   ▉▉▉▉▉▉▉▉▉▉     0.205   → −44%


   PEAK LOCATION  t*   at  ε = 6,  t_w = 10⁴:

          correlated:  t* ≈   918    (earlier)
          i.i.d.    :  t* ≈  2555    (~3× delayed)


   KWW STRETCHING  β   at  ε = 6,  t_w = 10⁴:

          correlated:  β ≈ 0.30   (broader spectrum)
          i.i.d.    :  β ≈ 0.33   (narrower spectrum)


   NON-GAUSSIANITY PEAK  α₂*   at  ε = 6,  t_w = 10⁴:

          correlated:  α₂* ≈ 8.0
          i.i.d.    :  α₂* ≈ 9.8   → +22%


   SUBDIFFUSION EXPONENT  α:

          correlated:  α ≈ 0.64   (more trapped locally)
          i.i.d.    :  α ≈ 0.66   (slightly freer hops)
```

**6. Strong glassiness pushes $e^{-1}$ out of the window.** For $\epsilon \geq 3$ and $t_w \geq 100$, $\widetilde{Q}$ never reaches $e^{-1}$ within $t_{\max} = 10^4$. The system has effectively **frozen on accessible timescales**.

**7. $\kappa$ continuously tunes heterogeneity.** Sweeping $\kappa$ from 0 (i.i.d.) to 0.9 (strongly correlated) at $t_w = 10^4$ smoothly suppresses $\chi_4^\star$ by ~40%. At intermediate ages ($t_w = 10^3$) the response is **non-monotonic** — moderate correlations can *increase* $\chi_4^\star$, strong correlations suppress it. **This rules out a simple "single effective disorder strength" interpretation.**

**8. KWW stretching becomes broader with correlations.** $\beta_{\text{KWW}} \approx 0.30$ (correlated) vs. $0.33$ (i.i.d.) at $\epsilon = 6$, $t_w = 10^4$. **Correlated disorder produces a broader spectrum of relaxation times**, even though the mean relaxation is faster.

**9. Non-Gaussianity confirms ensemble-dependent intermittency.** $\alpha_2^\star \approx 9.8$ (i.i.d.) vs. $\approx 8.0$ (correlated) at $\epsilon = 6$, $t_w = 10^4$ — a ~22% enhancement that mirrors the $\chi_4^\star$ enhancement.

**10. Subdiffusion is slightly more constrained under correlated disorder.** $\alpha \approx 0.64$ (correlated) vs. $0.66$ (i.i.d.). **Correlated disorder traps beads more locally, even as it speeds contact-network turnover.** A real tension: locally more constrained, contact-pattern-wise more flexible.

<div class="pgbreak"></div>

## 8. Interpretation

In the correlated ensemble, $\eta_{ij} = \kappa\,\sigma_i\sigma_j + \text{noise}$ injects **coherence** — it tells groups of monomers what to do *together*. This biases the polymer toward a smaller, more compatible set of preferred contact patterns.

A single mechanism explains all four key dynamical effects:

| Effect | Mechanism |
|---|---|
| Faster contact decorrelation | With fewer "good" patterns to choose from, the system can shuffle between them faster. |
| Less dynamical heterogeneity | Trajectories from the same aged config more reliably end up in the same coherent pattern, so they diverge less. |
| More quenched (sample-to-sample) variability | Different disorder realizations now look more different (different $\sigma$ sequences → different preferred patterns), even if individual realizations are dynamically more uniform. |
| More constrained beads (lower $\alpha$) | Each bead is locally trapped by its sequence-imposed preferences, so MSD subdiffusion is slightly slower. |

Correlations don't just renormalize disorder strength. They **reshape the topology of the energy landscape** — reducing the number of mutually-incompatible deep basins (which is what drives dynamical heterogeneity in i.i.d. disorder) but *sharpening* each individual basin.

### Visual comparison — the energy landscape

```
   ENERGY-LANDSCAPE CARTOON  (a 1-D schematic of a
   many-dimensional rugged surface)


   (A) i.i.d. — many shallow competing basins

        ╲_╱╲_╱╲_╱╲_╱╲_╱╲_╱

        • trajectories scatter into many basins
          → HIGH χ₄  (big trajectory disagreement)
        • intermittent hopping is easy
          → HIGH α₂  (more non-Gaussian)
        • beads leak basin-to-basin
          → higher α  (α ≈ 0.66,  slightly freer)


   (B) Correlated — fewer, deeper, sharper basins

        ╲________╱╲________╱╲________╱

        • trajectories funnel into a few basins
          → LOW χ₄  (trajectories agree more)
        • intermittent hopping suppressed
          → LOW α₂  (more Gaussian)
        • beads trapped deep in basin
          → lower α  (α ≈ 0.64,  more caged)


   KEY INSIGHT
   Correlations REDUCE the count of competing basins
   but DEEPEN each remaining basin.
   That is why dynamical heterogeneity DROPS while
   local subdiffusion TIGHTENS — one statement,
   two seemingly opposite effects.
```

---

## 9. Implications

**For glass physics.** Identifies a previously-overlooked control parameter — correlation topology — that's distinct from the usual one (disorder variance). A new axis in the theory of glassy materials.

**For protein / IDP biology.** Real intrinsically disordered proteins have sequence correlations (charged blocks, hydrophobic patches). This work suggests these correlations qualitatively shape their aging dynamics, not just the strength of their interactions.

**For sequence design.** A synthetic heteropolymer with target glassy properties can be tuned via $\kappa$ (and $\pi$) to dial in heterogeneity, timescale spread, and the dynamical-vs-quenched variability balance — **without changing the amount of disorder**.

<div class="pgbreak"></div>

## 10. Next steps — and where lab partnerships come in

Everything above is purely computational. To turn this into a major result rather than an isolated simulation paper, the predictions need experimental confrontation and the computational reach has to expand. Three classes of next step, roughly ordered by how dependent they are on lab infrastructure.

### A. Experimental confrontation

The cleanest test of the central claim — *correlation structure matters at fixed disorder amplitude* — needs a system where sequence statistics can be controlled and aging dynamics measured directly.

**Single-molecule FRET on intrinsically disordered proteins.** smFRET tracks inter-residue distances on individual molecules over time. After a thermal or chemical quench (temperature jump, denaturant dilution), per-molecule trajectories give a direct experimental analog of $\chi_4$: how much do trajectories from the same starting condition disagree? Pair this with **sequence-shuffled IDP variants** — same amino-acid composition, scrambled order — and you have a near-perfect experimental version of the i.i.d. vs. correlated comparison. **Needs:** a lab doing high-throughput smFRET on IDPs.

**Aging in liquid–liquid phase-separated condensates.** Many IDPs form condensates that themselves age (liquid → gel → glass over hours). Sequence statistics strongly affect condensate properties. Test: do block-copolymer-like sequences (high $\sigma$-correlation) and random-permuted sequences with identical composition show different aging timescales and KWW stretching exponents? **Needs:** a phase-separation / condensate biophysics lab with FRAP or microrheology capability.

**Sequence-defined synthetic polymers.** Modern polymer chemistry (Lutz-style sequence-controlled polymerization, DNA-coated patchy colloids) can produce polymers with statistically designed monomer sequences. Dielectric spectroscopy or rheology on a controlled series — same composition, varied $\sigma$-correlation length — would directly probe ensemble-dependent $\beta_{\text{KWW}}$. **Needs:** a synthesis lab plus low-frequency mechanical/dielectric spectroscopy.

### B. Scaling up the simulations

Current study: $N = 30$ monomers, $t_{\max} = 10^4$ sweeps. Real disordered proteins are 100–300 residues, and the most interesting aging effects (subaging, asymptotic exponents) only emerge at much longer times.

- **GPU port of the Metropolis kernel** — single-bead moves are embarrassingly parallel within a sweep with proper indexing. Estimated 50–100× speedup, opening $N = 100$–$300$ at the current wall time.
- **Replica-exchange variants** to accelerate the late-age regime where moves at $\beta = 1$ are almost always rejected.
- **Inherent-structure analysis** (gradient-descent quenches at each $t_w$) to directly map the energy-landscape descent the paper currently only infers.

**Needs:** HPC allocation or a lab GPU cluster. Modest by physics-lab standards (~50k GPU-hours for the full $N = 100$ scan).

### C. Real-data reanalysis

There are existing datasets that may already contain ensemble-dependent aging signatures nobody has looked for.

- **NMR relaxation dispersion** on IDPs encodes multi-timescale dynamics. Reanalyze existing measurements for stretched-exponential $\beta$ correlated with sequence $\sigma$-correlation.
- **All-atom MD trajectory archives** (D.E. Shaw IDP runs, Folding@Home, Anton trajectories) — extract $Q(t_w, t)$ and $\chi_4$ directly and check the $\sigma$-correlation prediction on real biological sequences.
- **Bioinformatic predictor.** Given an arbitrary IDP sequence, compute its $\sigma$-correlation length from a chosen flavor decomposition (charge, hydrophobicity), and predict relative aging timescale. Compare against measured condensate aging rates.

**Needs:** collaborators with NMR datasets, MD trajectory access, or proteomics databases.

### Why lab affiliation matters

A lab affiliation gives this work three things it can't get on its own: **(1)** experimental partners who can falsify or confirm the predictions on real biology; **(2)** computational scale beyond a single workstation; **(3)** the credibility and platform to argue that the i.i.d. assumption deserves to be retired.

Strongest natural homes:

- **Biophysics groups** working on IDPs or biomolecular condensates.
- **Glass-physics theory groups** extending random-heteropolymer / replica theory.
- **Sequence-controlled polymer chemistry labs** that can synthesize and characterize the predicted designs.

Even a brief hosted visit — running the smFRET-on-shuffled-IDP experiment, or applying the $\sigma$-correlation analysis to an existing all-atom trajectory database — would convert this from a methodology paper into a result with biological and material-design reach.

---

## Appendix A — Parameter cheat-sheet

| Symbol | Name | Default |
|---|---|---|
| $N$ | Chain length | 30 |
| $A$ | LJ attraction coefficient | 3.8 |
| $R$ | LJ repulsion coefficient | 2.0 |
| $h$ | Harmonic bond stiffness | 1.0 |
| $\beta_0$ | Pre-quench inverse temp | 0.05 ($T = 20$) |
| $\beta$ | Post-quench inverse temp | 1.0 ($T = 1$) |
| $\epsilon$ | Disorder amplitude | $\{0, 3, 6\}$ |
| $\kappa$ | Correlation strength | 0.7; scanned $\in \{0, 0.2, 0.4, 0.6, 0.8, 0.9\}$ |
| $\pi$ | Backbone persistence | 0.9 |
| $\delta$ | MC step size | 0.25 |
| $r_c$ | Contact cutoff | $\approx 1.259$ |
| $t_w$ | Waiting times | $\{0, 100, 300, 1000, 3000, 10000\}$ |
| $t_{\max}$ | Max lag | $10^4$ |
| $n_{\text{disorder}}$ | Disorder realizations | 20 (main), 12 (extended) |
| $n_{\text{traj}}$ | Trajectories per realization | 8 (main), 6 (extended) |

---

## Appendix B — Mental model

1. The polymer lives on a **rugged energy landscape** because random pair couplings can't all be satisfied (frustration).
2. After a temperature **quench**, the polymer slowly explores this landscape — relaxing more slowly the older it gets (**aging**).
3. The **structure of the disorder** (not just its variance) determines how rugged the landscape is and how the trajectories spread out (**dynamical heterogeneity**).
