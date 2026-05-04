---
title: "IMP Heteropolymer Aging — Working Notes"
tags: [glass-physics, heteropolymer, aging, dynamical-heterogeneity, IMP-model, monte-carlo]
date: 2026-04-20
---

# IMP Heteropolymer Aging — Working Notes

> [!abstract] Thesis
> In a frustrated heteropolymer, ==**how the random couplings are correlated along the backbone matters as much as how strong they are**==. This opens a new axis of glassy-aging behavior even when the disorder amplitude $\epsilon$ is held fixed.

---

## Contents

1. [[#0. Glossary]]
2. [[#1. Why this matters]]
3. [[#2. The IMP heteropolymer]]
4. [[#3. Diagram]]
5. [[#4. The two disorder ensembles]]
6. [[#5. Aging protocol]]
7. [[#6. Observables]]
8. [[#7. Findings]]
9. [[#8. Interpretation]]
10. [[#9. Implications]]

---

## 0. Glossary

> [!example] ==**Quench**==
> A sudden drop in temperature. In simulations, instantly change the inverse temperature $\beta = 1/(k_B T)$ from a small value (hot) to a large value (cold). After a quench, the system is *out of equilibrium*: it was happy at the old temperature and now has to slowly settle into the new one. Everything interesting in this work happens after the quench.

> [!example] ==**Aging**==
> The dynamics depend on how long it has been since the quench, not just on the current time. A 1-second-old glass relaxes faster than a 1-hour-old glass. The "age" is the **waiting time** $t_w$. Aging breaks time-translation invariance.

> [!example] ==**Heteropolymer**==
> A polymer (chain of beads) where the beads are not all the same. Different beads have different chemistries and therefore different interactions. Real-world examples: proteins, intrinsically disordered proteins.

> [!example] ==**Frustration**==
> No single configuration can simultaneously satisfy every local preference. Classic example: three Ising spins on a triangle each wanting to anti-align. In a heteropolymer, no 3D shape can make every "wants-to-stick" pair stick *and* every "wants-to-repel" pair stay apart.

> [!example] ==**Quenched disorder**==
> The random coupling $\eta_{ij}$ is drawn once when we build the polymer and never changes. Contrast with *annealed* disorder, where the random variables also fluctuate. Quenched disorder is the standard model for sequence-level randomness in proteins.

> [!example] ==**Glassy dynamics**==
> Slow, history-dependent, multi-timescale relaxation observed in window glass, supercooled liquids, gels, and frustrated polymers. Hallmarks: aging, subdiffusion, dynamical heterogeneity, broken time-translation invariance, stretched-exponential correlations.

> [!example] ==**Dynamical heterogeneity**==
> Different parts of the system relax at very different rates. Some regions are temporarily frozen while others are mobile. Trajectories from the same starting configuration disagree about how the system evolves. Quantified by $\chi_4$.

> [!example] ==**Subdiffusion**==
> Motion slower than normal diffusion. Normal: $\langle \Delta r^2 \rangle \sim t$. Subdiffusion: $\langle \Delta r^2 \rangle \sim t^\alpha$ with $\alpha < 1$. A signature of trapping or caging.

> [!example] ==**Ensemble**==
> A recipe for drawing the random couplings $\eta_{ij}$. We compare two: i.i.d. Gaussian, and a sequence-structured correlated ensemble.

> [!example] ==**Marginal variance**==
> The variance of a single random variable, ignoring its relationship to the others. Both ensembles here have the same marginal variance $\mathbb{E}[\eta_{ij}^2] = 1$, so any difference must come from correlations *between* different $\eta_{ij}$'s.

> [!example] ==**Lennard–Jones potential**==
> $V(r) = R/r^{12} - A/r^6$. The $r^{-12}$ piece is short-range repulsion, the $-r^{-6}$ piece is long-range attraction. Minimum at $r_{\min} = (2R/A)^{1/6}$.

> [!example] ==**Metropolis Monte Carlo**==
> Algorithm for sampling configurations at temperature $T$. Propose a small random change, compute $\Delta H$, accept with probability $\min(1, e^{-\beta \Delta H})$.

> [!example] ==**KWW (Kohlrausch–Williams–Watts) function**==
> Stretched exponential $f(t) = \exp[-(t/\tau)^\beta]$. $\beta = 1$ is normal exponential. $\beta < 1$ signals a broad distribution of relaxation timescales — fast and slow modes mixed.

> [!example] ==**Four-point susceptibility $\chi_4$**==
> Statistical quantity measuring how much trajectories disagree with each other. Its peak height $\chi_4^\star$ is the standard glass-physics diagnostic for dynamical heterogeneity.

> [!example] ==**Inverse temperature $\beta$**==
> $\beta = 1/(k_B T)$. Large $\beta$ = cold; small $\beta$ = hot. Used because it makes Boltzmann factors $e^{-\beta H}$ cleaner.

> [!example] ==**MC sweep**==
> One unit of simulation time, defined as $N$ attempted single-bead moves. Each bead gets, on average, one move attempt per sweep.

---

## 1. Why this matters

Glasses, gels, and frustrated polymers all show ==**aging**==: their dynamics slow down the longer you wait after preparing them. The standard minimal model in this space is the IMP heteropolymer, which puts a random attractive coupling between every pair of monomers. Everyone uses the simplest random ensemble — independent Gaussians, no correlations between pairs. Real biological sequences have correlations along the backbone. Does that matter for the dynamics?

We find: ==**yes, dramatically**==. The correlation structure of the disorder is an independent control knob, comparable in importance to the disorder amplitude.

Three things to keep in mind throughout:

1. Glasses age. A freshly-quenched glass relaxes fast at first; the older it gets, the slower it relaxes.
2. Heteropolymers are minimal models for glasses — a chain of beads with mismatched, frustrated interactions captures the essential ingredients.
3. The random couplings $\eta_{ij}$ have a structure we choose. Standard practice has been i.i.d. We push back on that.

> [!info] Two meanings of "correlation" — keep them separate
> - **Input correlation** — built into the random couplings $\eta_{ij}$. The i.i.d. ensemble has none; the correlated ensemble injects some via the backbone field $\sigma_i \sigma_j$. ==**This is the new control knob we tune.**==
> - **Output correlation** — what $\chi_4$ diagnoses. It correlates the overlap $Q$ across trajectories to measure how much they *disagree* — i.e., dynamical heterogeneity. This is what we *observe*.
>
> The whole study is: turn the input knob (1) and watch what happens to the output (2). Plus a side motivation — real protein sequences have input-correlation, so the i.i.d. default that 30+ years of literature has used is a convenience nobody had tested.

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

The first three terms together are a clean Lennard-Jones polymer. The fourth term ==**adds the disorder**==: some pairs become extra-attractive ($\eta_{ij} > 0$), others extra-repulsive ($\eta_{ij} < 0$). Because the chain is a connected object in 3D, ==**it can't satisfy all of these preferences at once**==. That's the **frustration**, and it's what makes the model glassy.

Parameters used:

- $A = 3.8$, $R = 2.0$, $h = 1.0$, $N = 30$
- LJ minimum: $r_{\min} = (2R/A)^{1/6} \approx 1.007$
- Contact cutoff: $r_c = 1.25\,r_{\min} \approx 1.259$

Two non-bonded beads ($|i-j|>1$) count as "in contact" when $r_{ij} < r_c$.

---

## 3. Diagram

```
                   Backbone sequence: σ = (+ − − + + − + + − +)

                                          ╭ ─ ─ ─ ─ ─ ─ ─ ─ ╮
                                          │  η₃,₈ (random)  │
                                          │                 ▼
    (+)─────(−)─────(−)─────(+)─────(+)─────(−)─────(+)─────(+)─────(−)─────(+)
     1   ↑   2       3   ↑   4       5       6   ↑   7       8       9      10
         │               │                       │
         │               ╰ ─ ─ ─ ╮               │
         │   η₂,₅          η₄,₉  │               │
         │                       ▼               │
         ╰ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ╯
                                                  η₁,₇

    ───  Harmonic backbone bond:  h · r²ᵢ,ᵢ₊₁         (rigid chain spine)
    ─ ─  Quenched LJ + disorder:  R/r¹² − A/r⁶ + √ε · ηᵢⱼ / r⁶
         (acts on EVERY non-bonded pair; only a few drawn here)

                      (+) and (−) denote σᵢ ∈ {+1, −1}
```

Solid lines are the backbone — every adjacent pair is glued together by a harmonic spring. Dashed arrows show a few of the many non-bonded couplings $\eta_{ij}$. In reality, every non-adjacent pair has one. For $N = 10$ beads, that's $\binom{10}{2} - 9 = 36$ random couplings; for $N = 30$, it's $\binom{30}{2} - 29 = 406$.

The $(+)$ and $(-)$ are the **backbone flavors** $\sigma_i$. In the i.i.d. ensemble these are irrelevant — every $\eta_{ij}$ is its own independent Gaussian. In the **correlated ensemble**, like-flavor pairs ($\sigma_i \sigma_j = +1$) get a positive bias and unlike-flavor pairs get a negative bias.

After the quench, the chain crumples into a globule. The random $\eta_{ij}$'s decide which non-bonded contacts get rewarded and which get punished. ==**Frustration → rugged energy landscape → glassy aging.**==

---

## 4. The two disorder ensembles

Both ensembles have ==**identical mean (0) and identical variance (1)**==. The only difference is in the correlation structure between *different* $\eta_{ij}$'s.

### (A) i.i.d. ensemble (standard)

$$\eta_{ij} \sim \mathcal{N}(0,1) \quad \text{independently for each pair } (i,j)$$

Every pair gets its own independent draw from a standard Gaussian.

### (B) Correlated ensemble (the new construction)

First, generate a backbone "flavor" sequence $\sigma_i \in \{+1, -1\}$ from a two-state Markov chain with persistence $\pi$:

$$P(\sigma_{i+1} = \sigma_i) = \pi$$

The autocorrelation along the backbone decays geometrically:

$$\langle \sigma_i \sigma_{i+k}\rangle = (2\pi - 1)^{|k|}$$

Then construct the pair couplings as a weighted mix of structured and noise:

$$\eta_{ij} \;=\; \kappa\,\sigma_i\sigma_j \;+\; \sqrt{1 - \kappa^2}\,\xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0, 1)$$

where $\kappa \in [0, 1]$ controls how much of $\eta_{ij}$ comes from the structured part vs. noise.

The $\sqrt{1-\kappa^2}$ weighting preserves the marginal variance:

$$\mathbb{E}[\eta_{ij}^2] \;=\; \kappa^2 \underbrace{\mathbb{E}[(\sigma_i\sigma_j)^2]}_{=\,1} + (1-\kappa^2)\underbrace{\mathbb{E}[\xi_{ij}^2]}_{=\,1} + 2\kappa\sqrt{1-\kappa^2}\underbrace{\mathbb{E}[\sigma_i\sigma_j\,\xi_{ij}]}_{=\,0} \;=\; 1$$

So both ensembles look ==**identical "one bond at a time"**==, but the correlated one has hidden multi-bond structure.

Defaults: $\kappa = 0.7$, $\pi = 0.9$. Backbone correlation length $\approx 1/\ln(2\pi - 1) \approx 5.8$ monomers.

> [!success] Internal control
> When $\epsilon = 0$, the disorder term vanishes from the Hamiltonian, so the dynamics ==**must be identical**== regardless of which ensemble was used to generate $\eta$. Any difference at $\epsilon > 0$ is therefore real physics, not numerics. This is our null test.

---

## 5. Aging protocol

1. **Initialize** the chain as a 3D random walk with bond length $r_{\min}$.
2. **Pre-equilibrate** for 2,000 MC sweeps at high temperature: $\beta_0 = 0.05$ ($T = 20$). The chain is hot and explores freely.
3. ==**Quench**== at $t = 0$: instantly drop the inverse temperature to $\beta = 1.0$ ($T = 1$). The chain is now far out of equilibrium.
4. **Wait** for $t_w \in \{0, 100, 300, 1000, 3000, 10000\}$ MC sweeps.
5. **Measure** how the chain evolves over a "lag time" $t$, comparing configurations at times $t_w$ and $t_w + t$.

An aging system depends on **two times**, not one:

- $t_w$ — the **age** (time since the quench)
- $t$ — the **lag** (time elapsed during the measurement)

In a system at equilibrium, only $t$ matters. In an aging system, ==**both matter independently**== — a 10-second measurement looks different on a 1-minute-old sample than on a 1-hour-old sample.

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

The fraction of contact pairs that exist at *both* times $t_w$ and $t_w + t$. The normalized version $\widetilde{Q}(t_w, t) = Q(t_w, t)/Q(t_w, 0)$ starts at 1 and decays toward 0 as the chain forgets its old contacts.

Aging signature: ==**older $t_w$ → slower decay of $\widetilde{Q}$**==. The system "remembers" its previous configuration longer the older it is.

### 6b. Four-point susceptibility $\chi_4(t_w, t)$ — dynamical heterogeneity

$$\chi_4(t_w, t) \;=\; N_{\text{pairs}}\left(\langle Q^2\rangle_{\text{traj}} \;-\; \langle Q\rangle_{\text{traj}}^2\right)$$

Run many independent trajectories from the **same** aged disorder configuration. If they all relax similarly, the variance is small; if they take wildly different paths, the variance is large. ==**$\chi_4$ measures *between-trajectory disagreement* — the canonical signature of dynamical heterogeneity.**==

$\chi_4$ typically shows a peak at a characteristic lag $t^\star$ where trajectories disagree maximally:

- Peak height: $\chi_4^\star(t_w) = \max_t \chi_4(t_w, t)$
- Peak location: $t^\star(t_w) = \arg\max_t \chi_4(t_w, t)$

### 6c. Mean-square displacement (MSD)

$$\langle\Delta r^2\rangle \;=\; \frac{1}{N}\sum_i \|\mathbf{r}_i(t_w + t) - \mathbf{r}_i(t_w)\|^2$$

How far does a typical bead move?

- Normal liquid: $\langle\Delta r^2\rangle \sim t$ (linear → diffusive)
- Glass / aging system: $\langle\Delta r^2\rangle \sim t^\alpha$ with $\alpha < 1$ (subdiffusive)

We find $\alpha \approx 0.64$–$0.68$.

### 6d. Non-Gaussianity parameter $\alpha_2$

$$\alpha_2(t_w, t) \;=\; \frac{3}{5}\frac{\langle r^4\rangle}{\langle r^2\rangle^2} \;-\; 1$$

For a Gaussian displacement distribution, $\langle r^4\rangle / \langle r^2\rangle^2 = 5/3$ exactly, so $\alpha_2 = 0$. Positive values mean the distribution has ==**heavy tails**==: most beads barely move, but a few jump a lot. This intermittent "hopping" is a hallmark of glassy dynamics.

### 6e. Structural displacement $D_4$

$$D_4(t_w, t) \;=\; \left\langle\bigl(d_{ij}^2(t_w + t) - d_{ij}^2(t_w)\bigr)^2\right\rangle_{i<j}$$

Mean-squared change in pairwise squared distances. Sensitive to *all* rearrangements, not just contact-breaking ones.

### 6f. KWW stretched-exponential fit

$$\widetilde{Q}(t_w, t) \;\approx\; \exp\!\left[-\left(\frac{t}{\tau_{\text{KWW}}}\right)^{\beta_{\text{KWW}}}\right]$$

| Stretching exponent | Meaning |
|---|---|
| $\beta_{\text{KWW}} = 1$ | Simple exponential — one timescale |
| $\beta_{\text{KWW}} < 1$ | ==**Stretched exponential**== — broad spectrum of timescales |

We find $\beta_{\text{KWW}} \approx 0.30$–$0.41$ — strongly stretched, glass-like.

### 6g. Partial relaxation time $\tau_c(t_w)$

The first lag $t$ at which $\widetilde{Q}(t_w, t)$ crosses a threshold $c$. For $c = e^{-1} \approx 0.368$, this is the standard "$\alpha$-relaxation time." We also track $c \in \{0.6, 0.8\}$. If the threshold is never reached within $t_{\max} = 10^4$, the value is right-censored.

### 6h. Aging exponent $\mu$

Fit $\tau_c(t_w) \sim t_w^\mu$ in log–log space.

| $\mu$ | Regime |
|---|---|
| $\mu = 1$ | Simple aging — relaxation time grows linearly with age |
| $\mu < 1$ | Subaging |
| $\mu > 1$ | Superaging |

---

## 7. Findings

> [!success] 1 — At $\epsilon = 0$, both ensembles agree exactly
> Internal control passes. Disorder term vanishes, dynamics are identical. ✓

> [!success] 2 — At $\epsilon > 0$, the ensembles diverge
> Correlated $\widetilde{Q}$ decays ==**faster**== than i.i.d. $\widetilde{Q}$ at the same $(\epsilon, t_w)$. Counterintuitive — structured couplings might be expected to create more rigidity, not less.

> [!success] 3 — i.i.d. disorder enhances late-age dynamical heterogeneity
> At $t_w = 10^4$:
> - $\epsilon = 3$: $\chi_4^\star \approx 0.283$ (i.i.d.) vs. $\approx 0.206$ (correlated) → ==**38% enhancement**==
> - $\epsilon = 6$: $\chi_4^\star \approx 0.295$ (i.i.d.) vs. $\approx 0.205$ (correlated) → ==**44% enhancement**==

> [!success] 4 — $\chi_4^\star(t_w)$ is non-monotonic in age
> Peak heterogeneity rises from early ages, peaks around $t_w \sim 10^3$, then *decreases* by $t_w = 10^4$. Reading: at intermediate ages the system is exploring; at late ages it gets locked into a narrower set of pathways and trajectories agree more.

> [!success] 5 — The heterogeneity timescale $t^\star$ is delayed under i.i.d. disorder
> At $\epsilon = 6$, $t_w = 10^4$: $t^\star \approx 2555$ (i.i.d.) vs. $t^\star \approx 918$ (correlated) — ==**nearly threefold**==.

> [!warning] 6 — Strong glassiness pushes $e^{-1}$ out of the window
> For $\epsilon \geq 3$ and $t_w \geq 100$, $\widetilde{Q}$ never reaches $e^{-1}$ within $t_{\max} = 10^4$. The system has effectively **frozen on accessible timescales** — the operational signature of glassiness.

> [!success] 7 — $\kappa$ continuously tunes heterogeneity
> Sweeping $\kappa$ from 0 (i.i.d.) to 0.9 (strongly correlated) at $t_w = 10^4$ smoothly suppresses $\chi_4^\star$ by ~40%. At intermediate ages ($t_w = 10^3$) the response is **non-monotonic** — moderate correlations can *increase* $\chi_4^\star$, strong correlations suppress it. ==This rules out a simple "single effective disorder strength" interpretation.==

> [!success] 8 — KWW stretching becomes broader with correlations
> $\beta_{\text{KWW}} \approx 0.30$ (correlated) vs. $0.33$ (i.i.d.) at $\epsilon = 6$, $t_w = 10^4$. ==Correlated disorder produces a broader spectrum of relaxation times==, even though the mean relaxation is faster.

> [!success] 9 — Non-Gaussianity confirms ensemble-dependent intermittency
> $\alpha_2^\star \approx 9.8$ (i.i.d.) vs. $\approx 8.0$ (correlated) at $\epsilon = 6$, $t_w = 10^4$ — a ~22% enhancement that mirrors the $\chi_4^\star$ enhancement, confirming the heterogeneity is real bead-level intermittency.

> [!success] 10 — Subdiffusion is slightly more constrained under correlated disorder
> $\alpha \approx 0.64$ (correlated) vs. $0.66$ (i.i.d.). ==Correlated disorder traps beads more locally, even as it speeds contact-network turnover.== A real tension: locally more constrained, contact-pattern-wise more flexible.

---

## 8. Interpretation

In the correlated ensemble, $\eta_{ij} = \kappa\,\sigma_i\sigma_j + \text{noise}$ injects ==**coherence**== — it tells groups of monomers what to do *together*. This biases the polymer toward a smaller, more compatible set of preferred contact patterns.

This single mechanism explains all four key dynamical effects:

| Effect | Mechanism |
|---|---|
| Faster contact decorrelation | With fewer "good" patterns to choose from, the system can shuffle between them faster. |
| Less dynamical heterogeneity | Trajectories from the same aged config more reliably end up in the same coherent pattern, so they diverge less. |
| More quenched (sample-to-sample) variability | Different disorder realizations now look more different (different $\sigma$ sequences → different preferred patterns), even if individual realizations are dynamically more uniform. |
| More constrained beads (lower $\alpha$) | Each bead is locally trapped by its sequence-imposed preferences, so MSD subdiffusion is slightly slower. |

Correlations don't just renormalize disorder strength. They ==**reshape the topology of the energy landscape**== — reducing the number of mutually-incompatible deep basins (which is what drives dynamical heterogeneity in i.i.d. disorder) but *sharpening* each individual basin.

---

## 9. Implications

**For glass physics.** Identifies a previously-overlooked control parameter — correlation topology — that's distinct from the usual one (disorder variance). A new axis in the theory of glassy materials.

**For protein / IDP biology.** Real intrinsically disordered proteins have sequence correlations (charged blocks, hydrophobic patches, etc.). This work suggests these correlations qualitatively shape their aging dynamics, not just the strength of their interactions.

**For sequence design.** A synthetic heteropolymer with target glassy properties can be tuned via $\kappa$ (and $\pi$) to dial in heterogeneity, timescale spread, and the dynamical-vs-quenched variability balance — ==**without changing the amount of disorder**==.

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
