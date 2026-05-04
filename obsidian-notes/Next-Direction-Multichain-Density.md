---
title: "Next direction — multi-chain A/B copolymer + density-field analysis"
tags: [next-steps, multichain, density, structure-factor, brownian-dynamics, feedback]
date: 2026-04-22
---

# Next direction — multi-chain A/B copolymer + density-field analysis

> [!abstract] Source
> Whiteboard meeting feedback, 2026-04-22. Three board photos archived alongside this note. The current single-chain edge-disorder model is being upgraded to a many-chain bead-type model with proper density-field observables and real-time dynamics.

---

## 1. What's wrong with the current setup (per feedback)

> [!warning] Five concrete critiques
> 1. **No density.** A single $N=30$ chain in vacuum has no density field, so we can't compute the structure factor $S(k)$ — which is the standard glass/melt observable and the thing that actually identifies length scales of segregated domains.
> 2. **Edge-disorder is the wrong knob.** Real heteropolymers have *bead types* (A, B), not pair-specific random couplings. A/B identity is a node property, and that's what should drive the physics.
> 3. **MMC is not real dynamics.** Markov-chain accept/reject introduces frustration artifacts and gives no physical timescale. Need Brownian dynamics for real-time kinetics.
> 4. **Single chain is a toy.** No melt, no condensate, no real polymer-physics phase behavior. We need many chains in a box.
> 5. **No phase-transition characterization.** The segregation/clustering of A vs. B is the actual interesting transition; we never measured it.

---

## 2. The new model — A/B multi-chain copolymer

> [!info] Switch from edge-disorder to bead-type
> Each bead $j$ on chain $i$ carries a binary type indicator (from the whiteboard):
>
> $$\sigma_{ij}^{(A)} = \begin{cases} 1 & \text{if the } j\text{-th bead on the } i\text{-th chain is type A} \\ 0 & \text{otherwise} \end{cases}$$
>
> and analogously $\sigma_{ij}^{(B)} = 1 - \sigma_{ij}^{(A)}$. The simulation box contains many such chains.

The **interactions** become type-specific (Lennard–Jones with three couplings):

$$V_{\alpha\beta}(r) = \frac{R_{\alpha\beta}}{r^{12}} - \frac{A_{\alpha\beta}}{r^{6}}, \qquad \alpha,\beta \in \{A,B\}$$

with $A_{AA}, A_{BB}, A_{AB}$ tuned so like-types attract more strongly than unlike-types — the standard Flory–Huggins-style asymmetry that drives microphase separation. The block/random/correlated character of each chain's sequence becomes the analog of our old $\kappa$ knob, but now built into the bead-type pattern along the chain rather than into pair-specific $\eta_{ij}$.

---

## 3. The new observables — density fields and structure factor

> [!info] From the whiteboard

**Real-space density of A-type beads:**

$$\phi_A(\mathbf{r}) = \nu_A \sum_i \sum_j \sigma_{ij}^{(A)} \, \delta(\mathbf{r} - \mathbf{r}_{ij})$$

where $\nu_A$ is the per-bead volume (or analogous normalization) and the sum runs over chains $i$ and beads $j$.

**Fourier-space density:**

$$\hat\phi_A(\mathbf{k}) = \nu_A \sum_i \sum_j \sigma_{ij}^{(A)} \, e^{i\mathbf{k}\cdot\mathbf{r}_{ij}}$$

**Two-point density correlation (real space):**

$$\langle \phi_A(\mathbf{r})\,\phi_A(\mathbf{r}')\rangle = C_{AA}(\mathbf{r}-\mathbf{r}')$$

This decays over a characteristic distance $\xi_{AA}$ — the **A-domain correlation length**.

**Two-point density correlation (Fourier space) = structure factor:**

$$\hat C_{AA}(\mathbf{k}) = \langle \hat\phi_A(\mathbf{k})\,\hat\phi_A(-\mathbf{k})\rangle = \nu_A^2 \sum_{i_1, j_1}\sum_{i_2, j_2} \sigma_{i_1 j_1}^{(A)}\,\sigma_{i_2 j_2}^{(A)}\, e^{i\mathbf{k}\cdot(\mathbf{r}_{i_1 j_1} - \mathbf{r}_{i_2 j_2})}$$

**Length scale extraction:**

$\hat C_{AA}(k)$ has a peak at some $k_\star$. The corresponding length is

$$\xi_{AA} = \frac{1}{|k_\star|}$$

==**This is the size of segregated A-domains.**== The peak existing at all signals microphase separation.

> [!success] What we actually pull out
> - **$\xi_{AA}(t)$** — domain-size growth in time (segregation kinetics)
> - **$k_\star$** — preferred wavevector (lengthscale of the pattern)
> - **Peak shape of $\hat C_{AA}(k)$** — order of the transition (sharp peak = first-order-like, broad = critical)
> - **A/B clustering metrics** beyond just AA (also $C_{AB}$, $C_{BB}$, contrast $\hat C_{AA} - \hat C_{AB}$)

---

## 4. Visual schematic (from boards)

```
   MULTI-CHAIN BOX                  A----B
   ┌─────────────────┐             ●----●
   │  ╭─╮  ╲╱╲   ╮  │              r - r'  →  k
   │ ╱   ╲ ╱   ╲╱  ╲│
   │╲    ╱╲   ╱╲   ╱│             σᵢⱼ^(A) = 1 if bead (i,j) is A
   │ ╲╱ ╱  ╲ ╱  ╲╱  │                       0 otherwise
   │   ╱    ╳    ╲  │
   │  ╲    ╱╲    ╱  │             φ_A(r) = ν_A Σᵢ Σⱼ σᵢⱼ^(A) δ(r-rᵢⱼ)
   │   ╲__╱  ╲__╱   │
   └─────────────────┘             φ̂_A(k) = ν_A Σᵢ Σⱼ σᵢⱼ^(A) exp(ik·rᵢⱼ)
   (many chains, A & B beads)


   STRUCTURE FACTOR  Ĉ_AA(k)             DOMAIN GROWTH ξ_AA(t)
   
        Ĉ_AA                                ξ_AA
         │                                   │
         │       ╱╲                          │              ╱
         │      ╱  ╲                         │           ╱
         │     ╱    ╲___                     │       ╱
         │    ╱         ╲___                 │   ╱
         │___╱               ╲__             │╱
         └────────────────── k               └─────────── t
              ↑                                
              k★                              
                                              ξ_AA = 1/|k★|
   peak at k★ → length scale of           grows with time
   segregated domains                      → microphase separation
```

---

## 5. Why Brownian dynamics, not MMC

> [!warning] MMC artifacts the feedback flagged
> - **Accept/reject creates pseudo-frustration.** When moves are routinely rejected at low temperature, MMC dynamics differ qualitatively from the underlying Langevin dynamics (slower relaxation, modified barrier crossing). At a quench, this can masquerade as physical glassiness when it's really sampling-protocol slowness.
> - **No physical timescale.** "MC sweeps" are not seconds. We can't compare to experimental relaxation times directly.
> - **No hydrodynamic / inertial regime.** BD captures real diffusive motion of the polymer in a solvent.

**Brownian dynamics equation** (overdamped Langevin):

$$\zeta\,\frac{d\mathbf{r}_{ij}}{dt} = -\nabla_{ij} U(\{\mathbf{r}\}) + \boldsymbol{\eta}_{ij}(t), \qquad \langle\eta_\alpha(t)\eta_\beta(t')\rangle = 2\zeta k_B T\,\delta_{\alpha\beta}\delta(t-t')$$

where $\zeta$ is the friction coefficient and $\boldsymbol{\eta}$ is thermal noise satisfying fluctuation–dissipation. ==**This gives us real time.**==

**Implementation choices (in order of effort):**
1. Existing BD frameworks: **HOOMD-blue**, **LAMMPS** (DPD/BD), **OpenMM** (custom integrator), **ESPResSo**.
2. Custom BD integrator on top of current code (Euler–Maruyama, $\Delta t \sim 10^{-3}\tau$).

---

## 6. Phase-transition characterization

> [!info] What "nature of phase transition" means here
> Microphase separation of A/B copolymers can be **continuous** (Leibler-type, ODT — order-disorder transition) or **first-order** (strong-segregation regime). The order is read off from:
> - **Sharpness of $\hat C_{AA}(k_\star)$** — divergent peak = continuous; jump = first-order.
> - **Hysteresis** in heating/cooling cycles.
> - **Distribution of order parameter** $\hat\phi_A(k_\star)$ — Gaussian (continuous) vs bimodal (first-order).
> - **Finite-size scaling** of the peak height.
>
> Sequence correlations ($\kappa$ analog — i.e., block vs random copolymer at fixed A-fraction) are predicted to shift the ODT temperature and may even change the order. Connecting the i.i.d. vs correlated story to the ODT character is the **theoretical hook** that makes this lab-affiliation-worthy.

---

## 7. Compute strategy

**Single-workstation regime is over.** Multi-chain BD with density-field analysis at the scale needed (hundreds of chains, $10^6$–$10^7$ BD steps, multiple $\kappa$, multiple compositions) needs real compute.

| Option | What it gets us |
|---|---|
| **AWS spot GPU instances** | Run HOOMD-blue / OpenMM at scale; ~$0.50/GPU-hr on spot pricing. |
| **Stanford Sherlock / Farmshare** | Free CPU + GPU access for academic; right place to start. |
| **NSF-ACCESS allocations** | Free large-scale HPC if affiliated with a PI. |
| **Field-theoretic treatment (SCFT)** | Self-consistent field theory analytically computes $\hat C_{AA}(k)$ for given block architecture — sanity check / fast scan, then BD for the dynamics. |

> [!tip] Recommended split
> SCFT first to map the equilibrium phase diagram (cheap, analytical). Then BD on AWS or Sherlock to study the **kinetics** of segregation after a quench (the actual aging story).

---

## 8. Implementation roadmap

> [!success] Concrete code & analysis steps

**Phase 1 — model rewrite (1–2 weeks)**
- [ ] Add `bead_type` array to chain initialization (replace edge-disorder).
- [ ] Add type-specific LJ parameters: $A_{AA}, A_{BB}, A_{AB}, R_{\alpha\beta}$.
- [ ] Multi-chain support: list of chains in a periodic cubic box (use minimum-image convention).
- [ ] Sequence generators: pure-A, pure-B, alternating, random copolymer (composition $f_A$), block copolymer (block length $\ell$), correlated copolymer (Markov with persistence $\pi$).

**Phase 2 — Brownian dynamics (1–2 weeks)**
- [ ] Replace `metropolis_sweep` with overdamped Langevin integrator (Euler–Maruyama).
- [ ] Bonded forces: harmonic + FENE option.
- [ ] Non-bonded forces: type-pair LJ with cutoff and neighbor lists.
- [ ] Validate against MMC equilibrium statistics at the same temperature.

**Phase 3 — density-field observables (1 week)**
- [ ] Grid the simulation box (e.g., $64^3$ Fourier grid).
- [ ] Compute $\phi_A(\mathbf{r})$, $\phi_B(\mathbf{r})$ on the grid each frame.
- [ ] FFT → $\hat\phi_\alpha(\mathbf{k})$.
- [ ] Spherically average $|\hat\phi_\alpha(\mathbf{k})|^2$ → $S_{\alpha\alpha}(k)$.
- [ ] Locate $k_\star$, extract $\xi_{AA}(t)$.

**Phase 4 — phase-transition characterization (1–2 weeks)**
- [ ] Temperature scan to locate ODT.
- [ ] Cumulant / Binder analysis on $\hat\phi_A(k_\star)$.
- [ ] Sequence-correlation scan ($\kappa$ analog): does block vs random shift the ODT?
- [ ] Quench-then-watch: $\xi_{AA}(t) \sim t^{1/z}$ — extract growth exponent (testing Lifshitz–Slyozov $z=3$ or other coarsening laws).

**Phase 5 — the new headline plot**
The microphase-separation $\xi_{AA}(t)$ curves under different sequence statistics, showing how sequence correlation tunes the **kinetics of pattern formation**. This is the multi-chain analog of the original $\chi_4^\star(t_w)$ ensemble-comparison story.

---

## 9. What survives from the current paper

- **The central claim** — sequence correlations tune dynamics independently of marginal statistics — translates directly: random vs block copolymer at the same $f_A$ should give different $\xi_{AA}(t)$, different $k_\star$, different ODT.
- **The KWW / non-Gaussianity / subdiffusion analysis** still applies, now per-chain inside a melt.
- **The $\kappa$-scan methodology** maps onto a sequence-correlation scan.
- **The two-time / aging framework** survives — quench from disordered (high-$T$) state and watch domains coarsen.

What gets retired:
- Single-chain $N=30$ in vacuum (becomes the "minimal cartoon" version, kept in the appendix).
- Edge-disorder $\eta_{ij}$ (replaced by bead-type $\sigma_{ij}^{(A)}$).
- MMC dynamics (replaced by BD).
- Contact-overlap $Q$ as the primary observable (replaced by $\hat C_{AA}(k)$ and $\xi_{AA}(t)$).

---

## 10. Open questions to bring back to the meeting

> [!question] Things I still want clarified
> 1. **Box geometry & boundary:** periodic cubic, slab, or droplet? Periodic is standard; droplet is more biologically relevant for condensates.
> 2. **Composition:** symmetric ($f_A = 0.5$) or asymmetric? Asymmetric gives lamellar → cylindrical → spherical pattern progression.
> 3. **Solvent treatment:** implicit (BD in vacuum with friction $\zeta$) or explicit (DPD, MARTINI)?
> 4. **Field-theoretic vs particle-based:** SCFT only gives equilibrium; we want kinetics. Is the suggestion to do both (SCFT for phase diagram, BD for dynamics) or pick one?
> 5. **Connection to original $\kappa$ study:** do we frame this as "edge-disorder model was a stepping stone, the real result is in copolymer melts," or as parallel results in different model systems?
> 6. **Lab affiliation target:** is this aimed at a specific group (polymer-physics theory / IDP biophysics / SCFT methods)?

---

## Appendix — Whiteboard transcription

**Board 1 (left):**
- Sketch of polymer chains in a box (multi-chain melt)
- Two-bead schematic: A — B, separated by $\mathbf{r} - \mathbf{r}' \to \mathbf{k}$
- $\sigma_{ij}^{(A)}$ definition (indicator)
- $\phi_A(\mathbf{r}) = \nu_A \sum_i \sum_j \sigma_{ij}^{(A)} \delta(\mathbf{r} - \mathbf{r}_{ij})$
- $\hat\phi_A(\mathbf{k}) = \nu_A \sum_i \sum_j \sigma_{ij}^{(A)} \exp(i\mathbf{k}\cdot\mathbf{r}_{ij})$
- $\langle \phi_A(\mathbf{r})\phi_A(\mathbf{r}')\rangle = C_{AA}(\mathbf{r}-\mathbf{r}')$
- Decay sketch: $C_{AA}$ vs $|\mathbf{r}-\mathbf{r}'|$ with correlation length $\xi_{AA}$ marked

**Board 2 (extension):**
- $\xi_{AA}$ vs $t$ plot (growth)
- $\hat C_{AA}(\mathbf{k}) = \langle \hat\phi_A(\mathbf{k})\hat\phi_A(-\mathbf{k})\rangle = \nu_A^2 \sum\sum\sigma\sigma\exp[i\mathbf{k}\cdot(\mathbf{r}-\mathbf{r})]$
- Plot of $\hat C_{AA}(k)$ with peak at $k_\star$
- $\xi_{AA} = 1/|k_\star|$
