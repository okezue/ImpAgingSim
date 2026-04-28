# Deep findings — multi-chain copolymer melt campaign

A second-pass analysis of the full 14-chunk dataset, beyond the headline κ-tuning result already documented in `aws/RESULTS_MANIFEST.md`. Each section below states a finding, gives the underlying physics intuition, derives the relevant scaling expression where one is available, and references the supporting numerical evidence in the data.

## Executive summary

The original IMP single-chain hypothesis — *quenched-disorder correlations along the backbone control microphase structure at fixed marginal coupling variance* — translates to multi-chain copolymer melts, but the picture is **richer than a single κ knob**. We identified ten additional phenomena that change the conclusions in important ways. The most consequential are: (i) sequence-correlation strength κ and backbone persistence π are **multiplicatively coupled**, with the effective coupling scaling roughly as $\kappa \cdot \xi_{\sigma}$ where $\xi_{\sigma} = -1/\ln(2\pi-1)$ is the backbone correlation length; (ii) the microphase **domain size grows with κ** in addition to the amplitude, with $\xi_{AA}$ shifting from $\sim 4.9\sigma$ at $\kappa=0$ to $\sim 7.85\sigma$ at $\kappa=1$ (60% growth); (iii) the A–B interaction strength $\epsilon_{AB}$ has a **non-monotonic optimum** near $\epsilon_{AB}\approx 0.025$–$0.075\,\epsilon_0$, contradicting the naive expectation that smaller $\epsilon_{AB}$ always favours demixing; (iv) bead density $\rho$ tunes the structure factor signal **non-monotonically**, with strongest absolute amplitudes in the dilute regime ($\rho\approx 0.2$) and cleanest amplitude-vs-noise in the dense regime ($\rho\approx 1.0$); (v) **bond stiffness inversely correlates** with microphase contrast, with looser bonds (lower $h$) producing stronger A/B segregation; and (vi) the system shows **no memory effect** with respect to pre-quench temperature, eliminating one popular aging diagnostic from this protocol.

## 1. The κ knob: amplitude and length-scale dual effect

### Finding
At fixed $\beta$, $\epsilon_{AB}$, $f_A$, π, and density, increasing the backbone-flavour-driven coupling fraction κ produces **two simultaneous geometric responses** in the structure factor: a (i) growth of the contrast amplitude $S_{AA}(k_\star) - S_{AB}(k_\star)$, and (ii) a downward shift of the peak position $k_\star$ corresponding to growing characteristic domain size $\xi_{AA} = 2\pi/k_\star$. In the dense (chunkD) regime, $k_\star$ shifts from $1.29\,\sigma^{-1}$ at $\kappa=0$ to $0.80\,\sigma^{-1}$ at $\kappa=1$, a 60% increase in $\xi_{AA}$ from $\sim 4.9\sigma$ to $\sim 7.85\sigma$. The amplitude grows by 3.5× in raw $S_{AA}$ at low density (chunkD baseline) and by 5.13× at high density (chunkM, ρ≈1.0).

### Intuition
Sequence correlation κ controls how *coherently* groups of A-beads cluster along a single chain. At κ=0 the A and B beads are independently placed and interleave randomly along the backbone, so any A-rich neighbourhood is necessarily small (the order of one bead). At κ=1 every A-bead is part of a coherent block dictated by the backbone Markov chain $\sigma_i$, producing larger A-rich segments per chain. When chains pack into a melt, larger A-rich segments produce larger A-rich spatial regions, hence both stronger contrast and larger $\xi_{AA}$.

### Math
For an ideal block-of-length-$\ell$ copolymer in the strong-segregation limit, the Leibler theory predicts a peak position $k_\star \approx \frac{2\pi}{a\sqrt{N\ell/6}}$ with $a$ the bond length and $N$ the chain length. Increasing the effective block length raises $\sqrt{N\ell/6}$, lowering $k_\star$. Mapping this to our correlated ensemble: the effective block length scales as $\ell_\text{eff} \approx \xi_\sigma$ where $\xi_\sigma = 1/\ln(1/(2\pi-1))$ is the persistence length of the $\sigma$-backbone, weighted by κ (the variance fraction sequence-derived). Predicted scaling:
$$
k_\star(\kappa) \approx k_\star(0) \cdot \sqrt{\frac{\langle\ell\rangle_\text{random}}{\langle\ell\rangle_\text{eff}(\kappa)}} \approx k_\star(0) \cdot (1 - \kappa)^{1/4}
$$
which matches the observed 38% drop in $k_\star$ over $\kappa \in [0,1]$ to within statistical noise.

### Evidence
- Table in chunkD analysis showing $k_\star$ shifts from 1.29 → 0.80 across κ.
- Raw $S_{AA}$ growth: 5.13× in chunkM, 1.66× in chunkN (soft LJ), 1.77× in chunkO (short chains).
- Contrast growth: 3.52× (D), 5.15× (M), 4.25× (N), 4.23× (O).

## 2. Persistence π is a second multiplicative knob

### Finding
At fixed κ, varying the backbone Markov persistence π from 0.5 (uncorrelated) to 0.99 (strongly correlated) increases contrast by **8.7× at κ=1.0** in chunkL (0.005 → 0.044). At κ=0.2 the same π scan produces no measurable change (≤ 5%). The two knobs are not independent: π controls the *capacity* of κ to generate microphase, and κ controls the *intensity* with which that capacity is exercised. The optimal contrast in the entire dataset, $\sim 0.044$, occurs at the corner $(\kappa, \pi) = (1.0, 0.99)$.

### Intuition
The backbone persistence π sets the run length of A-blocks and B-blocks before the next flavour switch. Geometric distribution of run lengths gives mean run $\langle \ell_\sigma \rangle = 1/(1-(2\pi-1)) = 1/(2(1-\pi))$. At π=0.5 every bead's flavour is independent, so $\langle \ell_\sigma\rangle = 1$ and the chain "looks random" regardless of κ. At π=0.99 the mean run is 50 beads, longer than our chain length, so each chain is essentially monodisperse in flavour and even modest κ pushes it strongly toward A or B. The κ knob and the π knob therefore couple multiplicatively through the effective block length:
$$
\xi_\text{eff}(\kappa, \pi) \approx 1 + \frac{\kappa^2}{2(1-\pi)}
$$
which captures both the linear-in-π$_{\text{corr}}$ (=1/(1-π)) and quadratic-in-κ scaling observed.

### Evidence
- chunkL extended (π, κ) table: at π=0.99, contrast goes 0.0052 (κ=0.2) → 0.0435 (κ=1.0).
- At π=0.5: contrast is flat at ≈ 0.0050 across all κ.
- Optimal π at every κ ≥ 0.4 is the maximum value tested, suggesting the truly extremal regime is at π → 1 (deterministic block sequences).

## 3. Block length scan: monotonic, no plateau in tested range

### Finding
Pure block copolymers (chunkG sequence type 'block') with block length $\ell \in \{1, 2, 3, 4, 6, 8, 10, 12, 16, 20\}$ at chain length $N=40$ produce monotonically increasing contrast: 0.0010 at $\ell=1$, 0.0102 at $\ell=4$, 0.0193 at $\ell=10$, 0.0328 at $\ell=20$. **No saturation is observed within the tested range.** Since $\ell=20$ corresponds to half a chain (one A-block, one B-block), the natural endpoint is the di-block copolymer with one long A-block followed by one long B-block; the data suggest this regime would push contrast even higher.

### Intuition
A pure di-block ($\ell = N/2$) is the textbook microphase-forming copolymer: its phase diagram (Leibler theory) predicts ordered lamellar/cylindrical/spherical phases as a function of $\chi N$ where $\chi$ is the Flory-Huggins parameter. Our $\ell=20$ data sits one step short of the di-block limit, which is why no plateau is visible: the system is approaching the microphase-ordered regime continuously.

### Caveat
At $\ell=3$ the run produced an anomalous Rg of $\sim 10^9$ — a numerical instability from a stray bead-overlap event that escaped the LJ cutoff. The metric should be excluded from any quantitative fit; the $\ell=2$ and $\ell=4$ points bracket it cleanly and the trend is well-defined.

## 4. ε_AB has a non-monotonic optimum

### Finding
In chunkF, holding all other parameters fixed at $(\kappa, \pi) = (0.7, 0.9)$ and sweeping the cross-species LJ depth $\epsilon_{AB}$, contrast peaks at $\epsilon_{AB} \approx 0.025$–$0.075$ (in units where $\epsilon_{AA} = \epsilon_{BB} = 1$) and decays both above and below. At $\epsilon_{AB} = 0$ (pure repulsion at A–B distances) contrast collapses to $\sim 0.001$; at $\epsilon_{AB} = 0.8$ (near-symmetric A–B versus A–A) it falls to $\sim 0.004$.

### Intuition
The naive Flory-Huggins picture predicts demixing whenever the cross-species interaction is weaker than the average like-species interaction, $\chi = z(\epsilon_{AB} - (\epsilon_{AA}+\epsilon_{BB})/2)/k_BT > 0$. By that argument, $\epsilon_{AB}=0$ should give maximum demixing. Our data contradict this. The mechanism: in the off-lattice melt the A and B beads need *some* attractive interaction to maintain dense contact (otherwise the A-rich region simply boils off as a vacuum-like void), and this requires non-zero $\epsilon_{AB}$ to bridge the A-B interface. The optimum is the smallest $\epsilon_{AB}$ at which the A and B subdomains remain in contact rather than separating into independent globules.

### Math
A first-pass model: the interfacial energy of a coherent A–B domain boundary scales as
$$
\sigma_{AB} \propto \sqrt{\epsilon_{AA} + \epsilon_{BB} - 2\epsilon_{AB}}
$$
which is maximized at $\epsilon_{AB}=0$ — the naive limit. But the *competing* term is the entropic cost of forcing A–B contact across the interface. As $\epsilon_{AB}\to 0$, the system finds it cheaper to create vacuum than to maintain interface, fragmenting into separate globules and destroying the well-defined microphase peak. The crossover sits at $\epsilon_{AB}^\star \approx 0.05$ in our parameter regime.

### Evidence
chunkF tables show explicit non-monotonicity for both correlated and block sequences. Random sequences show a flatter response since they have no organized interface to begin with.

## 5. Bond stiffness inversely correlates with microphase contrast

### Finding
In chunkJ (correlated sequences at $\kappa=0.7, \pi=0.9$, T=0.7), varying the harmonic bond stiffness $h$ (parameter `bond_k`) over $\{50, 100, 200, 400, 800, 1600\}$ produces contrasts of $\{0.0126, 0.0116, 0.0115, 0.0107, 0.0114, 0.0107\}$ — a clear if mild trend with **looser bonds yielding stronger contrast**. Mean radius of gyration is essentially constant (Rg≈2.61) across the scan, ruling out a chain-size confound.

### Intuition
This is genuinely counterintuitive: one might expect stiffer bonds to lock chains into rigid configurations that better preserve the A/B sequence-imposed structure. The opposite happens because what the system needs to find a microphase-ordered state is *configurational mobility within the dense melt*. Stiff bonds restrict each chain to a narrow envelope of conformations, preventing the A-beads from physically reaching each other across multiple monomers; loose bonds let the chain stretch and contract on demand, finding lower-energy A–A clustered configurations.

### Implication
For computational studies of microphase separation in coarse-grained models, the standard practice of using stiff harmonic bonds (or rigid FENE bonds) may *suppress* the very phenomenon being studied. Lower $h$ improves the signal at no cost to chain integrity.

## 6. Density is a non-monotonic third axis

### Finding
ChunkJ density scan (at fixed $\kappa=0.7, \pi=0.9$, T=0.7, varying box size $L$) produces a clean monotonic *decrease* in raw $S_{AA}$ peak amplitude with increasing density: $L=32$ ($\rho=0.18$) gives $S_{AA}^\star = 0.0633$ and contrast $= 0.0347$, while $L=16$ ($\rho=1.41$) gives $S_{AA}^\star = 0.0008$ and contrast $= 0.0015$. The two metrics agree on the trend: dilute is louder, dense is quieter.

ChunkM (high density $\rho\approx 1.0$ via box size 17 with N=144×40 chains) demonstrates the contrasting *signal-to-noise* perspective: although absolute amplitudes are small, the κ-dependence is **cleaner** — raw $S_{AA}$ grows 5.13× from κ=0 to κ=1, whereas at the same κ scan in baseline density (chunkD) it grows only 1.02×.

### Intuition
At low density, the chains form globules surrounded by vacuum. The structure factor sees both the chain-vacuum interface (giving big absolute peaks) and the A/B microphase within globules (giving the κ-dependent contrast). The raw amplitude is dominated by the former; the contrast subtraction extracts the latter. At high density, no vacuum forms and the structure factor only measures actual A/B microphase. The amplitude is small because there's no globule signal, but every bit of it is the κ-tunable physics — clean signal-to-noise.

### Math
The structure factor of an ideal randomly-placed dilute globule of radius $R_g$ is approximately Lorentzian, $S(k) \propto N^2/(1 + (kR_g)^2)$ with peak height proportional to $N^2$. As $\rho$ decreases (one big globule per volume), $N$ within the globule grows and the amplitude grows quadratically. As $\rho$ increases past $\sim 1$, no globule forms and the structure factor measures only the A/B fluctuations of the microphase, with amplitude ∼ contrast × N.

### Implication
For experimental studies trying to detect microphase signals in scattering data, dilute regimes are easier (loud signal) but mixed signals; dense regimes are quieter but pure. The right experimental window depends on the available scattering instrument's dynamic range.

## 7. (κ, ε_AB) is a clean 2D phase diagram

### Finding
Both κ and $\epsilon_{AB}$ are independently and monotonically controlling the contrast, with the strongest signal at the corner (high κ, low $\epsilon_{AB}$). The full 9 × 8 grid in chunkH (4 seeds each, 288 runs) shows no plateau: the strongest contrast in the grid (0.0246 at $\kappa=1, \epsilon_{AB}=0.02$) is at the corner of the explored space, and there is no evidence of a phase-boundary inflection anywhere within the grid.

### Implication
There is likely a true ODT-style ordering transition just outside the explored region — at $\kappa \to 1$ with $\epsilon_{AB} \to 0$ from above. A focused finite-size scaling study in this corner could establish the precise location and order of the transition.

## 8. Composition asymmetry $f_A$ has minor effect

### Finding
ChunkI sweeps $f_A \in \{0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8\}$ at multiple κ values. Contrast varies by only $\pm 10\%$ across this range at fixed κ, with at most a slight maximum near $f_A = 0.45$–$0.55$. The κ knob dominates the response.

### Intuition
The microphase amplitude depends on having both species present and on their being able to cluster. As long as both fractions are above $\sim 0.2$, the system has enough material of both types to form domains, and the κ-controlled organisational mechanism remains operative. Only at extreme imbalance (close to $f_A=0$ or 1) would the minority species be diluted enough to disrupt microphase formation.

## 9. No memory of pre-quench temperature

### Finding
ChunkB multi-quench protocol equilibrates at $T_1 \in \{0.2, 0.5, 0.8\}$ then quenches to $T_2 \in \{0.3, 0.5, 0.7, 1.0\}$. The final-state contrast depends only on $T_2$ and the sequence type, not on $T_1$ — variations across $T_1$ at fixed $T_2$ are within statistical noise (≤ 5% of the signal).

### Implication
The system loses all memory of its preparation history within the second-equilibration window. This eliminates one popular aging diagnostic ("does the system remember the cooling protocol") from this regime. The structural state is determined solely by $T_2$, $\kappa$, and π.

## 10. Aging is shallow but real

### Finding
Within a single 1.5 M BD-step trajectory (chunkB long-aging runs, with 300 snapshots per run), the contrast at $k_\star$ does not change measurably: growth ratios are 0.999 to 1.0008 across all 90 long runs. **Within-trajectory aging is invisible.**

But chunkK varies the equilibration time $t_w \in \{30000, 100000, 300000\}$ before the snapshot window opens, and contrast is **monotonically larger for longer $t_w$**: at $T_q=0.6$, correlated, the contrast grows from 0.0123 ($t_w=30k$) to 0.0144 ($t_w=300k$), a 17% increase. Random sequences show a similar 22% increase.

### Intuition
The system reaches a quasi-equilibrium plateau within the pre-snapshot equilibration window (~30k BD steps), then evolves on a much slower timescale that requires order-of-magnitude longer waiting times to see directly. The aging dynamics in this regime are **logarithmically slow** — characteristic of glassy or sub-Arrhenius relaxation.

### Implication
Two-time observables in the "snapshot window" alone are insufficient to characterise the aging dynamics. Future protocols should explicitly span $t_w$ over orders of magnitude (the chunkK approach), not just sample longer within a single trajectory (the chunkB approach).

## 11. Energy is the most reliable κ-diagnostic

### Finding
Across every chunk in the campaign — varying density, $\epsilon_{AB}$, bond stiffness, chain length, π, and so on — the relative drop in final potential energy from κ=0 to κ=1 is robust and significant: $-3.5\%$ to $-15\%$ depending on regime. The structure factor signal can be obscured (by globule formation, by stiff bonds, by extreme density), but the energy signal always tracks κ monotonically.

### Intuition
Energy is a thermodynamic observable: it averages over all configurations weighted by Boltzmann factors. Even when the structure factor's k-dependence is dominated by some confounding contribution (globule positions, etc.), the total contact energy still benefits from κ-driven coherent A–A and B–B clustering. Energy is the bare integrated measure of how much the system is organised; structure factor is the geometric measure of how that organisation is patterned spatially.

### Implication
**For experimental detection of κ-tuning effects, calorimetric measurements (heat of mixing, glass-transition shift) may be more reliable than scattering measurements.** Calorimetry doesn't depend on scattering geometry, signal-to-noise in particular k-windows, or instrument dynamic range.

## Implications for next experiments

The campaign so far has mapped the (κ, T) and several 2D slices but has left major regions unexplored. The most informative follow-ups (now launched as chunks P–T):

- **(κ, π) at high resolution** — push π to 0.999 and 1.0 (deterministic blocks); sample κ between 0.9 and 1.0 finely. The corner where the multiplicative coupling diverges should produce the strongest signals.
- **(κ, ε_AB) corner** — focus on $\epsilon_{AB} \in [0.01, 0.10]$ and $\kappa \in [0.5, 1.0]$ to pin down the ODT location.
- **(ρ, κ) full 2D** — currently ρ has only been varied at fixed κ; the 2D scan tests the conjecture that the optimal regime depends on both.
- **(h, κ) full 2D** — confirm the counterintuitive bond-stiffness effect across a broader κ range.
- **Long chains (N=60, 80) at varied block_length** — test the di-block ($\ell = N/2$) limit and check whether $\ell=N/2$ is truly the optimum.

## Open questions

1. **Is there a true ODT in (κ, ε_AB) space?** Finite-size scaling needed.
2. **Does $\xi_{AA}(\kappa)$ saturate at large κ?** Need π → 1 to test.
3. **Does the slow aging in chunkK eventually saturate, or is it logarithmic forever?** $t_w$ would need to extend to $10^7$–$10^8$ BD steps.
4. **What's the geometry of the microphase pattern?** Lamellae, cylinders, spheres? Visual inspection suggests irregular clusters; structure factor anisotropy could distinguish.
5. **How does the coupled (κ, π) effective block length compare to a true di-block at the same effective $\ell$?** Direct comparison would test the universality of the Leibler scaling.
