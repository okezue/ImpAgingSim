# Briefing for the first in-person meeting with Prof. Spakowitz (September 2026)

Everything below is in the repository (`analysis/*/README.md`) and on Zenodo
(10.5281/zenodo.20499120 for the melt; 10.5281/zenodo.22819768 and
10.5281/zenodo.22819770 for the two applications). Numbers quoted are from those notes.

## 1. Status against the NYC-call work plan

| task | status |
|---|---|
| 1. Locate the mixed -> demixed transition at fixed `kappa = 0.5` | done: `eps_AB,c = 0.85 +/- 0.03` at `T* = 0.7`; identical for 144 and 288 chains; RPA spinodal 0.905 (`alpha ~ 2.0`) |
| 2. Extend to the dynamic structure factor | done: exact box-mode recorder; `S(q,tau)`, `F(q,tau)` for every shell; theory-comparison dataset at `T* = 1.5` and `2.0` (8 seeds, 5M steps, 1-tau resolution) |
| 3. Compute | Sherlock kit written and documented but unused; all 1,592 runs were done on xAI GPUs (OpenMM CUDA). Sherlock remains an option |
| 4. Background reading | ChemE 466 lectures 10-11 and the 2017 random copolymer paper: to be discussed |

## 2. The answer to Task 1, and how it maps onto the phase diagram

- The chi = 0 reference (`eps_AB = 1`) is mixed and matches the ideal single-chain form
  factor, so the sweep is anchored.
- Demixing sets in at `eps_AB = 0.85` (half-rise of `S_psi(q*)`), the same in two box sizes.
  The RPA instability, fitted on the mixed side, sits at 0.905: the real melt demixes about
  50% deeper in incompatibility than mean-field theory says, the expected sign of the
  fluctuation correction for random copolymers.
- Along `kappa` at `T* = 1.5`: no transition for `kappa <= 0.25`; `delta_eps_c` = 0.20, 0.15,
  0.115 for `kappa` = 0.5, 0.75, 1, with plateau amplitudes 120, 500, 1,020. RPA matches at
  `kappa = 0.5` and predicts instability too early for `kappa >= 0.75`. This is the (lambda,
  chi) phase boundary of the 2017 paper as the simulation sees it.
- Of the two fluctuation spikes he described: the homogeneous -> random microphase one is
  present for `kappa >= 0.5`; the random microphase -> structured (lamellar) one appears only
  for `kappa = 1` at strong incompatibility in this box (shell anisotropy 2.7, near the
  lamellar value 3). At `kappa = 0.5` it does not exist in the accessible range.

Open question from the notes now answered: the fluctuation observable. The peak-intensity
time variance and the seed variance both peak at the transition in an ergodic melt; the
per-mode non-Gaussianity ratio (2 for Gaussian amplitudes, 1 for a frozen pattern) is the
cleanest order-parameter-like indicator, and shell anisotropy is the lamellar diagnostic.

## 3. The answer to Task 2, and what to bring to the theory

- At `T* = 0.7` the composition field is arrested: no mode below `q ~ 0.5/sigma` relaxes
  even in a 6,250-tau window. The dynamics comparison has to be done at `T* >= 1.5`.
  (Question for him: was `T* = 0.7` deliberate, or should the whole program move to an
  ergodic temperature? The chi mapping `delta_eps / T*` collapses the temperatures, so
  nothing is lost.)
- At chi = 0 and `T* = 1.5` the composition relaxation is a single exponential
  (`beta = 0.99`) with `tau ~ q^-3.9` for `q R_g > 1.8`, the Rouse prediction, and a
  crossover toward diffusion below `q ~ 0.4`.
- Toward the transition only the low-q modes slow down; `Gamma(q) S(q)` stays roughly
  constant while `S(q_min)` grows eightfold. That is the structure a dynamic RPA of the
  composition field has (thermodynamic slowing with an incompatibility-independent kinetic
  coefficient), so it is the first quantitative target.
- Two puzzles for the theory: (i) below `q ~ 0.3` the relaxation time levels off near
  3,000 tau in both the 144- and 288-chain boxes instead of continuing as `q^-2`;
  (ii) warming from `T* = 1.5` to 2.0 speeds every mode by 2.2x, more than the 1.33x of a
  purely thermal mobility, so the friction constant in the theory has an activated part.
- Files ready for fitting: `analysis/dsf_theory/*/F_peak_by_condition.csv` (1-tau lag grid,
  mean and SEM over 8 seeds) and `tau_by_shell_condition.csv`; full mode-resolved data on
  Zenodo. If his theory needs lags below 1 tau or wavevectors above `1.5/sigma`, one more
  campaign with a shorter mode interval is a few GPU-hours.

## 4. Coarsening and arrest (relevant to the aging story)

Over 100,000 tau at `T* = 0.7`, quenches just past the crossover coarsen as
`S_psi(q*) ~ t^0.5` with `q* ~ t^-0.1` for three decades until the domains reach the box,
while a deeper quench arrests after ~10,000 tau at a finite domain spacing (~16 sigma).
Stronger incompatibility gives a smaller, frozen pattern. Every "plateau" seen in the
short Stage 1/2 runs was a snapshot of a still-growing pattern.

## 5. The application (the chromatin note), and what we found

Both applications the note mentions were built and run.

**Protein / sequence-programmed condensation** (360 runs; B loves B, A-B neutral, free
volume). Whether a shared condensate forms is decided by sequence blockiness, not sticker
strength: the boundary in the `eps_BB x kappa` plane is a horizontal line at `kappa ~ 0.5`.
Scattered patterns never phase separate even at `eps_BB = 3` (each chain collapses onto its
own stickers); blocky patterns condense at `eps_BB = 0.5` with open chains. Strong
attraction fragments the condensate kinetically.

**Chromatin / epigenetic memory** (160 runs; marks turn over at `k_off`, are written at
`k_on + k_fb H(n_B)` with `n_B` the marked neighbors within 1.5 sigma). Answers to the
headline questions:
1. *What happens to transient blobs?* Without feedback they nucleate and dissolve; their
   lifetime (two-time B-density correlation) is 16-290 tau and tracks the mark memory time
   `1/(k_on + k_off)`. At weak attraction they die by diffusion before the marks forget; at
   stronger attraction mark turnover is the limiting step.
2. *When is a blob stabilized?* When `k_fb / k_off >~ 10` (marginal at 3): a diagonal boundary
   in the turnover x gain plane, the same at both attractions tested.
3. *Then what?* The stabilized blob accretes marks until the marked fraction reaches a plateau
   set by the ratio alone (0.75 at 10, 0.93 at 30, complete conversion above 100). The
   minimal positive-feedback model has no intrinsic domain-size control.

The open questions in his note, as the data answers them:

| question | what the model says |
|---|---|
| What sets the flip rate, uniform or environment-dependent? | With a uniform `k_off` and an environment-dependent `k_on`, the ratio of the two is the only control parameter of stabilization; the biology's environment dependence of *erasure* (demethylases recruited to open chromatin) is the natural next ingredient and would provide the missing size control |
| Does writer recruitment depend on blob size, lifetime or local B fraction? | Local B count alone (Hill on `n_B`) is enough to produce the transient/stabilized dichotomy; size enters only through the surface-to-volume balance of writing vs turnover |
| Explicit or implicit solvent? | Free volume (bead density 0.2) as implicit solvent was sufficient to get B condensation with open A; no explicit solvent needed for these questions |
| How many mark types? | Two states suffice for the first two questions; the third finding (no size control) is where a second, antagonistic mark becomes necessary |
| What observable distinguishes a stabilized from a transient blob? | The two-time density correlation at the blob wavevector, exactly as he anticipated: equal-time cluster statistics of the transient and marginal regimes overlap, while blob lifetime separates them (below the mark memory time vs above it vs beyond the window) |

## 6. Proposed next steps (his call)

1. **Theory comparison at `T* = 1.5`/2.0** on the dataset in `analysis/dsf_theory/`; decide
   the lag grid and wavevector range his theory needs.
2. **Size control in the chromatin model**: three candidate ingredients, each a one-line
   change to `melt/marks.py` and a ~2 GPU-hour campaign: (a) a finite writer pool (global
   conservation), (b) an antagonistic A-stabilizing mark with its own feedback, (c)
   sequence-encoded nucleation and boundary sites (the underlying `kappa` pattern as a
   template). Each predicts a different domain-size law; that is the discriminating
   experiment.
3. **Replication dilution** as periodic halving of marks (an effective `k_off` burst), to
   connect `k_fb / k_off >~ 10` to cell-cycle timescales.
4. **Protein mode**: the liquid-vs-arrested character of the condensate as a function of
   sticker strength, from `F_BB(q*, t)`, to connect with the aromatic-patterning literature.

## 7. Practical

- Weekly check-ins from September; publication as the goal for the year. Two possible
  papers are visible in the data: the transition/dynamics program (with his theory), and
  the epigenetic-memory model.
- Compute: xAI GPUs are available and fast (a 250k-step run in 30 s; the whole 1,592-run
  program took under a day of wall time). Sherlock is documented as a fallback.
