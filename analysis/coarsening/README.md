# Coarsening of the arrested microphase at T* = 0.7 (campaign `coarsening_T0p7`)

Twelve runs: `kappa = 0.5`, `eps_AB` in {0.8, 0.7, 0.5} (`delta_eps` = 0.2, 0.3, 0.5, all on
the demixed side of the `eps_AB = 0.85` crossover), four seeds, 144 x 40 beads, `T* = 0.7`,
20,000,000 production steps (100,000 tau, eighty times Stage 1), box modes every
4,000 steps. Code commit `b35e4be`, 8 x H100 (~45 min per run).
`scripts/coarsening_analysis.py` bins `S_psi(q, t)` in logarithmic time; `q*(t)` is the
intensity-weighted wavevector of the three strongest shells so it moves continuously
between the discrete box shells.

## Result 1: near the crossover the pattern coarsens for three decades

At `eps_AB = 0.8` and `0.7` the peak amplitude grows from 12 at 100 tau to 200 at
`10^5` tau with `S_psi(q*, t) ~ t^0.49` over the last decade, and the peak wavevector
falls monotonically from 0.6 to `q_min = 0.29` (`q* ~ t^-0.10`), i.e. the domain
spacing grows from ~10 sigma to the box size. Only then does `S` saturate, at ~200:
the box, not the dynamics, ends the growth. Every plateau reported for this
temperature in Stages 1 and 2 (40-55 after 5,000 tau) was therefore a snapshot of a
still-growing pattern, as the positive drifts there already suggested.

## Result 2: a deeper quench arrests at a finite domain size

At `eps_AB = 0.5` the amplitude rises faster at first but saturates near
`S_psi(q*) ~ 100` after `~10^4` tau (`t^0.11` over the last decade) with `q*` stuck at
0.38 (domain spacing ~16 sigma), well above `q_min`; the lowest-shell intensity
`S_psi(q_min, t)` keeps creeping upward but stays a factor 2 below the two milder
quenches. Stronger incompatibility therefore produces a smaller, frozen pattern
rather than a larger one: the melt arrests before the domains reach the box, which
is the dynamical-arrest picture of the earlier aging study now seen directly in the
composition field.

## Result 3: contrast grows faster than domain size

Over the growth regime `S_psi(q*)` rises as `t^0.5` while `q*` falls only as
`t^-0.1`. If the pattern were coarsening at fixed contrast the amplitude would scale
with the domain volume, `S ~ q*^-3 ~ t^0.3`; the extra growth means the composition
contrast inside the domains is still increasing. Segregation strengthens and the
domains grow at the same time, and neither is finished at `10^5` tau for the milder
quenches.

## Files

| file | contents |
|---|---|
| `coarsening_binned.csv` | per (eps_AB, log time bin): `S_psi(q*)`, `q*`, `S_psi(q_min)`, mean and SEM over four seeds |
| `coarsening_fits.json` | growth exponents over the last decade |
| `coarsening.*` | the three panels |
| `coarsening_T0p7/` | campaign tables, `transition_summary.json` (trailing-window statistics), `manifest.json`, figures |
