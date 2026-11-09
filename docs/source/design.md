# Simulation design: the values the names do not carry

A parser in {mod}`polyphys.manage.parser` records what an artifact's name
carries, plus a few values it fixes itself. The rest of a study's design is
set in the LAMMPS input scripts that ran its simulations: the interaction
strengths, the cutoffs, the bonds and angles, the masses, and the
thermostat. This page lists those values for each study that has a parser.

The values come from the production inputs of the research code that ran
the simulations, and from the generator scripts that wrote each run's
input. Where those inputs differ between sets of simulations, the page gives
each value and the set that used it. The page describes the scripts; it has
not been checked against the copy of the input kept in each run's
directory.

## Notation and units

Every input uses LAMMPS reduced units (`units lj`): lengths are in the
monomer diameter, energies in the Lennard-Jones (LJ) strength
$\epsilon = 1$, masses in the monomer mass, and the thermostat temperature
is $k_B T = \epsilon$.

$a_m = 1$, $a_l$, $a_c$
: The diameters of a monomer (the small monomer of a two-size chain), a
  large monomer (keyword `al`), and a crowder (keyword `ac`).

$\epsilon_{ij}$, $\sigma_{ij}$, $r_{c,ij}$
: The strength, size, and cutoff of the LJ interaction between species $i$
  and $j$.

$\epsilon_w$, $\sigma_w$, $r_{c,w}$
: The strength, size, and cutoff of the LJ interaction between the
  cylindrical wall and one species.

$K$, $R_0$, $r_0$, $\theta_0$
: The constants of a bond or angle potential, in the LAMMPS form given with
  each potential below.

## Shared by every study

Each section below lists only where a study departs from these choices.

Pair interactions
: Every pair of particles interacts through the LJ potential of
  `pair_style lj/cut` ([LAMMPS](https://docs.lammps.org/pair_lj.html)),

  $$
  E(r) = 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r}\right)^{12}
  - \left(\frac{\sigma_{ij}}{r}\right)^{6}\right], \qquad r < r_{c,ij},
  $$

  with $\epsilon_{ij} = 1$, $\sigma_{ij} = (a_i + a_j)/2$, and
  $r_{c,ij} = 2^{1/6}\sigma_{ij}$. The cutoff sits at the minimum of the
  potential, so every pair, crowders included, is purely repulsive and
  feels the forces of the Weeks-Chandler-Andersen potential. The energies
  are not shifted to zero at the cutoff unless a section says so.

Bonds
: A bond between species $i$ and $j$ is a finitely extensible nonlinear
  elastic (FENE) bond, `bond_style fene` ([LAMMPS](https://docs.lammps.org/bond_fene.html)),

  $$
  E(r) = -\frac{1}{2} K R_0^2 \ln\left[1 - \left(\frac{r}{R_0}\right)^2\right]
  + 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
  - \left(\frac{\sigma}{r}\right)^{6}\right] + \epsilon,
  $$

  with $K = 30/\sigma_{ij}^2$, $R_0 = 1.5\,\sigma_{ij}$, $\epsilon = 1$,
  and $\sigma = \sigma_{ij}$. With `special_bonds fene`, two bonded
  particles do not also interact through the pair potential. There is no
  angle potential except in the H-NS studies.

Thermostat and integrator
: A Langevin thermostat, `fix langevin` ([LAMMPS](https://docs.lammps.org/fix_langevin.html)), acts
  on every particle at $T = 1$ with the damping time $10$: each particle
  of mass $m$ and velocity $v$ feels the friction $-(m/10)\,v$ and a matching
  random force. It is switched on before equilibration and stays on in
  production: the production inputs either keep it or set it again with
  the same values, and integrate with `fix nve`. Production therefore
  samples constant temperature, not constant energy. Equilibration
  integrates with `fix nve/limit 0.01`, which caps each particle's
  displacement per step at 0.01.

Recentering
: `fix recenter` keeps the center of mass of the chain at the origin, along
  the axis of a cylinder and in all three directions in a cube, by shifting
  every particle by the same amount at each step.

Preparing the start configuration
: Several studies run a separate preparation input before the production
  input. Despite its name ("minimize"), it is a short Langevin run, not an
  energy minimization. In the TransFoci studies its bonds are harmonic,
  `bond_style harmonic`
  ([LAMMPS](https://docs.lammps.org/bond_harmonic.html)),
  $E = K (r - r_0)^2$ with $K = 80/\sigma_{ij}^2$ and $r_0 = \sigma_{ij}$;
  the production input replaces them with FENE bonds.

## TransFociCyl

`TransFociCyl`: a ring of $n_s$ small and
$n_l$ large monomers, with crowders, in a cylinder that is periodic along
its axis.

| Value | Setting |
| --- | --- |
| Wall potential | `fix wall/region` on the cylinder of radius $r$, style `lj126` ([LAMMPS](https://docs.lammps.org/fix_wall_region.html)): $E(d) = 4\epsilon_w[(\sigma_w/d)^{12} - (\sigma_w/d)^6]$ for $d < r_{c,w}$, where $d$ is the distance from a particle to the cylinder surface |
| Wall–small monomer | $\epsilon_w$ = keyword `epss`, $\sigma_w = 1$, $r_{c,w} = 2^{1/6}$ |
| Wall–large monomer | $\epsilon_w$ = keyword `epsl`, $\sigma_w = (1 + a_l)/2$, $r_{c,w} = 2^{1/6}\sigma_w$ |
| Wall–crowder | $\epsilon_w = 1$, $\sigma_w = (1 + a_c)/2$, $r_{c,w} = 2^{1/6}\sigma_w$ |
| `epss`, `epsl` | 5 and 5 in every simulation |
| Bonds | small–small and small–large FENE bonds, each with its own $\sigma_{ij}$ |
| Masses | small monomers 1, large monomers $m_l$ (keyword `ml`), crowders 1 ($a_c = 1$ in every simulation) |

## TransFociCub

`TransFociCub`: the TransFociCyl chain and
crowders in a periodic cube, so there is no wall. The other values are those
of TransFociCyl, with one exception: in the set of simulations with
$a_l \in \{5, 6\}$ from May 2024, every pair is cut at
$r_{c,ij} = \sigma_{ij}$, not at $2^{1/6}\sigma_{ij}$, and its energy is
shifted to zero at the cutoff. Those pairs are still purely repulsive, but
only the part of the repulsion with $r < \sigma_{ij}$ is kept.

## HnsCub and HnsCyl

`HnsCub` and
`HnsCyl`: a ring of $N$ monomers, $n_h$ H-NS
proteins, and crowders, in a periodic cube or in a cylinder that is periodic
along its axis. Each protein is a core of size 1 bonded to two poles. The
particle types are the monomer (1), the H-NS pole (2), the H-NS core (3),
and the crowder (4). Every pair energy is shifted to zero at its cutoff.

| Pair | Setting |
| --- | --- |
| Monomer–pole | $\epsilon$ = 1 in equilibration and 29 in production, $\sigma = 0.535$, $r_c = 0.685 = 1.28\,\sigma$: the one attractive pair |
| Pole–pole | $\epsilon = 1$, $\sigma = 0.178$, $r_c = 0.2$ |
| Pole–core, pole–crowder | no coefficients written, so LAMMPS mixes them by its default geometric rule ([LAMMPS](https://docs.lammps.org/pair_modify.html)): $\epsilon = 1$, $\sigma = \sqrt{0.178\,a_j}$, $r_c = \sqrt{0.2\,r_{c,jj}}$, purely repulsive. A comment in the input says these pairs do not interact. |
| Core–crowder | $\epsilon = 1$, $\sigma = (1 + a_c)/2$, $r_c = 2^{1/6}\sigma$, whatever the name's `epshc`: the generator writes `epshc` into the name and into an input variable the input never reads |
| Other pairs | the shared choices |

| Value | Setting |
| --- | --- |
| Wall (HnsCyl) | style `lj126`, as in TransFociCyl: monomers $\epsilon_w = 5$, $\sigma_w = 1$, $r_{c,w} = 2^{1/6}$; poles $\epsilon_w = 1$, $\sigma_w = 0.535$, $r_{c,w} = 0.685$; cores $\epsilon_w = 1$, $\sigma_w = 1$, $r_{c,w} = 2^{1/6}$; crowders $\epsilon_w = 1$, $\sigma_w = (1 + a_c)/2$, $r_{c,w} = 2^{1/6}\sigma_w$ |
| Bonds | monomer–monomer FENE; pole–core `bond_style harmonic`, $K = 120$, $r_0 = 0.4$ |
| Angles | `angle_style harmonic` ([LAMMPS](https://docs.lammps.org/angle_harmonic.html)), $E = K(\theta - \theta_0)^2$: monomer–monomer–monomer $K$ = (keyword `kbmm`)/2, $\theta_0 = 180°$; pole–core–pole $K = 50$, $\theta_0 = 180°$ |
| Masses | monomers, poles, and cores 1; crowders $a_c^3$ |
| Thermostat | damping time 1, not 10; for crowders multiplied by $m/a_c = a_c^2$, so their friction $m/a_c^2 = a_c$ grows with their size |
| Time steps | equilibration 0.01; production 0.005 |

## References

The LAMMPS documentation, <https://docs.lammps.org/>, for the forms of the
potentials and the meaning of the commands:
[`pair_style lj/cut`](https://docs.lammps.org/pair_lj.html),
[`pair_modify`](https://docs.lammps.org/pair_modify.html) (mixing rules and
the energy shift; the default mixing rule is geometric),
[`bond_style fene`](https://docs.lammps.org/bond_fene.html),
[`bond_style harmonic`](https://docs.lammps.org/bond_harmonic.html),
[`angle_style harmonic`](https://docs.lammps.org/angle_harmonic.html),
[`fix langevin`](https://docs.lammps.org/fix_langevin.html),
[`fix wall/region`](https://docs.lammps.org/fix_wall_region.html), and
[`fix spring`](https://docs.lammps.org/fix_spring.html).
