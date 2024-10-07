# Study a heterogeneous ring polymer in Lammps

In a simulation directory, there are 4 files: **initial_config.data**,
**min_init_config.lmp**, **input.lmp**, and **submit.sh**. The data file is
generated by means *Moltemplate* or *amirhsi-defined* approac and is a
hetergenoues ring polymer mostly in a helical conformation. The bash file
contains information about running LAMMPS on a PC in serial or on a HPC
facility in parallel. This bash script calls LAMMPS twice; it calls LAMMPS
first to minimize the initial helical conformation and second to equilibrate
the system and produce thermodynamically equlibrated states. These two calls
are explained in details below.

The bond energies are ussually too highs (due to highly overlapped atoms) in
ininital (usually helical) conformations; as a result, a integration technique
is apply to overcome this problem; for instance, by choosing **fix nve/limit**
integration scheme in *LAMMPS*. For these initital configuration, it is not
recommanded to bond monomers with the FENE potential as if one starts with the
FENE bond potential, they will probably see that LAMMPS crashes due to lost
atoms or too long bonds. For this, it is better to bond the monomers in their
initial (usually helical) confromations with the harmonic potential, relax
these monomers with a Langevin thermosthat, and then dump out a data file of a
relaxed conformation after, say, $10^5$ steps with $dt=0.01\tau$. This first
call of LAMMPS is done by the **minimize_initial_config.lmp** and the output
data file is called **minimized.initial.config.data**.

Next, LAMMPS is called for the second time by **input.lmp**. In this call, the
bond potenital is changed to FENE, crowders added (if needed) and then the
system is equilibarted for another, say, $1\times 10^6$ steps. Once the system
is equilibrated again with **fix nve/limit** and **fix langevin**. Then, the
limit on the **fix nve** is **dropped** and the system is aqulibrated for
another $4\times 10^6$ steps with $dt=0.01\tau$ ($dt=0.02\tau$ does not work)
with **fix nve** and **fix langevin**. After this **two-step** equilibration
phase. The production phase is done by running the system for $10^8$ steps with
$dt=0.005\tau$ ($a_c\le 1.0\sigma$) or $dt=0.001\tau$ ($a_c<1.0$).

**CATUION:** The time step of system at the equilibration phase $dt$ is changed
from $dt=0.01\tau$ to $dt=0.02\tau$ on the date $20220701$

**CAUTION:** The **pre-equilibration** phase for the polymer with **harmonic**
bond potential has been applied for simulations after the date $20220524$.

**CAUTION:** The **two-step** equilibration phase for the crowded system has
been applied for simulations done after the date $20220701$. The rationale for
apply two-step equilibration phase is:

1. The bead-and-spring polymers under study are long and need a long time for
   relaxing their internal degrees of freedom, allowing for an appropriate
   sampling of the phase space and ensuring **ergodicity assummption**.
2. The goal is to study the equilibrium properties of the system, not the
   dynamic ones, thus allowing the system to evolve without any limitions on
   the motion of its constituents (i.e., replacing **fix nve/limit** with
   **fix nve**) does not impact the production phase. Although it is highly
   possible that many initial (after relaxing the **fix nve/limit**)
   configurations of the system are lost and computaional cost increases, but
   configurations in the production phasr are statistically more indepdent and
   less correlated with themselves and with the status of the system at the
   equilibration phase.

## Run file changes: Cylindrical space

- Rule: If a changed applied in a given data (for instance, *A new strategy
  defined for the equilibration phase in the presence of crowders*), then this
  change automatically applied to the simulations run after that day.

|Folder name|describtion|
|:-|:-|
|20220613-ns400nl5al5d20ac1-harmonic_fene|Now data trj files are named based on their "whole" or "segment" naming style in their `pholyphys.manage.parser` module, a single chain first has harmonic bond for smooth relaxing and then swiched to FENE and crowders added|
|20220701-d21nl5ns400al1ac1phic0-new_equlibration_phase_new_submits| A new strategy defined for the equilibration phase in the presence of crowders. The number of time steps in the equilibration phase changed. `def-byha` used on slurm. $dt=0.02\tau$ used for the `fix nve/limit` phase of the equilibration. $dt=0.01\tau$ used for the `fix nve` phase of the equilibration.|
|20220701-d21nl5ns400al3ac1-new_equlibration_phase_new_submits|`rrg-byha` used.|
|20220704-d21nl5ns400al1ac1phic0.1_0.4-used_already_minimized_data-def_byha_used-submit_changes| Instead of separately minimize each helical confomration with `foci-ring-cylinder-init_config_minimize-harmonic.lmp` on the cluster, the data file generated via Moelteplate  first minimized on PC and then comy pased for all the simluation. Date file from `20220701-d21nl5ns400al1ac1phic0-new_equlibration_phase_new_submits` as the input for `foci-ring-cylinder-init_config_minimize-harmonic.lmp` (seed  `i=1` used for langevin thermostat), so the same file used in all the simulations as the input data file. `neighbor`, `neighbor_modify`, and `comm_modify` commands changed from *multi* to *single* since all the atom types have the same size. `def-byha` used. Particle limit for useing 8 core or 16 cores changed from $25000$ to $30000$.|

## Run file changes: Cubic/free space