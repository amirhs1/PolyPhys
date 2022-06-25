# How to extract different information from simulations

 1. Set the name of the simulations directory based on either of these two conventions:

     1. ``#ns#nl#al#D#ac#-all_simulations``
     2. ``#ns#nl#al#D#ac#-cont_simulations``

Here *D* is the size (diamter) of cylinderical simulation box (an open-end simulation box), *al* is the size of large monomers, *nl* is the number of large monomers, *ns* is the number of small monomers, and *ac* is the size of crowders. The volume fraction of crowders *phic* is needed if the simulation direcotry is not a *full* simulation group, i.e. the simulation group does not coverge the full range of *phic* values in the direoctory.

 1. Run **sim_reporter_v---.sh** to to create a report for this siumlation set.
 2. Run **trj_extractor_v1_mutiple_all_single_bug.sh** to extract **bug.lammpstrj**, **all.lammpstrj**,  **all.data**, and **log.lammps** from simulation directories, and move these files to a directory with this name pattern: **N#D#ac\$phic#-trjs**. **Additionally** this command make **run_files-N#D#ac\$phic#** direcotry and move all the Lammps files except the **N#D#ac$phic#-summary.txt** report file to this directory. This script also copy **probe_runfiles-first_time.sh** from the **4-probe** directory to **trjdir**.:wq

**Caution:** The same porcess should be done for **cont_simulations**; however, in this case, the **trj_extractor_v#...** are replace with **trj_extractor_v#_cont...**.

**Cuation** Some scripts used in the **sum_rule** project are deleted in the **trans_foci** project.
