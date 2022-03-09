# How to extract different information from simulations

 1. Set the name of the simulations directory based on either of these two conventions:

     1. ``N*D*ac*phic*-all_simulations``
     2. ``N*D*ac*phic*-cont_simulations``

Here *N* is the number of monomers in a chain, *D* is the diamter of cylinderical simulation box (an open-end simulation box), and *ac* is the size of crowders. The volume fraction of crowders *phic* is needed if the simulation direcotry is not a *full* simulation group, i.e. the simulation group does not coverge the full range of *phic* values in the direoctory.

 1. Run **sim_reporter_v---.sh** to to create a report for this siumlation set.
 2. Depending on the type of all files (one all file or multiple all files), run **trj_extractor_v\$_mutiple_all_files** or **trj_extractor_v\$_one_all_file** to extract **bug.lammpstrj**, **all.lammpstrj**,  **all.data**, and **log.lammps** from simulation directories, and move these files to a directory with this name pattern: **N#D#ac\$phic#-trjs**. **Additionally** this command make **run_files-N#D#ac\$phic#** direcotry and move all the Lammps files except the **N#D#ac$phic#-summary.txt** report file to this directory.

**Caution:** The same porcess should be done for **cont_simulations**; however, in this case, the **trj_extractor_v\$...** are replace with **trj_extractor_v\$_cont...**.
