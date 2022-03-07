# How to use GNU Parallel to run **extracton** scripts form the **original** Lammps output directories

1. If there are two trajectory directories (**N_D_ac_phic_-trjs** and **N_D_ac_phic_-trjs_cont**), move all the trajectory files in **N_D_ac_phic_-trjs_cont** to **N_D_ac_phic_-trjs**; else directly go to the directory with pattern **N_D_ac_phic_-trjs**.
2. Check the **PolyPhys** package to make sure which packages you need.
3. Run **trjs_directories.sh** to create the simulations directories and move the files into them.
4. Add the packages to the ***submit-probe_all_in_one.sh** file.
5. Set the number of cores, time, and memeory in the **submit-probe_all_in_one.sh** file.
6. Check the name of Python file in the **submit-probe_all_in_one.sh** file.
7. Run the **sbatch-probe_all_in_one.sh** file.
8. After the rum finsihed, check slurm file and change its name to **N_D_ac_phic_-probe-report.out**
9. Run **probe_directories.sh** to create a new directory **N_D_ac_phic_-probe** and move all the **csv** and **txt** files to that directory.