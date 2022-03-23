# How to use GNU Parallel to run **extracton** scripts form the **original** Lammps output directories

1. If there are two trajectory directories (**N_D_ac_phic_-trjs** and **N_D_ac_phic_-trjs_cont**), move all the trajectory files in **N_D_ac_phic_-trjs_cont** to **N_D_ac_phic_-trjs** and create as single directory ****N_D_ac_phic_-all_trjs**; otherwise, go to the directory with pattern **N_D_ac_phic_-trjs**.
2. Run **create_trjs_directories.sh** to create the simulations directories and move the corresponding toplogy and trjaectory files into them, and to copy relevant scripts and files for doing the **probe** phase in python.
3. Run **probe_runfiles-first_time.sh** to copy run files from **PolyPhys** to **N___-all_trjs** directory.
    1. *Caution: check how many **probe_all_in_one_---.py** you have.*
4. Check the **PipeLine.py** file to make sure which packages you need.
5. Add the packages to the ***submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
6. Set the number of cores, time, and memeory in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file based on what group/type of particles you want to analyze.
7. Check the name of Python file in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
8. Run the **sbatch-probe_all_in_one.sh**, **sbatch-probe_all_trj_segments.sh**, or **sbatch-probe_all_bug_segments.sh** file.
9. After the rum finsihed, check slurm file and change its name to **N_D_ac_phic_-probe-slurm_report.out**
10. Run **create_probe_directories.sh** to create a new directory **N_D_ac_phic_-probe** and move all the **csv**, **npy**, or **txt** files to their corresponding simulation directories.
