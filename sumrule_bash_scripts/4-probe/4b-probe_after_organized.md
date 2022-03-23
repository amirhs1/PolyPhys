# How to use GNU Parallel for *probe* scripts in the organized *trjs_bug* and *trjs_all* directories

1. If there are two trajectory directories (**N_D_ac_phic_-trjs** and **N_D_ac_phic_-trjs_cont**), move all the trajectory files in **N_D_ac_phic_-trjs_cont** to **N_D_ac_phic_-trjs** and create as single directory ****N_D_ac_phic_-all_trjs**; otherwise, go to the directory with pattern **N_D_ac_phic_-trjs**.
2. Run **probe_runfiles-after.sh** to copy run files from **PolyPhys** to **N___-all_trjs** directory.
    1. *Caution: check how many **probe_all_in_one_---.py** you have.*
3. Check the **PipeLine.py** file to make sure which packages you need.
4. Add the packages to the ***submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
5. Set the number of cores, time, and memeory in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file based on what group/type of particles you want to analyze.
6. Check the name of Python file in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
7. Run the **sbatch-probe_all_in_one.sh**, **sbatch-probe_all_trj_segments.sh**, or **sbatch-probe_all_bug_segments.sh** file.
8. After the rum finsihed, check slurm file and change its name to **N_D_ac_phic_-probe-slurm_report.out**
9. Run **create_probe_directories.sh** to create a new directory **N_D_ac_phic_-probe** and move all the **csv**, **npy**, or **txt** files to their corresponding simulation directories.