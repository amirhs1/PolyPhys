# How to use GNU Parallel to run **probe** scripts in the **trj**   directory

**Caution**: Sometimes instead of **N#D#ac#-___** pattern, you have **N#D#ac#phic#-___** pattern.

1. If there are two trajectory directories (**N#D#ac#-trjs** and **N#D#ac#-trjs_cont**), move all the trajectory files in **N#D#ac#-trjs_cont** to **N#D#ac#-trjs** and create as single directory ****N#D#ac#-all_trjs**; otherwise, go to the directory with pattern **N#D#ac#-trjs**.
2. Run **probe_runfiles.sh** to copy run files from **PolyPhys** to **N___-all_trjs** directory.
    1. *Caution: check how many **probe_all_in_one_---.py** you have.*
3. Run **probe_gnuparallel-___.sh** to copy the **polyphys** package and **probe___.py** file to each **simulation directory**.
4. Check the **PolyPhys** package to make sure which packages you need.
5. Add the packages to the ***submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
6. Set the number of cores, time, and memeory in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file based on what group/type of particles you want to analyze.
7. Check the name of Python file in the **submit-probe_all_in_one.sh**, **submit-probe_all_trj_segments.sh**, or **submit-probe_all_bug_segments.sh** file.
**Cuation:** Check the name of **python**, its **report**, the **pattern** of **directories** and the **report** of the **gnuparallel** code.
8. Run the **sbatch-probe_all_in_one.sh**, **sbatch-probe_all_trj_segments.sh**, or **sbatch-probe_all_bug_segments.sh** file.
9. After the rum finsihed, check slurm file and change its name to **N#D#ac#-probe-slurm_report-probe_all_in_one.out**, **N#D#ac#-probe-slurm_report-probe_all_trj_segments.out**, **N#D#ac#-probe-slurm_report-probe_all_bug_segments.out**, or the like.
**TODO:** this step can be combined with the next one in as if-statement or switch statement.
10. Run **create_probe_directories.sh** to create a new directory **N#D#ac#-probe** and move all the **csv**, **npy**, or **txt** files to their corresponding simulation directories.