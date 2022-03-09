# How to use GNU Parallel for *probe* scripts in the organized *trjs_bug* and *trjs_all* directories

1. Move all the files in **N___-trjs_bug** to their **corresponding direcotry** in the **N___-trjs_all**.
2. Run **probe_runfiles.sh** to copy run files from **polyphys** to **N___-trjs_all** directory.
    1. *Caution: check how many **probe_all_in_one_---.py** you have.*
3. Check the **PipeLine.py** file to make sure which packages you need.
4. Add the packages to **submit-probe_all_in_one.sh** files.
5. Run **probe__gnuparallel.sh** to setup the gnuparallel files and directories in simulation directories.
6. Set the number of cores, time, and memeory, **python module** in the **submit-probe_all_in_one.sh** files.
7. Check the name of *Python file* in the **submit-probe_all_in_one.sh** files.
8. Run the **sbatch-probe_all_in_one.sh** file.
9. Check slurm file and change its name to **N___-probe_report.out**
    1. *Caution: Use the slurm files to estimate runtime in future.*
