# How to restart a Lammps simulation

Follow these instructions to restart a crashed simulation; the main idea is to find the *restart* file that is written at the closest time step to the crashed time step.

1. Detects the simulations that are crashed.
2. Create a new folder for each of them **SIMULATION_NAME_res**.
3. Copy the last **restart** file that match the **all___.lammpstrj.gz** to that folder.
   1. **Caution:** the last time step of a **all___.lammpstrj.gz** file depends on the number of time steps in the equilibrium process; for instance, if equilibrium timestep is 10000 and the sampling time step is 14*50000, then the last time step of each **all___.lammpstrj.gz** file is j*50000+10000 where j is between 1 and 14. This means for a complete j-th **all___.lammpstrj.gz** file, you need j*50000+10000 restart file if you write restart files every 10000 time step.
4. Copy a template **restart-input-___.lmp** file that matches your problem to that folder.
5. Modify the **restart-input-___.lmp** to match this new simulation. Check the following:
    a. The name of the restart file given as an input to **read_restart** command.
    b. The input variables that set the random number generator (fix langevin, etc) and simulation box (change_box, region, fix ID atom_ID wall/region, etc) should be given as input to the restart files. These variables is including but not limited to **r**, **lz**, **epsilon1**, **sig2**, and **randseed**.
    **Caution**: The value of **i** should match the *ens* number of the simulation.
    c. The values of variables that you used to create output files or/and manage the lammps scripts. variables is including but not limited to **n_crowd**, **i**, and **n_bug**.
    **Caution:** Check the style and name of **write_data** is different in various types of simulations.
    d. Set the **neighbor** and **nieghbor_modify**before **read_restart** command.
    e. **no processor** related command is allowed in a restart job.
    f. the **fix wall** command.
    g. the **run** and **timestep** commands.
    h. the **dump**, **dump_modify**, and **undump** commands, espcailly those used in the loop over **dump** and **run**.
    i. check **write_data**.
6. rename **restart-input-___.lmp** to **restart.lmp**.
7. Copy **submit.sh** file to restart directories and edit each submit file individually in the directories
   1. **Caution:**Change the name of lammps input file in the submit script to **restart.lmp**.
8. Run **loop_on_sbatch_restart.sh** to do the restart simulations.
9. Copy **all___.lammpstrj.gz** files, **bug___.lammpstrj**, and **log.lammps** to the **complete** directories, using **copy_files_from_res_example.sh** script as a template.
10. Go to each **complete** directory and do the following steps and run **restart_bug_trj_merger.sh** in each complete folder or run it as a loop.

## How to collect trajectories for the exctraction phase

1. Run **res_new_dir_gen.sh** to create new for the simulations that are restarted.
2. Edit and run **restart_files_collector.sh** trajectory, log, and data files from incomplete and restart directories to the new directories.
   1. **Caution** : edit the *restart_files_collectors.sh* to match the new of the incomplete and restart directories.
3. Edit and run **restart_trj_bug_merger.sh** to merge the *bug* incomplete and restart trajectories, creating full trajectories.
   1. **Caution** : edit the **ens** and **dir** in *restart_trj_merger.sh* to match the new of the incomplete and restart directories.
4. Edit and run **restart_trj_all_merger.sh** to merge the *all* incomplete and restart trajectories, creating full trajectories.
   1. **Caution**: edit the **ens** and **dir** in *restart_trj_merger.sh* to match the new of the incomplete and restart directories.
   2. **caution**: check whether the *all....lammpstrj* files are looped over **j** or not, based on this, edit the name of incomplete/restart all files in the merger).
5. Open the *complete* simulation directories to mannually rename the *all* incomplete and restart trajectories.
