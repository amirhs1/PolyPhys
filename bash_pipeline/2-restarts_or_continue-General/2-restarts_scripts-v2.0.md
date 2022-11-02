# How to restart or continue a simulation in LAMMPS

## Step 1: Restarting a simulation: See *restart_simulations* directory

Follow these instructions to restart a crashed simulation; the main idea is to find the *restart* file that is written at the closest time step to the crashed time step or the starting time step of the crashed loop.

1. Detects the simulations that are crashed.
2. Create a new folder for each of them **SIMULATION_NAME_res** (See the **res_dir_generator.sh** to see what you need to do).
3. Copy the last **restart** file that match the **all___.lammpstrj.gz** to that folder.
   1. **Caution**: the last time step of a **all___.lammpstrj.gz** file depends on the number of time steps in the equilibrium process; for instance, if equilibrium timestep is 10000 and the sampling time step is 14*50000, then the last time step of each **all___.lammpstrj.gz** file is j*50000+10000 where j is between 1 and 14. This means for a complete j-th **all___.lammpstrj.gz** file, you need j*50000+10000 restart file if you write restart files every 10000 time step.
4. See a example **restart-input-___.lmp** file and the **LAMMPS** webstie to see how to modify an old input and creat a restart input file.
5. Modify the **restart-input-___.lmp** to match this new simulation. Check the following:
   1. the name of the restart file given as an input to **read_restart** command.
   2. The input variables that set the random number generator (fix langevin, etc) and simulation box (change_box, region, fx ID atom_ID wall/region, etc) should be given as input to the restart files. For example in the **SumeRuleCyl** project, follow these steps:
      1. The input variables are including but not limited to **r**, **lz**, **epsilon1**, **sig2**, and **randseed**.
      **Caution**: the value of **i** should match the **ens** number of the simulation.
      2. The values of variables that you used to create output files or/and manage the lammps scripts. These variables include but are not limited to **n_crowd**, **i**, and **n_bug**.
      **Caution:** Check the style and name of **write_data** is different in various types of simulations.
      3. Set the **neighbor** and **nieghbor_modify**.
      4. **no processor** related command is allowed in a restart job.
      5. the **fix wall** command.
      6. the **run** command.
      7. the **dump**, **dump_modify**, and **undump** commands, espcailly those used in the loop over **dump** and **run**.
      8. check **write_data**.
         **Caution: Use the appropriate retart-input-....lmp files; see the explanation in each file.**
         **Caution: Use the appropriate esm_restart_....sh files; see the explanation in each file.**
6. rename **restart-input-___.lmp** to **input.lmp**.
7. Copy **submit.sh** file to restart directories and edit each submit file individually in the directories
   1. **Caution:**Change the name of lammps input file in the submit script to **restart.lmp**.
8. Run **loop_on_sbatch_restart.sh** to do the restart simulations.
9. Copy **all___.lammpstrj.gz** files, check **j** value, **bug___.lammpstrj**, and **log.lammps** to the **complete** directories, using **copy_files_from_res_example.sh** script as a template.10. Go to each **complete** directory and do the following steps and run **restart_bug_trj_merger.sh** in each complete folder or run it as a loop.

## Step 2: Generating a complete simulation directory: See *merge_outputs* directory

1. Run **org_dir_generator.sh** to rename *incomplete* and *res* directories, and create a *complete* directory.
   1. **Caution:** The name of append name for *incomplete* and *restart* directories are different in different projects.
2. Edit and run **copy_files_from_---_example** (where **---** is either **res** or **incomplete** keywords.) to copy trajectory, log, and data files from incomplete and restart directories to the new directories. (**Caution** : edit the **restart_files_collectors.sh** to match the new of the incomplete and restart directories).
(**Caution** : check the name of **bug**, **log**, and **all** files since they should be different for **incomplete** and **restart** files.).
3. Edit and run **restart_trj_bug_merger.sh** to merge the **bug** incomplete and restart trajectories, creating full trajectories.
   1. **Caution:** Edit the **ens** and **dir** in **restart_trj_merger.sh** to match the new of the incomplete and restart directories.
4. Edit and run **restart_trj_all_merger.sh** to merge the **all** incomplete and restart trajectories, creating full trajectories.
   1. **Caution:** edit the **ens** and **dir** in **restart_trj_merger.sh** to match the new of the incomplete and restart directories.
   2. **caution:** check whether the **all....lammpstrj** files are looped over **j** or not, based on this, edit the name of incomplete/restart all files in the merger.
5. Open the **complete** simulation directories to mannually rename the **all** incomplete and restart trajectories.

## Step 3: Continuing a finished simulation

1. The steps are similar to restarting a simulation; however, here is important to double check that the file names (the bug, all, data, log, and restart filenames) do not coincide with the old ones.
2. It is sometimes needed to modify the name of a file after the continued simulation is finished. In these cases, use the scripts in **cont_scripts** directory.
