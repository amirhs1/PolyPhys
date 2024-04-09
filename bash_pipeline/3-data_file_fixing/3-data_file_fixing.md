# Setting up the *PipeLine* packages

**Check the project** and **uncomments relevant lines** in **fix_PairIJ_all_data_files_v---.sh** and **fix_PairIJ_all_data_files_v**.

1. Go to the **___-trjs** directory with this pattern: **N#D#ac#-trjs**
2. Check and (if needed) fix **....all.data** files via running **fix_PairIJ_all_data_files_v---.sh**. This bash script fixes some lines in the **all** data (topology) files, so they can be read by **MDAnalysis** package. This step can be dropped if the bug/issue in **MDAnalysis** is resolved.
**Caution**: the *fix_PairIJ_all_data_files_v#.sh* assumes the pattern of the parent directory is **N#D#ac#-trjs** and the pair style is **lj-cut**.
3. copy **___nc0___ens1.all.data** to **nc0_fake_all.data**.
**Caution**: Some templates are availble in [archive_technique](./archive-technique/) directory. All of these templates start **data_templete** but itis better to use the origin ones from the simulations.
4. Check the **allLine** parameter based on the **project** and **geometry**.
5. Chech the **for** loop based on the **project** and **geometry**.
6. Run the **fake_bug_data.sh** script to fake **N___.bug.data** for **bug trajectories**. Since we want to use the data file to extract the mass and topology information, such a faking is ok.
**Caution**: set the **operating system** for **sed**.
**Caution: this doesn not work on both operating systems**: **fake_bug_data.sh** copy the first 10 lines of the real **all.data** file to the generated **bug.data** file to ensure that the later has the correct information about the simulation box. See the bash script for detailed explanation.

The following step are need if the ablove approach does not work - See the [archive_technique](./archive-technique/) directory:

1. Go to the **___-trjs** directory.
2. Check and (if needed) fix **....all.data** files via running **fix_PairIJ_all_data_files_v---.sh**. This bash script fixes some lines in the **all** data (topology) files, so they can be read by **MDAnalysis** package. This step can be dropped if the bug/issue in **MDAnalysis** is resolved.
3. copy **___nc0___ens1.all.data** to **data_template-biaxial-bug-N___D___L___.data**.
4. Go to the *$HOME/github/sumrule_pipeline*:
    1. Run **git pull git@github.com:amirhs1/polyphys.git master*** to pull from repository.
    2. Enter the **passcode** for **Graham** cluster.
    3. Run **rsync -axvH --no-g --no-p --exclude='.*' PACKAGE_NAME /destination** to copy directory, excluding hidden files.

5. Request a computing node with **salloc --ntasks=1 --cpus-per-task=1 --mem=3G --time=02:00:00 --account=def-someuser**.
6. Active **daskEnv** by running **source ~/daskEnv/bin/activate**
7. Run the **datafile_generator.py** in the **____-trjs** directory. This python script creates a **bug** data (topology) file for bug from bug trajectory file. The positions and velocities of atoms in this data file are not important since they are overwritten by the data from the trajectory file once read inby MDAnalyis; only, the topolgy and box information are important.
