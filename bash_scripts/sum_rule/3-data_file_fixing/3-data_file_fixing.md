# Setting up the *PipeLine* packages

1. Go to the **___-trjs** directory.
2. Check and (if needed) fix **....all.data** files via running **fix_PairIJ_all_data_files_v---.sh**. This bash script fixes some lines in the **all** data (topology) files, so they can be read by **MDAnalysis** package. This step can be dropped if the bug/issue in **MDAnalysis** is resolved.
3. copy **___nc0___ens1.all.data** to **fake_nc0_all.data**.
4. Run the **fake_bug_data.sh** script to fake **N___.bug.data** for **bug trajectories**. Since we want to use the data file to extract the mass and topology information, such a faking is ok.

The following step are need if the ablove approach does not work:

1. Go to the **___-trjs** directory.
2. Check and (if needed) fix **....all.data** files via running **fix_PairIJ_all_data_files_v---.sh**. This bash script fixes some lines in the **all** data (topology) files, so they can be read by **MDAnalysis** package. This step can be dropped if the bug/issue in **MDAnalysis** is resolved.
3. copy **___nc0___ens1.all.data** to **data_template-biaxial-bug-N___D___L___.data**.
4. Go to the *$HOME/github/sumrule_pipeline*:
    1. Run **git pull git@github.com:amirhs1/sumrule_pipeline.git main** to pull from repository.
    2. Enter the **passcode** for **Graham** cluster.
    3. Run **rsync -axvH --no-g --no-p --exclude='.*' PACKAGE_NAME /destination** to copy directory, excluding hidden files.

5. Request a computing node with **salloc --ntasks=1 --cpus-per-task=1 --mem=2G --time=0:30:0 --account=def-someuser**.
6. Active **daskEnv** by running **source ~/daskEnv/bin/activate**
7. Run the **datafile_generator.py** in the **____-trjs** directory. This python script creates a **bug** data (topology) file for bug from bug trajectory file. The positions and velocities of atoms in this data file are not important since they are overwritten by the data from the trajectory file once read inby MDAnalyis; only, the topolgy and box information are important.
