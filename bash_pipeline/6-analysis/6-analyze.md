# Analyze probe files to create ensembles and ensemble averages

 1. Copy the **polyphys** direcotry, **setup_gnuparallel.sh**, **analysis-PROJECT.py**, and
    **submit_analysis_parallel.sh** to the probe directort,where **PROJECT**
    is the name of the project.
 2. Check the name **project** and the **input_databases** in the python script.
    In the analysis phase, the name of project in the *analysis-PROJECT.py*
    script has eitheor of these two suffixes: *Segment* if the "bug" or
    "nucleoid" lammps trajectory files are more than one file and *Whole* if it
    is only one file.
 3. Run **setup_gnuparallel.sh** file.
 4. Check the detail of **submit-analysis_parallel.sh** and then submit it.
