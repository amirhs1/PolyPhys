# Analyze probe files to create ensembles and ensemble averages

 1. Replace **probe** extension in the names of **probe** directories with **bugWhole** if these directories contains a single Bug file; otherwise, replace it with **bugSegment**.
 2. Check the name **project** and the **input_databases** in the python script.
 3. Copy **polyphys** to the **analyze** directory.
 4. Submit **submit-analyze_v##.sh**

**CAUTION:** For running on a cluster, do the following:

1. Request a computing node with **salloc --ntasks=1 --cpus-per-task=4 --mem=8G --time=1:0:0 --account=def-someuser**.
2. Active **dasnEnv**.
