# How to compress and archive a file

1. Connect to the HPC data trnasferring node, using your username and password. For Sharcnet's Graham cluster, the node is `username@gra-dta1.computecanada.ca` or `username@gra-dta1.sharcnet.ca`.
2. Use the most recent version of **arch\_---\_v---.sh** to compress and tar the files and folders of interest:
    1. There are **four different archiver** bash scripts for all_simulations, extraction_bug, analyze_bug, and trjs_bug directories.
    2. There is a **loop_on_arch_simulations.sh** to run tar on all the **all_simulations** directory together.
3. Move the files/directories to **project** folder, using the most recent version of **archive\_rsync\_v---.sh**.
    1. Check files are copied completely and then remove the ones in *scratch* directory.