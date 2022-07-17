# How to compress and archive a file

1. Connect to the HPC data trnasferring node, using your username and password. For Sharcnet's Graham cluster, the node is `username@gra-dta1.computecanada.ca` or `username@gra-dta1.sharcnet.ca`.
2. Use the most recent version of [zip_tar-all_simulations-v1.0.sh](./zip_tar-all_simulations-v1.0.sh) to compress and tar the files and directoriess of interest:
3. Use [archive_rsync-v1.0](./archive_rsync-v1.0.sh) to rsync the directory of interest.

**CAUTION:** Steps 2 and 3 are for the *___-all_simulations* directory, modify accrodingly for other directories.
