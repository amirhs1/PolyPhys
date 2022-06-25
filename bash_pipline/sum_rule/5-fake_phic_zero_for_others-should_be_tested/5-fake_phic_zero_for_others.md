# Copy simulation with phic=0 to this simulation set

**Caution**: Sometimes instead of **N#D#ac#-___** pattern, you have **N#D#ac#phic#-___** pattern.

1. Go the simulation **N#D#ac#-probe** directory that has the **probe** data of **ncrowd=0 (or phi_c=0)**.
2. Create a **test** directory and copy the 8 simulation directories of **ncrowd=0** to that folder.
3. Move the **test** directory to the **N#D#ac#-probe** that **does not** have **the phic=0 simulation data** and you want to fake them for it (**Caution:** do **not** copy the **trj and data files**).
4. Copy **phic_zero_rename_dir_and_files.sh** and **fix_stamps_files.sh** to the test directory:
    1. Choose the version **sed** based on the operating system.
5. Run the above mentioned files in **the following order**: **phic_zero_rename_dir_and_files.sh** and **fix_stamps_files.sh**.
6. Move the faked phic=0 directories to **N#D#ac#-probe** that does not have the phic=0 simulation data.
