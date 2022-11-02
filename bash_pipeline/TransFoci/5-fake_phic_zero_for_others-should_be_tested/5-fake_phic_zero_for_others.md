# Copy simulation with phic=0 to this simulation set

**Caution**: Sometimes instead of **ns#nl#al#D#ac#-___** pattern, you have **ns#nl#al#D#ac#phic#-___** pattern.

1. Go the simulation **ns#nl#al#D#ac#-probe** directory that has the **probe** data of **ncrowd=0 (or phi_c=0)**.
2. Create a **test** directory and copy the 8 simulation directories of **ncrowd=0** to that folder.
3. Move the **test** directory to the **ns#nl#al#D#ac#-probe** that **does not** have **the phic=0 simulation data** and you want to fake them for it (**Caution:** do **not** copy the **trj and data files**).
4. Copy **fix_phic0_sims.sh** the test directory:
    1. Choose the version **sed** based on the operating system.
    2. **INSTALL** *rename* and *gawk* on **MacOSX*8 by **homebrew**.
    3. Check that the directory and file names **patterns** and **extensions** are correct for this **PROJECT**.
5. Run the **fix_phic0_sims.sh**
6. Move the faked phic=0 directories to **ns#nl#al#D#ac#-probe** that does not have the phic=0 simulation data.
