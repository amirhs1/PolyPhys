#cdf Copy simulation with phic=0 to this simulation set

**Caution:** Different projects have differant **PACE_PROJECT_PATTERN** project patterns. Below, **PACE_PROJECT_PATTERN=ns#nl#al#D#ac#** which is for the **TransFociCyl** project.

**Caution:** Sometimes instead of **SPACE_PROJECT_PATTERN-___** pattern, you
have **SPACE_PROJECT_PATTERNphic#-___** pattern.

1. Go the simulation **N#D#ac#-probe** directory that has the **probe** data of **ncrowd=0 (or phi_c=0)**.
2. Create a **test** directory and copy the 8 simulation directories of **ncrowd=0** to that folder.
3. Move the **test** directory to the **N#D#ac#-probe** that **does not** have **the phic=0 simulation data** and you want to fake them for it (**Caution:** do **not** copy the **trj and data files**).
4. Copy **fix_phic0_sims.sh** the test directory:
    1. Choose the version **sed** based on the operating system.
    2. **INSTALL** *rename* and *gawk* on **MacOSX*8 by **homebrew**.
    3. Check that the directory and file names **patterns** and **extensions** are correct for this **PROJECT**.
5. Run the **fix_phic0_sims.sh**
6. Move the faked phic=0 directories to **PACE_PROJECT_PATTERN-probe** that does not have the phic=0 simulation data.
