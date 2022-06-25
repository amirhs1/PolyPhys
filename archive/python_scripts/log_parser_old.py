def log_parser(
    logs: List[Tuple[str]],
    details_out: str,
    runtime_out: str,
    geometry: str = 'biaxial',
    group: str = 'bug',
    lineage: str = 'segment'
) -> None:
    """
    parses a LAMMPS `logs`, gather the information about simulations from
    their `logs`, and write down simulation details to `details_out` csv
    files and runtime information to `runtime_out` csv file.

    Parameters
    ----------
    logs: list of str
        The list of LAMMPS log files where each element of list is a tuple
        file with one element; a string that is the filepath for the log.
    details_out: str
        The name of the 'details' file in which the simulation details is
        written.
    runtime_out: str
        The name of the 'runtime' file in which the runtime information is
        written.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    lineage: {'segment', 'whole'}, default 'segment'
        Type of the input file.
    """
    for log in logs:
        sim_info = SumRule(
            log[0],
            geometry=geometry,
            group=group,
            lineage=lineage,
        )
        filename = sim_info.filename
        ens = sim_info.ensemble_id
        groupname = sim_info.lineage_name
        with open(details_out, mode='w') as detailsfile:
            # write simulation details
            for lineage_name in sim_info.genealogy:
                attr_value = getattr(sim_info, lineage_name)
                detailsfile.write(f"{attr_value},")
            for attr_name in sim_info.attributes:
                attr_value = getattr(sim_info, attr_name)
                detailsfile.write(f"{attr_value},")
        with open(log[0], 'r') as logfile,\
            open(details_out, 'a') as detailsfile,\
                open(runtime_out, 'a') as runfile:
            line = logfile.readline()
            # neigh_modify delay every check page one
            j = 1
            while line:
                if line.startswith('neighbor'):
                    words = line.split()
                    rskin = words[1].strip()  # rskin
                # neigh_modify delay every check page one
                if line.startswith('neigh_modify'):
                    words = line.split()
                    # picking the NUMs and Yes/No from neigh_modify
                    delay = words[2].strip()
                    every = words[4].strip()
                    check = words[6].strip()
                if line.startswith('Loop time'):
                    detailsfile.write(str(j))  # total time
                    detailsfile.write(",")
                    j += 1
                    # neighbor and neigh_modify occurs
                    # one time but other occurs 15 times.
                    detailsfile.write(rskin)  # rskin
                    detailsfile.write(",")
                    detailsfile.write(delay)  # delay
                    detailsfile.write(",")
                    detailsfile.write(every)  # every
                    detailsfile.write(",")
                    detailsfile.write(check)  # check
                    detailsfile.write(",")
                    words = line.split()
                    detailsfile.write(words[3].strip())  # total time
                    detailsfile.write(",")
                    n_cores = words[5].strip()
                    detailsfile.write(n_cores)  # # of cores
                    detailsfile.write(",")
                    detailsfile.write(words[8].strip())  # total timesteps
                    detailsfile.write(",")
                    n_atoms = words[11].strip()
                    detailsfile.write(n_atoms)  # total atoms
                    detailsfile.write(",")
                if line.startswith('Performance:'):
                    words = line.split()
                    detailsfile.write(words[3].strip())  # timesteps per second
                    detailsfile.write(",")
                if line.startswith('Section'):
                    _ = logfile.readline()
                    for i in range(6):  \
                            # Section rows: Pair, Bond, Neigh, Comm, Output,
                        # Modify, Other
                        # Section columns: min time, avg time, max time,
                        # %varavg, %total"
                        line = logfile.readline()
                        sect_min = line.split('|')[2].strip()
                        detailsfile.write(sect_min)
                        detailsfile.write(",")

                        sect_pct = line.split()[-1]  # Pair pct of total time
                        detailsfile.write(sect_pct)
                        detailsfile.write(",")
                    line = logfile.readline()
                    sect_min = line.split('|')[2].strip()
                    detailsfile.write(sect_min)
                    detailsfile.write(",")
                    sect_pct = line.split()[-1]  # Pair pct of total time
                    detailsfile.write(sect_pct)
                    detailsfile.write(",")
                if line.startswith('Dangerous'):
                    words = line.split()
                    detailsfile.write(str(int(words[-1])))  \
                        # # number of dangerous builds
                    detailsfile.write("\n")
                # runtime files
                if line.startswith('Total wall time'):
                    runfile.write(groupname)
                    runfile.write(",")
                    runfile.write(filename)
                    runfile.write(",")
                    runfile.write(n_cores)
                    runfile.write(",")
                    runfile.write(n_atoms)
                    runfile.write(",")
                    words = line.split()
                    runfile.write(words[-1])  # total wall time
                    runfile.write("\n")
                line = logfile.readline()