def log_datasets(
    log: str,
    geometry: str = 'biaxial',
    group: str = 'bug',
    lineage: str = 'segment'
) -> Tuple[str, str]:
    """
    generates the 'details' and 'runtime' 'csv' files with a pre-defined list
    of header names.

    The 'details' file contains the information about  input parameters
    and settings in a set of LAMMPS simulations.

    The 'runtime' file contains information about the runtimes of a set of
     LAMMPS simulations.

    Parameters
    ----------
    log: str
        The name of the LAMMPS log file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    lineage: {'segment', 'whole'}, default 'segment'
        Type of the input file.

    Return
    ------
    details_out: str
        The name of the details dataset.
    runtime_out: str
        The name of the runtime dataset.
    """
    sim_info = SumRule(
        log,
        geometry=geometry,
        group=group,
        lineage=lineage,
    )
    details_out = sim_info.ensemble + '-log_details.csv'
    with open(details_out, 'w') as detailsfile:
        # write simulation information
        for lineage_name in sim_info.genealogy:
            detailsfile.write(f"{lineage_name},")
        for attr_name in sim_info.attributes:
            detailsfile.write(f"{attr_name},")
        # neig_modify delay NUM every NUM check YES/NO:
        detailsfile.write(
            ',ens,run_seg,rskin,delay,every,check,'
        )
        detailsfile.write(
            'total_time_s,cores,timestep,atoms,ts_per_sec,'
        )
        # Section columns: min time, avg time, max time, %varavg, %total"
        # Section rows: Pair, Bond, Neigh, Comm, Output, Modify, Other
        detailsfile.write(
            'pair_avg_s,pair_pct,bond_avg_s,bond_pct,neigh_avg_s,'
            'neigh_pct,comm_avg_s,'
        )
        detailsfile.write(
            'comm_pct,output_avg_s,output_pct,modify_avg_s,'
            'modify_pct,other_avg_s,other_pct,dangerous\n'
        )
    runtime_out = sim_info.ensemble + "-log_runtime.csv"
    with open(runtime_out, 'w') as runfile:
        runfile.write('groupname,filename,n_cores,n_atoms,wall_time\n')
    return details_out, runtime_out
