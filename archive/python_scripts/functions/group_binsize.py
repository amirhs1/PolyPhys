def group_binsize(ensembles, property_name, geometry, to_file=True):
    """
    ens_avg_group averages over the property_name in all the simulations of each ensemble in the ensembles.
    
    Caution: 
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group. 
    
    Parameters:
    ensembles (dict): a dict of ensembles (dataframes) where the number/names of columns in each ensemble is euqal ot the number/names of simulations in that ensemble. The values in each column are the measurements of a given property in each simulation.
    property_name (str): name of the property measured in each of the simulations.
    geometry (str): shape of the simulation box.
    to_file (bool): whether save to file as csv or not
    
    Return:
    ens_avgs (dict): a dict of ensemble-averaged groups where keys are names of ensemble-averaged groups and values are dataframes of ensemble-averaged
    
    Requirements:
    Pandas, Numpy, PipeLine
    
    """
    # Averging over ensembles with simailar initial paramters
    ens_avgs = {} # ensemble-average groups: each group is a bunch of ensemble-averaged simulations.
    for ens_name, ensemble in ensembles.items():
        ens_avg = np.zeros(len(ensemble)) # the length is the number of measurements on the property which is the same in all the simulations in an ensemble.
        ens_count = 0
        for sim_name in ensemble.columns:
            if sim_name.split('ens')[0] == ens_name:
                ens_count += 1
                ens_avg += ensemble[sim_name].values
        ens_avg = ens_avg / ens_count # ensemble-averaged property
        ens_avg = pd.DataFrame(np.array(ens_avg),columns=[f"{ens_name}-ens_avg"]) # converting to a dataframe with ens_avg as values of a column and the ens_name as the name of that column.
        ens_avg.set_index(ensemble.index.values,inplace=True)
        if to_file:
            cell_attrs = PipeLine.cellAttributes(ensemble.columns[0],geometry,warning=False)
            pattern = re.compile('([a-zA-Z\_]+)')
            words = pattern.split(cell_attrs.filename)
            bin_factor = float(words[words.index('binFactor')+1]) # crowders size/diameter
            output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}nc{cell_attrs.ncrowd}binFactor{bin_factor}-{property_name}-ens_avg.csv'
            ens_avg.to_csv(output)
        ens_avgs[ens_name] = ens_avg
    return ens_avgs

