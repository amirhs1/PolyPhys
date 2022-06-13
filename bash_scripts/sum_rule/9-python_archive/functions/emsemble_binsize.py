def ensemble_binsize(simulations, property_name, geometry, single=True, to_file=True, sep='_', **kwargs):
    """
    ensemble generates an ensemble (dataframe) of simulations (columns) for the given property_name in the geometry of interest. 
    
    Caution:
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    By default, the indexes of an ensemble (dataframe) are timesteps. For some properties such as histograms, however, the indexes are bin centers. As a result, simulations are tuples where each tuple can have one or more members. For histogram-like properties, each tuple has two members where the second member is the bin edges.
    
    Parameters:
    simulations (list of tuples): a sorted liste of tuples where each tuple at least has one member (the path to the simulation data for the property_name). For histogram-like proeprties, the second member of each tuple is the path to the bin_edges.
    property_name (str): name of the property measured in each of the simulations.
    geometry (str): shape of the simulation box
    single (bool): whether each tuple is a single file or more.
    to_file (bool): whether save to file or not.
    sep (str): a delimiter that used to extract a simulation name from a simulation path.
    
    Return:
    ensembles (dict): a dictionary of ensembles where keys are ensemble names and values are ensembles (dataframes). In each ensemble (value of dict/dataframe), the columns are the simulations of that ensemble.
    
    Requirements:
    Pandas, Numpy, my own cellAttributes class.
    """
    # Unique ensembles with simailar initial paramters    
    ens_names = [PipeLine.cellAttributes(simulation[0],geometry,warning=False).filename.split('ens')[0] for simulation in  simulations] # name of simulaion groups
    ensembles = dict.fromkeys(ens_names)
    if not single :
        bin_edges = dict.fromkeys(ens_names) # bin edges create bin centers which in turn used as indexes of histograms dataframes.
        
    for ens_name in ensembles.keys():
        ens_simulations = []
        for sim_tuple in simulations: #files are organized based on their name
            sim_name = sim_tuple[0].split("/")[-1].split(sep+property_name)[0]
            if sim_name.split('ens')[0] == ens_name:
                simulation = pd.read_csv(sim_tuple[0],names=[sim_name],**kwargs)
                ens_simulations.append(simulation)    
        ensembles[ens_name] = ens_simulations # each key (ensemble) has a list of dataframes (simulations)
        ensembles[ens_name] = pd.concat(ensembles[ens_name],axis=1) # convert the list of dataframes to columns 
        if not single :
            bin_edges[ens_name] = np.around(np.loadtxt(sim_tuple[1]),decimals=2) # warning bin edge rounded based on the binsize variable in extract trj functions
            print("Warning: bin edges, bin centers, and the like are rounded to 2 decimals.")
            bin_centers = np.around((bin_edges[ens_name][:-1] + bin_edges[ens_name][1:]) / 2.0,decimals=2)
            ens_index = pd.Index(bin_centers, dtype = float) # use bin_centers as index
            ensembles[ens_name].set_index(ens_index,inplace=True)
        if to_file:
            cell_attrs = PipeLine.cellAttributes(ens_name,geometry,warning=False)
            pattern = re.compile('([a-zA-Z\_]+)')
            words = pattern.split(cell_attrs.filename)
            bin_factor = float(words[words.index('binFactor')+1]) # crowders size/diameter
            output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}nc{cell_attrs.ncrowd}binFactor{bin_factor}-{property_name}.csv'
            ensembles[ens_name].to_csv(output)
    return ensembles 