def extract_all_mon_trj(all_data, all_trj, geom):
    """
    The function does all the analysis on the a whole trajectory (simulation),

    Parameters:
    fname (string): the name of trajectory or simulation.
    geom (string): the name of the geometry of the simulation box.

    Requirements:
    All the above-defined fucntions and classes, MDAnalysis.
    
    Caution:
    The histograms are for the cylindrical coordinate system in the geom=cylindrical.
    For other goemteries you need to redfined them.
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry=geom, splitter='lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    TIME_UNIT = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    LJ_NSTEP = cellAttrs.bdump # Sampling steps via dump command in Lammps
    LJ_DT = cellAttrs.dt # Timestep during the sampling process in Lammps
    SAMPLE_DT = LJ_NSTEP * LJ_DT * TIME_UNIT

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = SAMPLE_DT)
    chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain, resid 0 is crowders


    print ("Caution:The histograms are for the cylindrical coordinate system in the geom=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.1
    lmax = cellAttrs.dcyl / 2.
    edge_fname = 'rEdgesMon'
    redges , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_fname)

    # z direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.2
    lmax = cellAttrs.lcyl / 2.
    edge_fname = 'zEdgesMon'
    zedges , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_fname)

    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/24) # in degrees
    edge_fname = 'thetaEdgesMon'
    thetaedges , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_fname)
    
    if rhists_mon.any() != 0:
        raise Exception("The histogram collectors are not empty!")
    if zhists_mon.any() != 0:
        raise Exception("The histogram collectors are not empty!")
    if thetahists_mon.any() != 0:
        raise Exception("The histogram collectors are not empty!")

    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta     
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')