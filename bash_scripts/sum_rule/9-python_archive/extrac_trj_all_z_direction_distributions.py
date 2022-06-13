def extract_all_trj_cog_polymer(all_data, all_trj, geometry):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    z_cutoff = 100 * cellAttrs.dmon
    near_chrm_cog = cell.select_atoms(f'around {z_cutoff} resid 1', updating=True)
    crds = near_chrm_cog.select_atoms('resid 0', updating=True)
    chrm = near_chrm_cog.select_atoms('resid 1', updating=True)

    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    bin_size = 0.1 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    # z direction of the cylindrical coordinate system
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = z_cutoff
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # bin size 5 degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
               
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)
        
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(crds.positions[:,1], crds.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_crd = np.add(thetahists_crd, dummy_hist)
        
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
        
    lastname = 'Crd'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')

def extract_all_trj_box_center(all_data, all_trj, geometry):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    z_cutoff = 100 * cellAttrs.dmon
    near_center = cell.select_atoms(f'prop abs z <= {z_cutoff}', updating=True)
    crds = near_center.select_atoms('resid 0', updating=True)
    chrm = near_center.select_atoms('resid 1', updating=True)

    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    bin_size = 0.1 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    # z direction of the cylindrical coordinate system
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = z_cutoff
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # bin size 5 degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
               
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)
        
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(crds.positions[:,1], crds.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_crd = np.add(thetahists_crd, dummy_hist)
        
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
        
    lastname = 'Crd'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')
    
def extract_all_trj_polymer_cog_fsd(all_data, all_trj, geometry):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    chrm = cell.select_atoms('resid 1')
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # radial direction of the cylindrical coordinate system
    bin_size = 0.1 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    # z direction of the cylindrical coordinate system
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.lcyl
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # in degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
    
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        chrm_fsd =  fsd(chrm.positions)
        z_cutoff = 0.5 * (chrm_fsd+cellAttrs.dmon) # the radii of the three types of particle.
        around_chrm = cell.select_atoms(f'around {z_cutoff} resid 1', updating=True,periodic=False)
        crds = around_chrm.select_atoms('resid 0', updating=True,periodic=False)
        
        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)
        
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(crds.positions[:,1], crds.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_crd = np.add(thetahists_crd, dummy_hist)
        
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
        
    lastname = 'Crd'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')
    
def extract_all_trj_polymer_cog_fsd_eff_dcyl(all_data, all_trj, geometry):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    chrm = cell.select_atoms('resid 1')
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # radial direction of the cylindrical coordinate system
    bin_size = 0.1 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)
    # z direction of the cylindrical coordinate system
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.lcyl
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # in degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
    
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        
        chrm_fsd =  fsd(chrm.positions)
        z_cutoff = 0.5 * (chrm_fsd+cellAttrs.dmon) # the radii of the three types of particle.
        rpos_mon = np.linalg.norm(chrm.positions[:,:2], axis = 1)
        rpos_max = np.max(rpos_mon)
        around_chrm = cell.select_atoms(f'around {z_cutoff} resid 1', updating=True,periodic=False) # particles within the nucleoid in z direction
        crds_cyl = around_chrm.select_atoms('resid 0', updating=True,periodic=False) # crowders within the nucleoid 
        crds = crds_cyl.select_atoms(f'prop abs y <= {rpos_max} and (prop abs x <= {rpos_max})', updating=True,periodic=False) # # crowders within the nucleoid with x,y smaller than effective radius.
        
        # histogram in r direction
        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)
        
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(crds.positions[:,1], crds.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_crd = np.add(thetahists_crd, dummy_hist)
        
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
        
    lastname = 'Crd'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')