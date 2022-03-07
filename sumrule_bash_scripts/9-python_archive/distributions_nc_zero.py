def distributions_nc_zero(old_csvs, save_to, old_name='Mon',new_name='Crd',empty=True):
    """
    distributions_nc_zero creates fake distribution for crowders/monmers/species of interest in which all the frequencies can be zero, based on the input distributions in a system with no crowders. 
    
    Caution:
    For systems with no crowders, local distributions have been already created when the bug/polymer/chain is analyzed in "bug"-type analysis. For this reasons, the distribtuions of monomers and crowders in "all"-type analysis of these class of systems (sysstems with nc=0) are faked. For monomers, the old disstributions are just rename based on the standards of "all"-type analysis. For crowders, empty distributions (distirbutions in which all the frequencies are zero) the old distributions are just rename based on the standards of "all"-type analysis. 
      
    Parameters:
    old_csvs (str): address of the monomer distributions.
    save_to (str): address otv which the new files are saved.
    old_name (str): the species or sub-system name of the old distributions.
    new_name (str): the species or sub-system name of the new distributions.
    empty (bool): If True, the fake distributions are set to zero; otherwise, just save the old distributions with its new names 
    
    Return:
    Save crowder dstributions to the disk.
    
    Requirements:
    Pandas
    """  
    for old_csv in old_csvs:
        df_name = new_name.join(old_csv.split('/')[-1].split(old_name))
        if empty:
            df = empty_dataframe(old_csv)
        else:
            df = pd.read_csv(old_csv,header=0,index_col=0)
        df.to_csv(save_to+df_name)