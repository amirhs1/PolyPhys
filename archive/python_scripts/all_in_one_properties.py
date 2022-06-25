def all_in_one_properties(ens_evg_properties_files,df_name,geometry,filename=None, round_to=4,rename_cols=True):
    ens_avg_properties_df = df_of_properties(ens_evg_properties_files, df_name, geometry, to_file=False,index_col=0)
    if rename_cols:
            ens_avg_properties_df.rename(columns={"vfrc_m": "phi_m", "vfrc_crowd": "phi_c",'rho_crowd':'rho_c'},inplace=True)
    ens_avg_properties_df = ens_avg_properties_df.round(round_to)
    ens_avg_properties_df['fsd_normalized'] = 0.0
    ens_avg_properties_df['gyr_normalized'] = 0.0
    ens_avg_properties_df['rflory_normalized'] = 0.0
    ens_avg_properties_df['phi_c_normalized'] = (ens_avg_properties_df['dmon'] * ens_avg_properties_df['phi_c']) / ens_avg_properties_df['dcrowd'] 
    unique_simulations = [list(input_set) for input_set in list(ens_avg_properties_df.groupby(['nmon','dcyl','dcrowd']).size().index)]
    for nmon, dcyl, dcrowd in unique_simulations:
        condition = (ens_avg_properties_df['nmon'] == nmon) & (ens_avg_properties_df['dcyl'] == dcyl) & (ens_avg_properties_df['dcrowd'] == dcrowd) & (ens_avg_properties_df['phi_c'] == 0)
        fsd = ens_avg_properties_df[condition]['fsd'].values[0]
        gyr = ens_avg_properties_df[condition]['gyr'].values[0]
        rflory = ens_avg_properties_df[condition]['rflory'].values[0]
        cond_normal = (ens_avg_properties_df['nmon'] == nmon) & (ens_avg_properties_df['dcyl'] == dcyl) & (ens_avg_properties_df['dcrowd'] == dcrowd)
        ens_avg_properties_df.loc[cond_normal,'fsd_normalized'] = ens_avg_properties_df.loc[cond_normal,'fsd'] / fsd
        ens_avg_properties_df.loc[cond_normal,'gyr_normalized'] = ens_avg_properties_df.loc[cond_normal,'gyr'] / gyr
        ens_avg_properties_df.loc[cond_normal,'rflory_normalized'] = ens_avg_properties_df.loc[cond_normal,'rflory'] / rflory

    ens_avg_properties_df['phi_c_eff'] = (np.pi * ens_avg_properties_df['dcrowd'] ** 3 / 6) * ens_avg_properties_df['ncrowd'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dcrowd']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['phi_c_eff_normalized'] = (ens_avg_properties_df['dmon'] * ens_avg_properties_df['phi_c_eff']) / ens_avg_properties_df['dcrowd'] 
    ens_avg_properties_df['rho_c_normalized'] = ens_avg_properties_df['rho_c'] * ens_avg_properties_df['dcrowd'] ** 2
    
    ens_avg_properties_df['rho_c_eff'] = ens_avg_properties_df['ncrowd'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dcrowd']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['rho_c_eff_normalized'] = ens_avg_properties_df['rho_c_eff'] * ens_avg_properties_df['dcrowd'] ** 2
    
    ens_avg_properties_df['phi_m_eff'] = (np.pi * ens_avg_properties_df['dmon'] ** 3 / 6) * ens_avg_properties_df['nmon'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dmon']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['rho_m_eff'] = ens_avg_properties_df['nmon'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dmon']) ** 2) * ens_avg_properties_df['lcyl'])
    
    
    ens_avg_properties_df = ens_avg_properties_df.round(round_to)
    if filename != None:
        ens_avg_properties_df.to_csv(filename+'.csv',index=False)
    return ens_avg_properties_df