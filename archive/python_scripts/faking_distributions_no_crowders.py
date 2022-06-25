hist_names = ['rHists','rPhis','rRhos','thetaHists','zHists','zPhis','zRhos']
geometry = 'cylindrical'
## Faking distributions with no crowders
parent_path ='/Users/amirhsi_mini/'
database = parent_path+'trjs_analysis_bug/'
bug_groups = glob(database+"N*")
for bug_group in bug_groups:
    csv_files = glob(bug_group+'/N*nc0*.csv')
    bug_name = 'nc0'.join(bug_group.split('/')[-1].split('-bug-analysis'))
    Path(parent_path+bug_name).mkdir(parents=True, exist_ok=True)
    for hist_name in hist_names:
        if bug_name[-7:]=='ens_avg':
            bug_nc0_csvs = PipeLine.file_reader(csv_files,extensions=[hist_name+'-ens_avg.csv'])
        else:
            bug_nc0_csvs = PipeLine.file_reader(csv_files,extensions=[hist_name+'.csv'])
        bug_nc0_csvs = [csv[0] for csv in bug_nc0_csvs]
        PipeLine.distributions_nc_zero(bug_nc0_csvs, parent_path+bug_name+'/', old_name=hist_name, new_name=hist_name+'Mon', empty=False)
        PipeLine.distributions_nc_zero(bug_nc0_csvs, parent_path+bug_name+'/', old_name=hist_name, new_name=hist_name+'Crd', empty=True)
        #bug_nc0_csvs_ens_avg = PipeLine.file_reader(database,extensions=[hist_name+'-ens_avg.csv'])
        #bug_nc0_csvs_ens_avg = [csv[0] for csv in bug_nc0_csvs_ens_avg]
        #PipeLine.distributions_nc_zero(bug_nc0_csvs_ens_avg, parent_path+bug_name+'nc0-ens_avg/', old_name=hist_name, new_name=hist_name+'Mon', empty=False)