# Importing necessary packages:
from glob import glob
from PipeLine import *

ens_files = glob("./N*/N*.csv")
geometry = 'cylinder'

rhist_name = 'rHists'
rhist_files = PipeLine.file_reader(ens_files,extensions=['_rHists.csv','_rEdges.csv'])
rhist_dicts = PipeLine.ensemble(rhist_files,rhist_name,geometry,dtype=int,single=False)
ens_evg_dict_rhist = PipeLine.group(rhist_dicts,rhist_name,geometry)

zhist_name = 'zHists'
zhist_files = PipeLine.file_reader(ens_files,extensions=['_zHists.csv','_zEdges.csv'])
zhist_dicts = PipeLine.ensemble(zhist_files,zhist_name,geometry,dtype=int,single=False)
ens_evg_dict_zhist = PipeLine.group(zhist_dicts,zhist_name,geometry)

thetahist_name = 'thetaHists'
thetahist_files = PipeLine.file_reader(ens_files,extensions=['_thetaHists.csv','_thetaEdges.csv'])
thetahist_dicts = PipeLine.ensemble(thetahist_files,thetahist_name,geometry,dtype=int,single=False)
ens_evg_dict_thetahist = PipeLine.group(thetahist_dicts,thetahist_name,geometry)

rFloryhist_name = 'rFloryHists'
rFloryhist_files = PipeLine.file_reader(ens_files,extensions=['_rFloryHists.csv','_rFloryEdges.csv'])
rFloryhist_dicts = PipeLine.ensemble(rFloryhist_files,rFloryhist_name,geometry,dtype=int,single=False)
ens_evg_dict_rFloryhist = PipeLine.group(rFloryhist_dicts,rFloryhist_name,geometry)

fsd_name = 'fsd_t'
fsd_files = PipeLine.file_reader(ens_files,extensions=['_fsd_t.csv'])
fsd_dicts = PipeLine.ensemble(fsd_files,fsd_name,geometry)
ens_evg_dict_fsd = PipeLine.group(fsd_dicts,fsd_name,geometry)

gyr_name = 'gyr_t'
gyr_files = PipeLine.file_reader(ens_files,extensions=['_gyr_t.csv'])
gyr_dicts = PipeLine.ensemble(gyr_files,gyr_name,geometry)
ens_evg_dict_gyr = PipeLine.group(gyr_dicts,gyr_name,geometry)

rFlory_name = 'rFlory_t'
rFlory_files = PipeLine.file_reader(ens_files,extensions=['_rFlory_t.csv'])
rFlory_dicts = PipeLine.ensemble(rFlory_files,rFlory_name,geometry)
ens_evg_dict_rFlory = PipeLine.group(rFlory_dicts,rFlory_name,geometry)

properties_name = 'properties'
properties_files = PipeLine.file_reader(ens_files,extensions=['_properties.csv'])
properties_df = PipeLine.properties(properties_files, properties_name, geometry)
ens_avg_properties = PipeLine.ens_avg_properties(properties_df, properties_name, geometry)


all_rho_r_dicts, all_phi_r_dicts = PipeLine.distributions_generator(rhist_dicts, properties_df, 'dmon', 'cylindrical', 'radial', '')
ens_avg_dict_rho_r = PipeLine.group(all_rho_r_dicts,'rRhos',geometry)
ens_avg_dict_phi_r = PipeLine.group(all_phi_r_dicts,'rPhis',geometry)

all_rho_z_dicts, all_phi_z_dicts = PipeLine.distributions_generator(zhist_dicts, properties_df, 'dmon', 'cylindrical', 'longitudinal', '')
ens_avg_dict_rho_z = PipeLine.group(all_rho_z_dicts,'zRhos',geometry)
ens_avg_dict_phi_z = PipeLine.group(all_phi_z_dicts,'zPhis',geometry)



log_files = glob('./N*/N*.log')
geometry = 'cylinder'
log_files = PipeLine.file_reader(log_files,extensions=['.log'])
details_out , runtime_out = PipeLine.log_outputs(log_files[0][0], geometry=geometry)
PipeLine.lammps_log_details(log_files, details_out , runtime_out)