from glob import glob
from PipeLine import *

bug_files = glob("../sumrule_data/N*-analyze_bug/N*.csv")
property_files = PipeLine.file_reader(bug_files,extensions=['-properties.csv'])
properties_all_in_one = PipeLine.all_dfs_in_one(property_files, 'properties', index_col=0)
property_files = PipeLine.file_reader(bug_files,extensions=['-properties-ens_avg.csv'])
properties_all_in_one_ens_avg = PipeLine.all_dfs_in_one(property_files,'properties', ens_avg=True, norm_func=PipeLine.cyl_sumrule_norms_ens_evg, index_col=0)