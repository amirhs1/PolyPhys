# Importing necessary packages:
from glob import glob
from PipeLine import *

csv_files = glob("./N*-analyze_segments/N*.csv")
#csv_files = glob("../N2000D25.0ac1.0phic0.325_0.4-analyze_segments/N*.csv")

hist_names = ['rHistsCrd','rHistsMon','zHistsCrd','zHistsMon','thetaHistsCrd','thetaHistsMon','rPhisCrd','rPhisMon','rRhosCrd','rRhosMon','zPhisCrd','zPhisMon','zRhosCrd','zRhosMon']
geometry = 'cylindrical'
seperator = '-'
for hist_name in hist_names:
    ext = [seperator+hist_name+'.csv']
    hist_files = PipeLine.file_reader(csv_files,extensions=ext)
    hist_dicts = PipeLine.dict_of_ens_dfs(hist_files,hist_name,geometry,single=True,sep=seperator,index_col=0,dtype=float,skiprows=1)
    ens_evg_dict_rhist = PipeLine.ensemble_avg(hist_dicts,hist_name,geometry)
    