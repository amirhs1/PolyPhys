from glob import glob
import pandas as pd
from PipeLine import *

csv_files = glob("./N*-extraction_all/*all*.csv")
properties_csvs = glob("./N*-analyze_bug/N*-properties.csv")

all_properties = []
for properties_csv in properties_csvs:
    df = pd.read_csv(properties_csv,index_col=0)
    all_properties.append(df) 
all_properties = pd.concat(all_properties)
all_properties.reset_index(inplace=True,drop=True)

rhist_ensembles_crd = PipeLine.simulation_from_segments(csv_files,'rHistsCrd', geometry,['_rHistsCrd.csv','_rEdgesCrd.csv'])
rhist_ensembles_mon = PipeLine.simulation_from_segments(csv_files,'rHistsMon', geometry,['_rHistsMon.csv','_rEdgesMon.csv'])
_, _ = PipeLine.distributions_generator(rhist_ensembles_crd, all_properties, 'dcrowd',geometry,'radial','Crd',normalized=False)
_, _ = PipeLine.distributions_generator(rhist_ensembles_mon, all_properties, 'dmon', geometry,'radial','Mon',normalized=False)

zhist_ensembles_crd = PipeLine.simulation_from_segments(csv_files,'zHistsCrd',geometry,['_zHistsCrd.csv','_zEdgesCrd.csv'])
zhist_ensembles_mon = PipeLine.simulation_from_segments(csv_files,'zHistsMon',geometry,['_zHistsMon.csv','_zEdgesMon.csv'])
_, _ = PipeLine.distributions_generator(zhist_ensembles_crd, all_properties, 'dcrowd', geometry,'longitudinal','Crd',normalized=False)
_, _ = PipeLine.distributions_generator(zhist_ensembles_mon, all_properties, 'dmon', geometry,'longitudinal','Mon',normalized=False)

thetahist_ensembles_crd = PipeLine.simulation_from_segments(csv_files,'thetaHistsCrd', geometry,['_thetaHistsCrd.csv','_thetaEdgesCrd.csv'])
thetahist_ensembles_mon = PipeLine.simulation_from_segments(csv_files,'thetaHistsMon', geometry,['_thetaHistsMon.csv','_thetaEdgesMon.csv'])