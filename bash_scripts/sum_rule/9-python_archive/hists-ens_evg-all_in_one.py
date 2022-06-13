from glob import glob
import pandas as pd

from PipeLine import *



properties_path = "../sumrule_data/properties-ens_avg-all_in_one.csv"
properties_ens_avg = pd.read_csv(properties_path,header=0)

xname = 'r'
yname = 'phi_r'
xnorm = 'dcyl'
ynorm = 'phi' # it can be 'rho', 'phi', 'hist'
ext = ['-rPhisCrd-ens_avg.csv','-rPhisMon-ens_avg.csv']
all_files = glob("../sumrule_data/*-analyze_all/N*.csv")

rPhis_ens_avg_csvs = PipeLine.file_reader(all_files,extensions=ext)
rPhis_ens_avg_df = PipeLine.organize_hists_ens_evg(rPhis_ens_avg_csvs,properties_ens_avg,xname,yname,xnorm,ynorm,filename='rPhis-ens_avg-all_in_one',round_to=6)