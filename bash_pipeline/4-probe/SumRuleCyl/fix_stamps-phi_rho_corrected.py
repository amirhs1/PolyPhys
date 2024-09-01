import os
from glob import glob
import numpy as np
import pandas as pd


stamps = glob('./stamps*/N*stamps.csv')
for file in stamps:
    print(file)
    stamp = pd.read_csv(file)
    vol_cell_m = np.pi * stamp['lcyl'] * (stamp['dcyl'] - stamp['dmon'])**2 / 4
    stamp['rho_m_bulk'] = stamp['nmon'] / vol_cell_m
    stamp['phi_m_bulk'] = np.pi * stamp['rho_m_bulk'] / 6 # dmon = 1.0
    vol_cell_c = np.pi * stamp['lcyl'] * (stamp['dcyl'] - stamp['dcrowd'])**2 / 4
    stamp['rho_c_bulk'] = stamp['ncrowd'] / vol_cell_c
    stamp['phi_c_bulk'] = np.pi * stamp['rho_c_bulk'] * stamp['dcrowd']**3 / 6
    output = file.split('/')[-1]
    space_dir = file.split('/')[1].split('-')[1]+'-probe'
    if not os.path.exists(f'./{space_dir}'):
        os.mkdir(f'./{space_dir}')
    stamp.to_csv(f'./{space_dir}/{output}',index=False)