# Importing necessary packages:
from glob import glob
import numpy as np
import pandas as pd

from polyphys.analyze import clusters
from polyphys.manage.parser import TransFociCub
probe_path = "/Users/amirhsi_mini/research_data/ns*-probe/al*.ring/al*ens*-bug-distMatTFoci.npy"
# probe_path = "./al*ens1*-bug-distMatTFoci.npy"
paths = sorted(glob(probe_path))
lj_cut = 2**(1/6)
for path in paths:
    bonding = {
        'dcrowd': [],
        'phi_c_bulk': [],
        'dmon_large': [],
        'rcut': [],
        'bond_freq_mean': [],
        'ensemble_id': [],
        'whole': []
    }
    bonding_df_list = []
    print(path)
    foci_dist_t = np.load(path)
    sim_info = TransFociCub(
        path,
        'whole',
        'cubic',
        'bug',
        'ring'
    )
    whole = sim_info.whole
    rcut_min = sim_info.dmon_large + sim_info.dcrowd
    rcut_max = np.round(lj_cut * (sim_info.dmon_large + sim_info.dcrowd), 3)
    rcuts = np.arange(rcut_min, rcut_max, 0.035)
    # (d_m_large+d_crowd) <= rcut <= 1.122460*(dm_mon_large+d_crowd), so
    # 0.89 <= r_cut_norm <= 1
    # dr_cut = 0.035 so the number of rcut increases with d_mon_large.
    for rcut in rcuts:
        print(rcut)
        bonds_t = np.empty([0, sim_info.nmon_large], dtype=int)
        bonding['dcrowd'].append(sim_info.dcrowd)
        bonding['phi_c_bulk'].append(round(sim_info.phi_c_bulk, 3))
        bonding['dmon_large'].append(sim_info.dmon_large)
        bonding['rcut'].append(rcut)
        bonding['ensemble_id'].append(sim_info.ensemble_id)
        bonding['whole'].append(sim_info.whole)
        n_frames = len(foci_dist_t)
        # The foci_dist_t has a especial structure, at a given time, foci_dist
        # is upper trinagle matrix where the diagonal elements are the index of
        # foci, and the upper traingle elements are the distance between large
        # monomers i and j. This happened becauses pair distance matrix is
        # symmetric.
        for foci_dist in foci_dist_t:
            np.fill_diagonal(foci_dist, 0)
            upper_tri = np.triu(foci_dist)
            dist_mat = upper_tri + upper_tri.T - np.diag(np.diag(upper_tri))
            dir_contacts = clusters.find_direct_contacts(dist_mat, rcut)
            bonds_stat = clusters.count_foci_bonds(dir_contacts)
            bonds_t = np.append(bonds_t, np.array([bonds_stat]), axis=0)
        bonding['bond_freq_mean'].append(bonds_t.mean(axis=0))
        bonding_dict = {}
        bond_max_size = sim_info.nmon_large
        for key, value in bonding.items():
            if key == 'bond_freq_mean':
                pass
            else:
                bonding_dict[key] = value
        n_wholes = len(bonding_dict['dcrowd'])
        for i in range(bond_max_size):
            bonding_dict[i] = []
        for i in range(n_wholes):
            hist_dummy = bonding['bond_freq_mean'][i]
            for j in range(bond_max_size):
                bonding_dict[j].append(hist_dummy[j])
        bonding_df = pd.DataFrame.from_dict(bonding_dict)
        bonding_df_list.append(bonding_df)
    output = whole + "-choice_of_threshold-bonding.csv"
    bonding_dfs = pd.concat(bonding_df_list)
    bonding_dfs.to_csv(output, index=False)
