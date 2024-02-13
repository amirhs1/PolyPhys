from glob import glob
from polyphys.manage.parser import HnsCyl
from polyphys.analyze import clusters
import pandas as pd
import numpy as np

lj_cut = np.round(2**(1/6), 6)
dmon = 1
dhcore = 1
dhpatch_cut = 0.2
dhpatch = np.round(dhpatch_cut/lj_cut, 3)
dist_m_hpatch = 0.5 * (dhpatch_cut + dmon)
sigma_mon_hpatch = np.round(dist_m_hpatch / lj_cut, 3)
rcut_mon_hpatch = sigma_mon_hpatch + 0.15*dmon
rcut_m_hpatch_min = sigma_mon_hpatch
rcut_m_hpatch_max = rcut_mon_hpatch
rcut_m_hpatch_range = np.round(
    np.arange(rcut_m_hpatch_min, rcut_m_hpatch_max+0.01, 0.01), 3)
rcut_m_hpatch_range.sort()

probe_path = "./N*-distMatTMonHnsPatch.npy"
paths = sorted(glob(probe_path))
bindings = {
    'dcrowd': [],
    'phi_c_bulk': [],
    'nhns': [],
    'rcut': [],
    'cis_threshold': [],
    'n_m_hpatch_bound': [],
    'n_hpatch_free': [],
    'n_hpatch_engaged': [],
    'n_hcore_free': [],
    'n_hcore_bridge': [],
    'n_hcore_dangle': [],
    'n_hcore_cis': [],
    'n_hcore_trans': [],
    'whole': [],
    'ensemble_id': []
}
m_m_gen_dist_dict = {
    'dcrowd': [],
    'phi_c_bulk': [],
    'nhns': [],
    'rcut': [],
    'cis_threshold': [],
    'loop_size_hist': [],
    'whole': [],
    'ensemble_id': []
}
binding_dfs = []
m_m_gen_dist_dfs = []
for p in paths:
    dist_m_hpatch_t = np.load(p)
    print(p)
    sim_info = HnsCyl(
        p,
        'whole',
        'cylindrical',
        'nucleoid',
        'ring'
    )
    for cis_threshold in [2, 4]:
        print(cis_threshold)
        for rcut in rcut_m_hpatch_range:
            print(rcut)
            results = {
                'n_m_hpatch_bound': [],
                'n_hpatch_free': [],
                'n_hpatch_engaged': [],
                'n_hcore_free': [],
                'n_hcore_bridge': [],
                'n_hcore_dangle': [],
                'n_hcore_cis': [],
                'n_hcore_trans': []
            }
            if sim_info.topology == 'linear':
                loop_size_hist = np.zeros(sim_info.nmon, dtype=int)
            elif sim_info.topology == 'ring':
                loop_size_hist = np.zeros((sim_info.nmon//2)+1, dtype=int)
            else:
                raise ValueError(
                    "The genomic distance is not defined for "
                    f"'{sim_info.topology}' topology."
                    )
            bindings['ensemble_id'].append(sim_info.ensemble_id)
            bindings['whole'].append(sim_info.whole)
            bindings['dcrowd'].append(sim_info.dcrowd)
            bindings['phi_c_bulk'].append(round(sim_info.phi_c_bulk, 3))
            bindings['nhns'].append(sim_info.nhns)
            bindings['rcut'].append(rcut)
            bindings['cis_threshold'].append(cis_threshold)
            m_m_gen_dist_dict['ensemble_id'].append(sim_info.ensemble_id)
            m_m_gen_dist_dict['whole'].append(sim_info.whole)
            m_m_gen_dist_dict['dcrowd'].append(sim_info.dcrowd)
            m_m_gen_dist_dict['phi_c_bulk'].append(
                round(sim_info.phi_c_bulk, 3))
            m_m_gen_dist_dict['nhns'].append(sim_info.nhns)
            m_m_gen_dist_dict['rcut'].append(rcut)
            m_m_gen_dist_dict['cis_threshold'].append(cis_threshold)
            n_frames = len(dist_m_hpatch_t)
            for mat in dist_m_hpatch_t:
                direct_contact_m_hpatch = clusters.find_direct_contacts(
                    mat, rcut, inclusive=False)
                results, loop_size_hist = clusters.hns_binding(
                    direct_contact_m_hpatch,
                    sim_info.topology,
                    cis_threshold,
                    binding_stats=results,
                    loop_length_hist=loop_size_hist)
            bindings['n_m_hpatch_bound'].append(np.sum(
                np.array(results['n_m_hpatch_bound'], dtype=int))/n_frames
                                             )
            bindings['n_hpatch_engaged'].append(np.sum(
                np.array(results['n_hpatch_engaged'], dtype=int))/n_frames
                                               )
            bindings['n_hpatch_free'].append(np.sum(
                np.array(results['n_hpatch_free'], dtype=int))/n_frames
                                            )
            bindings['n_hcore_free'].append(np.sum(
                np.array(results['n_hcore_free'], dtype=int))/n_frames
                                           )
            bindings['n_hcore_bridge'].append(np.sum(
                np.array(results['n_hcore_bridge'], dtype=int))/n_frames
                                             )
            bindings['n_hcore_dangle'].append(np.sum(
                np.array(results['n_hcore_dangle'], dtype=int))/n_frames
                                             )
            bindings['n_hcore_cis'].append(
                np.sum(np.array(results['n_hcore_cis'], dtype=int))/n_frames)
            bindings['n_hcore_trans'].append(
                np.sum(np.array(results['n_hcore_trans'], dtype=int))/n_frames)
            m_m_gen_dist_dict['loop_size_hist'].append(loop_size_hist/n_frames)
            loop_size_hist = {}
            max_loop_length = (sim_info.nmon//2)+1
            for key, value in m_m_gen_dist_dict.items():
                if key == 'loop_size_hist':
                    pass
                else:
                    loop_size_hist[key] = value
            nens = len(loop_size_hist['dcrowd'])
            for i in range(max_loop_length):
                loop_size_hist[i] = []
            for i in range(nens):
                hist_dummy = m_m_gen_dist_dict['loop_size_hist'][i]
                for j in range(max_loop_length):
                    loop_size_hist[j].append(hist_dummy[j])
            loop_size_hist_df = pd.DataFrame.from_dict(loop_size_hist)
            m_m_gen_dist_dfs.append(loop_size_hist_df)
            binding_rcut_df = pd.DataFrame.from_dict(bindings)
            binding_dfs.append(binding_rcut_df)
    output = f"{sim_info.whole}-HnsCyl-choice_of_threshold-loopLengthHistMon.csv"
    m_m_gen_dist_dfs = pd.concat(m_m_gen_dist_dfs)
    m_m_gen_dist_dfs.to_csv(output, index=False)
    output = f"{sim_info.whole}-HnsCyl-choice_of_threshold-bindingStat.csv"
    binding_dfs = pd.concat(binding_dfs)
    binding_dfs.to_csv(output, index=False)
