# use this file to probe trj bug wholes in cubic geometry
from glob import glob
from polyphys.manage import organizer
import numpy as np
import pandas as pd
from polyphys.manage.parser import HnsCub
from polyphys.analyze import clusters
import traceback

save_to = './'
wholes = organizer.sort_filenames(
    glob("./N*-nucleoid-distMatTMonHnsPatch.npy"), fmts=['npy'])
for idx in range(len(wholes)):
    whole = wholes[idx][0]
    print(whole)
    dist_m_hp = np.load(whole)
    sim_info = HnsCub(
            whole,
            'whole',
            'cubic',
            'nucleoid',
            'ring'
        )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    lj_cut = 2**(1/6)
    r_cutoff = np.round(
        0.5 * lj_cut * (sim_info.dmon + sim_info.dhns_patch), 3
        )
    n_frames, n_mon, n_hpatch = dist_m_hp.shape
    n_hcore = n_hpatch // 2
    # Process each time step
    foci_t = np.zeros((n_frames, n_hcore, n_hcore))
    dir_contacts_t = np.zeros(
        (n_frames, n_hcore, n_hcore), dtype=int
    )
    bonds_t = np.zeros((n_frames, n_hcore), dtype=int)
    clusters_t = np.zeros((n_frames, n_hcore+1), dtype=int)

    bonds_dangler_t = np.zeros((n_frames, n_hcore), dtype=int)
    clusters_dangler_t = np.zeros((n_frames, n_hcore+1), dtype=int)

    bonds_bridger_t = np.zeros((n_frames, n_hcore), dtype=int)
    clusters_bridger_t = np.zeros((n_frames, n_hcore+1), dtype=int)
    error_logs = pd.DataFrame(columns=['Frame_Index', 'Error_Message'])
    for idx in range(n_frames):
        dir_m_hp_dirty = clusters.find_direct_contacts(
            dist_m_hp[idx], r_cutoff, inclusive=False
            )
        dir_m_hp = clusters.enforce_single_patch_dir_contact(dir_m_hp_dirty)
        dir_m_hc = clusters.generate_mon_bind_direct(dir_m_hp, 2)
        dir_m_hc_dangler, dir_m_hc_bridger = \
            clusters.split_binder_matrix(dir_m_hc)
        try:
            # bound binder direct contact matirx
            dir_hc_hc = clusters.find_binder_clusters(dir_m_hc)
            bonds_stat = clusters.count_foci_bonds(dir_hc_hc)
            bonds_t[idx] = bonds_stat
            full_hc_hc = clusters.generate_contact_matrix(dir_hc_hc)
            clusters_stat = clusters.count_foci_clusters(full_hc_hc)
            clusters_t[idx] = clusters_stat
            # dangling binder direct contact matirx
            dir_hc_hc = clusters.find_binder_clusters(dir_m_hc_dangler)
            bonds_stat = clusters.count_foci_bonds(dir_hc_hc)
            bonds_dangler_t[idx] = bonds_stat
            full_hc_hc = clusters.generate_contact_matrix(dir_hc_hc)
            clusters_stat = clusters.count_foci_clusters(full_hc_hc)
            clusters_dangler_t[idx] = clusters_stat
            # bridging binder direct contact matirx
            dir_hc_hc = clusters.find_binder_clusters(dir_m_hc_bridger)
            bonds_stat = clusters.count_foci_bonds(dir_hc_hc)
            bonds_bridger_t[idx] = bonds_stat
            full_hc_hc = clusters.generate_contact_matrix(dir_hc_hc)
            clusters_stat = clusters.count_foci_clusters(full_hc_hc)
            clusters_bridger_t[idx] = clusters_stat
        except Exception as e:
            # Log the error with frame index
            error_message = \
                f"Error at frame {idx} in whole {idx}: {str(e)}\n" + \
                f"Traceback: {traceback.format_exc()}"
            error_logs = error_logs.append(
                {'Frame_Index': idx,
                 'Error_Message': error_message},
                ignore_index=True)
            continue

    np.save(save_to + sim_name + '-bondsHistTHnsCore.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTHnsCore.npy', clusters_t)
    np.save(save_to + sim_name + '-bondsHistDangleTHnsCore.npy',
            bonds_dangler_t)
    np.save(save_to + sim_name + '-clustersHistDangleTHnsCore.npy',
            clusters_dangler_t)
    np.save(save_to + sim_name + '-bondsHistBridgeTHnsCore.npy',
            bonds_bridger_t)
    np.save(save_to + sim_name + '-clustersHistBridgeTHnsCore.npy',
            clusters_bridger_t)
    error_logs.to_csv(save_to + sim_name + 'error_logs.csv', index=False)
