from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober
from polyphys.manage.parser import TransFociCyl

# analyzing bug files.
group = 'bug'
lineage = 'whole'
save_to = './'
bug_pairs = glob('./eps*' + group + '*')
bug_pairs = organizer.sort_filenames(
    bug_pairs,
    fmts=['.' + group + '.data', '.' + group + '.lammpstrj']
)
for (bug_topo, bug_trj) in bug_pairs:
    prober.trans_fuci_bug_cyl_transverse_size(
        bug_topo, bug_trj, lineage, save_to=save_to
    )

# analyzig all files
group = 'all'
topo_lineage = 'whole'
lineage = 'segment'
save_to = './'
all_trjs = glob('./eps*' + group + '*')
all_trjs = organizer.sort_filenames(
    all_trjs,
    fmts=['.' + group + '.lammpstrj']
)
all_trjs = [all_trj[0] for all_trj in all_trjs]
all_topo = glob('./eps*' + group + '*')
all_topo = organizer.sort_filenames(all_topo, fmts=['.' + group + '.data'])
all_topo = all_topo[0][0]
max_segment_id = len(all_trjs)
# analyzig all files
# it is assumed that the all trjs are numbers from 1 to max_segment_id
for all_trj in all_trjs:
    trj_info = TransFociCyl(
        all_trj, topo_lineage, 'cylindrical', group, 'ring'
        )
    # all the frames in the last segment are probed:
    if trj_info.segment_id == max_segment_id:
        prober.trans_foci_all_cyl_hist2d(
            all_topo, all_trj, lineage, save_to=save_to
        )
    else:

        prober.trans_foci_all_cyl_hist2d(
            all_topo,
            all_trj,
            lineage,
            save_to=save_to,
            continuous=True
        )
