# use this scripts if there is one "bug" trjs and several "all" trjs
from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober
from polyphys.manage.parser import SumRuleCyl

geometry = 'biaxial'
group = 'all'
topo_lineage = 'whole'
lineage = 'segment'
save_to = './'
all_trjs = glob('./N*' + group + '*')
all_trjs = organizer.sort_filenames(
    all_trjs,
    fmts=['.' + group + '.lammpstrj']
)
all_trjs = [all_trj[0] for all_trj in all_trjs]
all_topo = glob('./N*' + group + '*')
all_topo = organizer.sort_filenames(all_topo, fmts=['.' + group + '.data'])
all_topo = all_topo[0][0]
max_segment_id = len(all_trjs)
# analyzig all files
# it is assumed that the all trjs are numbers from 1 to max_segment_id
for all_trj in all_trjs:
    trj_info = SumRuleCyl(
        all_trj,
        geometry,
        group,
        topo_lineage,
        'linear'
    )
    # all the frames in the last segment are probed:
    if trj_info.segment_id == max_segment_id:
        prober.sum_rule_all_cyl_hist2d(
            all_topo, all_trj, lineage, save_to=save_to
        )
    else:
        prober.sum_rule_all_cyl_hist2d(
            all_topo,
            all_trj,
            lineage,
            save_to=save_to,
            continuous=True
        )
