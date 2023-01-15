# use this scripts if there are several "all" trjs.
from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober
from polyphys.manage.parser import SumRuleCyl


group = 'all'
top_lineage = 'whole'
lineage = 'segment'
save_to = './'
all_trjs = glob('./N*' + group + '*')
all_trjs = organizer.sort_filenames(
    all_trjs,
    fmts=['.' + group + '.lammpstrj']
)
all_trjs = [all_trj[0] for all_trj in all_trjs]
all_top = glob('./N*' + group + '*')
all_top = organizer.sort_filenames(all_top, fmts=['.' + group + '.data'])
all_top = all_top[0][0]
max_segment_id = len(all_trjs)
# analyzig all files
# it is assumed that the all trjs are numbers from 1 to max_segment_id
for all_trj in all_trjs:
    trj_info = SumRuleCyl(all_trj, top_lineage, 'cylindrical', group, 'linear')
    # all the frames in the last segment are probed:
    if trj_info.segment_id == max_segment_id:
        prober.sum_rule_all_cyl(all_top, all_trj, lineage, save_to=save_to)
    else:
        prober.sum_rule_all_cyl(
            all_top,
            all_trj,
            lineage,
            save_to=save_to,
            continuous=True
        )
