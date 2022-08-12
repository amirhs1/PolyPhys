# use this scripts if there is several "bug" trjs
from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober
from polyphys.manage.parser import SumRule


geometry = 'biaxial'
group = 'bug'
topo_lineage = 'whole'
lineage = 'segment'
save_to = './'
bug_trjs = glob('./N*' + group + '*')
bug_trjs = organizer.sort_filenames(
    bug_trjs,
    fmts=['.' + group + '.lammpstrj']
)
bug_trjs = [bug_trj[0] for bug_trj in bug_trjs]
bug_topo = glob('./N*' + group + '*')
bug_topo = organizer.sort_filenames(bug_topo, fmts=['.' + group + '.data'])
bug_topo = bug_topo[0][0]
max_segment_id = len(bug_trjs)
# analyzig all files
# it is assumed that the all trjs are numbers from 1 to max_segment_id
for bug_trj in bug_trjs:
    trj_info = SumRule(
        bug_trj,
        geometry=geometry,
        group=group,
        lineage=topo_lineage
    )
    # all the frames in the last segment are probed:
    if trj_info.segment_id == max_segment_id:
        prober.sum_rule_bug(
            bug_topo, bug_trj, geometry, lineage, save_to=save_to
        )
    else:
        prober.sum_rule_bug(
            bug_topo,
            bug_trj,
            geometry,
            lineage,
            save_to=save_to,
            continuous=True
        )
