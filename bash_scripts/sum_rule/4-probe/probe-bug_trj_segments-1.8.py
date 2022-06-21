from glob import glob
from polyphys.manage import organizer
from polyphys.manage.parser import SumRule
from polyphys.probe import prober


geometry = 'biaxial'
trj_lineage = 'segment'
save_to = "./"
bug_trjs = glob("./N*.bug.lammpstrj")
bug_trjs = organizer.sort_filenames(bug_trjs, fmts=['.bug.lammpstrj'])
bug_trjs = [bug_trj[0] for bug_trj in bug_trjs]
bug_topo = glob("./N*.bug.data")
bug_topo = organizer.sort_filenames(bug_topo, fmts=['.bug.data'])
bug_topo = bug_topo[0][0]
print(bug_topo)
for bug_trj in bug_trjs:
    print(bug_trj)
    trj_info = SumRule(bug_trj, geometry=geometry, group='bug',
                       lineage=trj_lineage)
    # all the frames in the last segment are probed:
    if trj_info.segment_id == len(bug_trjs):
        prober.sum_rule_bug(
            bug_topo, bug_trj, geometry, trj_lineage, save_to
        )
    # the last frame in the all other segments is ignored:
    else:
        prober.sum_rule_bug(
            bug_topo, bug_trj, geometry, trj_lineage, save_to, continuous=True
        )
