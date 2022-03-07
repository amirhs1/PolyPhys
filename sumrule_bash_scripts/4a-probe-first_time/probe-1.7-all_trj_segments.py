from glob import glob
from polyphys.manage import organizer
from polyphys.manage.parser import SumRule
from polyphys.probe import prober


geometry = 'biaxial'
trj_lineage = 'segment'
save_to="./"
all_trjs = glob("./N*.all.lammpstrj")
all_trjs = organizer.sort_filenames(all_trjs,extensions=['.all.lammpstrj'])
all_trjs = [all_trj[0] for all_trj in all_trjs]
all_topo = glob("./N*.all.data")
all_topo = organizer.sort_filenames(all_topo,extensions=['.all.data'])
all_topo = all_topo[0][0]
print(all_topo)
for all_trj in all_trjs:
    print(all_trj)
    trj_info = SumRule(all_trj, geometry=geometry, group='all',
                       lineage=trj_lineage)
    # all the frames in the last segment are probed:
    if trj_info.segment_id ==len(all_trjs):
        prober.probe_all(trj_topo, bug_trj, geometry, trj_lineage, save_to)
    # the last frame in the all other segments is ignored:
    else:
        prober.probe_all(trj_topo, bug_trj, geometry, trj_lineage, save_to,
                         continuous=True)
