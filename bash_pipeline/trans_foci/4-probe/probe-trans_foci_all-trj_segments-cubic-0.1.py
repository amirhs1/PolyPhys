from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober
from polyphys.manage.parser import TransFoci, TransFociCubic


# analyzig all files
geometry = 'cubic'
parsers = {
    'biaxial': TransFoci,
    'cubic': TransFociCubic
}
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
# it is assumed that the all trjs are numbers from 1 to max_segment_id
for all_trj in all_trjs:
    trj_info = TransFoci(
        all_trj,
        geometry=geometry,
        group=group,
        lineage=topo_lineage,
        parser=parsers[geometry]
    )
    # all the frames in the last segment are probed:
    if trj_info.segment_id == max_segment_id:
        prober.trans_foci_all(
            all_topo, all_trj, geometry, lineage, save_to=save_to
        )
    else:

        prober.trans_foci_all(
            all_topo,
            all_trj,
            geometry,
            lineage,
            save_to=save_to,
            parser=parsers[geometry],
            continuous=True
        )
