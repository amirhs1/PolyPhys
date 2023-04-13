# use this file to probe trj bug wholes in cylindrical geometry
from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober

# analyzing bug files.
group = 'nucleoid'
lineage = 'whole'
save_to = './'
bug_pairs = glob('./N*' + group + '*')
bug_pairs = organizer.sort_filenames(
    bug_pairs,
    fmts=['.' + group + '.data', '.' + group + '.lammpstrj']
)
for (bug_topo, bug_trj) in bug_pairs:
    prober.hns_nucleoid_cyl(bug_topo, bug_trj, lineage, save_to=save_to)
