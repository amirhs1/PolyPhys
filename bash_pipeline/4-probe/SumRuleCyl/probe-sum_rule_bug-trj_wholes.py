# use this scripts if there is one "bug" trjs and several "all" trjs
from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober

# analyzing bug files:
group = 'bug'
lineage = 'whole'
save_to = './'
bug_pairs = glob('./N*' + group + '*')
bug_pairs = organizer.sort_filenames(
    bug_pairs,
    fmts=['.' + group + '.data', '.' + group + '.lammpstrj']
)
for (bug_topo, bug_trj) in bug_pairs:
    prober.sum_rule_bug_cyl(bug_topo, bug_trj, lineage, save_to=save_to)