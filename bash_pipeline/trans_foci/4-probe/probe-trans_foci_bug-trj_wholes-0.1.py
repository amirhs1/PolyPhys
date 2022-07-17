from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober

# analyzing bug files.
geometry = 'biaxial'
group = 'bug'
lineage = 'whole'
save_to = './'
bug_pairs = glob('./eps*' + group + '*')
bug_pairs = organizer.sort_filenames(
    bug_pairs,
    fmts=['.' + group + '.data', '.' + group + '.lammpstrj']
)
for (bug_topo, bug_trj) in bug_pairs:
    prober.trans_fuci_bug(
        bug_topo,
        bug_trj,
        lineage,
        geometry,
        save_to=save_to
    )
