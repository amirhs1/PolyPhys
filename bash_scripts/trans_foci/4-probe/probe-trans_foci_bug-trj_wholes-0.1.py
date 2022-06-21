from glob import glob
from polyphys.manage import organizer
from polyphys.probe import prober

# analyzing bug files.
geometry = 'biaxial'
lineage = 'whole'
bug_pairs = glob("./eps*.bug.*")
bug_pairs = organizer.sort_filenames(bug_pairs)
for (bug_topo, bug_trj) in bug_pairs:
    prober.trans_fuci_bug(
        bug_topo,
        bug_trj,
        lineage=lineage,
        geometry=geometry
    )
