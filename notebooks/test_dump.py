from glob import glob
from polyphys.manage.organizer import invalid_keyword, sort_filenames
from polyphys.manage import dump

# analyzing bug files.
group = 'nucleoid'
lineage = 'whole'
save_to = './'
#path = "/Users/amirhsi/OneDrive - University of Waterloo/PhD Research/Jupyter/N200epshm29kbmm10nh0ac2nc0l50dt0.005ndump2000adump5000ens1.ring"
#path = "/Users/amirhsi_mini/OneDrive - University of Waterloo/PhD Research/Jupyter/N200epshm29kbmm20nh48ac2nc0l60dt0.005ndump2000adump5000ens1.ring"
#path = "/Users/amirhsi_mini/research_data/hns_cubic-trjs/N*/N*"
path = "/Users/amirhsi_mini/research_data/N200epshm29kbmm20nh48ac2nc0l60dt0" \
       ".005ndump2000adump5000ens1.ring/N*"
nuc_pairs = glob(path + '/N*' + group + '*')
nuc_pairs = sort_filenames(
    nuc_pairs,
    fmts=['.' + group + '.data', '.' + group + '.lammpstrj']
)
print(nuc_pairs)
trj = dump.Dump(nuc_pairs[0][1])
