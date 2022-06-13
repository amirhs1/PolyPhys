# Importing necessary packages:

from glob import glob
from PipeLine import *


fnames = glob("../sumrule_test/*.bug.*")
fnames = PipeLine.file_reader(fnames)

geom = 'cylinder'

for fname in fnames:
    PipeLine.analyze_trj(fname, geom)
    PipeLine.trj_rmsd(fname, geom)
    