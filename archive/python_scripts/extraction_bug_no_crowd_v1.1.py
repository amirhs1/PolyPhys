# Importing necessary packages:
from glob import glob
from PipeLine import *

fname = glob("./N*.bug.*")
fname = PipeLine.file_reader(fname) # This is a list with one member

geom = 'cylinder'
print(fname)
PipeLine.extract_trj_bug(fname[0], geom) # A list with one member, the member is a tuple of a trj and data pair.
PipeLine.bug_trj_rmsd(fname[0], geom)
