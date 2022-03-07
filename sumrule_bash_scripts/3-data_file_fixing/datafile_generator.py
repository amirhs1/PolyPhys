#!/bin/bash
# Importing necessary packages:
from glob import glob
from polyphys.construct.builder import data_file_generator

path = "./*.bug.lammpstrj"
fnames = glob(path)
dfile_template = glob("./data_template-cylinder_bug_n*.data")
data_file_generator(dfile_template, fnames)
