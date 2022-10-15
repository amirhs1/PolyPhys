#!/bin/bash
# Script to run a lammps MD simjulation
#List of variables and their names defined by command line method. Their description and the step they will be used in are described in next sections:
#variable ${in_filename}
#variable ${i} # i is the seed for the RNG.
#variable ${r} # the radius of the cylinder.
#variable ${lz} # the length of the cylinder in the longitudinal direction z
#variable ${sig2} # diameter of a crowder
#variable ${n_crowd} # Number of the crowders
#variable ${epsilon1} # t

mkdir images # Make a directory for saving images
mkdir restarts # Make a directory for saving restart files

lmp_serial -v in_filename data.chain.80 -v i 1 -v r 3 -v lz 52 -v sig2 0.3 -v n_crowd 1000 -v epsilon1 5.0 -i  *.lmp
