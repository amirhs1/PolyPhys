#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/sbatch*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/submit*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/probe-*.py .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/probe_gnuparallel*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/create_probe_directories.sh .
cp -R ../PolyPhys/polyphys .