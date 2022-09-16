#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/sum_rule/4-probe/sbatch*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/sum_rule/4-probe/submit*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/sum_rule/4-probe/probe-*histdd_transverse_size.py .
cp ../PolyPhys/bash_pipeline/sum_rule/4-probe/probe_gnuparallel*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/sum_rule/4-probe/create_probe_directories.sh .
cp -R ../PolyPhys/polyphys .