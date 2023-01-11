#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl-histdd_transverse_size/sbatch*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl-histdd_transverse_size/submit*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl-histdd_transverse_size/probe-*histdd_transverse_size.py .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl-histdd_transverse_size/probe_gnuparallel*histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl-histdd_transverse_size/create_probe_directories*histdd_transverse_size.sh .
cp -R ../PolyPhys/polyphys .