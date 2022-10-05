#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/histdd_transverse_size/sbatch*_histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/histdd_transverse_size/submit*_histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/histdd_transverse_size/probe-*_histdd_transverse_size.py .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/histdd_transverse_size/probe_gnuparallel*_histdd_transverse_size.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/histdd_transverse_size/create_probe_directories_histdd_transverse_size.sh .
cp -R ../PolyPhys/polyphys .
