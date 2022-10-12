#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/sbatch*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/submit*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/probe-*.py .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/probe_gnuparallel-*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/create_probe_directories-*.sh .
cp ../PolyPhys/bash_pipeline/trans_foci/4-probe/cubic/create_probe_directories-next_probes-*.sh .
cp -R ../PolyPhys/polyphys .