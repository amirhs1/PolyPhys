#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/4-probe/HnsCub/submit*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/HnsCub/probe-*.py .
cp ../PolyPhys/bash_pipeline/4-probe/HnsCub/probe_gnuparallel-*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/HnsCub/create_probe_directories-*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/HnsCub/create_probe_directories-next_probes-*.sh .
cp -R ../PolyPhys/polyphys .