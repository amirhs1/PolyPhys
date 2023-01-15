#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/4-probe/TransFociCyl-hist2d_transverse_size/submit*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/TransFociCyl-hist2d_transverse_size/probe*.py .
cp ../PolyPhys/bash_pipeline/4-probe/TransFociCyl-hist2d_transverse_size/probe_gnuparallel-*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/TransFociCyl-hist2d_transverse_size/create_probe_directories-*.sh .
cp -R ../PolyPhys/polyphys .
