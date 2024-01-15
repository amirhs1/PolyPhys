#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl/submit*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl/probe-*.py .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl/probe_gnuparallel*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCyl/create_probe_directories-SumRuleCyl.sh .
cp -R ../PolyPhys/polyphys .