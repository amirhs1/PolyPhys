#!/bin/bash
# copies necessary files for a job array on the slurm schedular for probing for the first time
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCubHeteroLinear/submit*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCubHeteroLinear/probe-*.py .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCubHeteroLinear/probe_gnuparallel-*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCubHeteroLinear/create_probe_directories-*.sh .
cp ../PolyPhys/bash_pipeline/4-probe/SumRuleCubHeteroLinear/create_probe_directories-next_probes-*.sh .
cp -R ../PolyPhys/polyphys .