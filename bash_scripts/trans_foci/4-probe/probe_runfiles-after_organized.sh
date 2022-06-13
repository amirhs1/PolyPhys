#!/bin/bash
# copies necessary files for a job array on the slurm schedular for prbing after file organized
cp ../PolyPhys/sumrule_bash_scripts/4-probe/sbatch*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4-probe/submit*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4-probe/probe-*.py .
cp ../PolyPhys/sumrule_bash_scripts/4-probe/probe_gnuparallel*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4-probe/create_probe_directories.sh .
cp -R ../PolyPhys/polyphys .