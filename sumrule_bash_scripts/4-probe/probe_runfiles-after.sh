#!/bin/bash
# copies necessary files for a job array on the slurm schedular
cp ../PolyPhys/sumrule_bash_scripts/4b-probe-after_organized/sbatch*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4b-probe-after_organized/submit*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4b-probe-after_organized/probe-*.py .
cp -R ../PolyPhys/polyphys .
cp ../PolyPhys/sumrule_bash_scripts/4b-probe-after_organized/probe_gnuparallel*.sh .