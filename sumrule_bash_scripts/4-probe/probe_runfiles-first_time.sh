#!/bin/bash
# copies necessary files for a job array on the slurm schedular
cp ../PolyPhys/sumrule_bash_scripts/4a-probe-first_time/sbatch*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4a-probe-first_time/submit*.sh .
cp ../PolyPhys/sumrule_bash_scripts/4a-probe-first_time/probe-*.py .
cp -R ../PolyPhys/polyphys .
cp ../PolyPhys/sumrule_bash_scripts/4a-probe-first_time/probe_gnuparallel*.sh .