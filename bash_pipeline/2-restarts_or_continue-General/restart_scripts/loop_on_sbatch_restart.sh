#!/bin/bash

for dir in N*/; do # SumRuleCyl
#for dir in al*/; do # TransFociCub
#for dir in eps*/; do # TransFociCyl
   cd "${dir}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done