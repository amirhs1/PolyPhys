#!/bin/bash
for dir in N*ens[1-8]_res/; do # SumRuleCyl
#for dir in N*ring_res/; do # HnsCub HnsCyl
#for dir in al*ring_res/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear_res/; do # SumRuleCubHeteroLinear
#for dir in eps*_res/; do # TransFociCyl
   cd "${dir}" || exit
   echo "${dir}"
   sbatch submit.sh
   sleep 5
   cd ..
done
