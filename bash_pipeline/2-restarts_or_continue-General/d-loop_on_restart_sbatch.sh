#!/bin/bash
for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear/; do # TransFociCub, SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit.sh
   sleep 5
   cd ..
done
