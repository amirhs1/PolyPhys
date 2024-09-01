#!/bin/bash

#for dir in al*[1-8].ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*[1-8].linear/; do # SumRuleCubHeteroLinear
#for dir in N*.ring/; do #HnsCub
#for dir in N*[1-8]/; do #SumRuleCul
for dir in eps*[1-8].ring/; do #TransFociCyl
    # Check if submit.sh exists in the directory
    if [[ -f "${dir}submit_system_minimize.sh" ]]; then
        # Use sed to replace 'rrg-byha' with 'def-byha' in submit.sh
        sed -i 's/rrg-byha/def-byha/g' "${dir}submit_system_minimize.sh"
    fi
done