#!/bin/bash

# Search for directories starting with 'N'
for dir in N*ens[1-8]/; do
    # Check if submit.sh exists in the directory
    if [[ -f "${dir}submit.sh" ]]; then
        # Use sed to replace 'rrg-byha' with 'def-group' in submit.sh
        sed -i 's/rrg-byha/def-group/g' "${dir}submit.sh"
    fi
done