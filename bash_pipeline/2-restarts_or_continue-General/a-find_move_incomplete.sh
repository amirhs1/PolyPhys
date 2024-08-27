#!/bin/bash

# Check if the incomplete directory exists, if not, create it
if [ ! -d "incomplete" ]; then
    mkdir "incomplete"
fi

# Loop through directories starting with "N"
for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear/; do # TransFociCub, SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
    # Find the file starting with "slurm" in the current directory
    slurm_file=$(find "$dir" -maxdepth 1 -name "slurm*" -type f)

    # Check if the slurm file exists
    if [ -f "$slurm_file" ]; then
        # Read the last line of the file
        last_line=$(tail -n 1 "$slurm_file")

        # Check if the last line starts with the specific text
        if [[ "$last_line" != "Starting run at:"* ]] && [[ "$last_line" != "Program finished with exit code 0 at"* ]]; then
            # Move the directory to incomplete
            mv "$dir" "incomplete/"
        fi
    fi
done