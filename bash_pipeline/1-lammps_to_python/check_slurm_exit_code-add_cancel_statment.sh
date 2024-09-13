#!/bin/bash

# Loop through all directories that match the pattern "al*linear*"
for slurm_file in al*linear*/slurm*.out; do
    # Check if the file contains a line starting with "Program finished with exit"
    if grep -q "Program finished with exit code" "$slurm_file"; then
        # Extract the exit code from the matching line
        exit_code=$(grep "Program finished with exit code" "$slurm_file" | awk '{print $6}')
        # Check if the exit code is not 0
        if [ "$exit_code" -ne 0 ]; then
            # Print the folder name if exit code is not 0
            folder=$(dirname $slurm_file)
            echo "Exit code $exit_code found in $folder"
        fi
    fi
done

