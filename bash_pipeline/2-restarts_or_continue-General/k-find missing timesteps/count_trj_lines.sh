#!/bin/bash

# Loop through directories matching the pattern "N*ens[1-4]"
for dir in N*ens[1-4]; do
    # Check if the directory exists and is a directory
    if [[ -d "$dir" ]]; then
        # Count the number of lines in the file "${dir}.bug.lammpstrj"
        line_count=$(wc -l < "${dir}/${dir}.bug.lammpstrj" 2>/dev/null)

        # Check if the line count is not "120542009"
        if [[ "$line_count" != "120542009" ]]; then
            # Add the directory name to "missing_runs.txt"
            echo "$dir" >> missing_runs.txt
        fi
    fi
done

