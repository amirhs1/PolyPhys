#!/bin/bash
# CHECK: the total line number in bug file below
# CHECK: the pattern of bug files 
# Loop through directories matching the pattern "N*ens[1-4]"
for dir in al*linear; do
    # Check if the directory exists and is a directory
    if [[ -d "$dir" ]]; then
        # Count the number of lines in the file "${dir}.bug.lammpstrj"
        line_count=$(wc -l < "${dir}/${dir}.bug.lammpstrj" 2>/dev/null)

        # Check if the line count is not "20700414"
        if [[ "$line_count" != "20700414" ]]; then
            # Add the directory name to "missing_runs.txt"
            echo "$dir" >> missing_runs.txt
        fi
    fi
done
echo "Finished!"
