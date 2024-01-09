#!/bin/bash

# Check if missing_runs.txt exists
if [[ -f "missing_runs.txt" ]]; then
    # Read each line in missing_runs.txt
    while IFS= read -r dir; do
        # Check if the directory exists
        if [[ -d "$dir" ]]; then
            # Define file names
            input_file="${dir}/${dir}.bug.lammpstrj"
            output_file="${dir}-missing_timestep.txt"

            # Find lines after "ITEM: TIMESTEP" and save to a temporary file
            awk '/ITEM: TIMESTEP/{getline; print}' "$input_file" > temp.txt

            # Initialize variables
            previous=0

            # Process the temporary file
            while read -r line; do
                # Calculate the difference
                difference=$((line - previous))

                # Check if the difference is not 5000
                if [ "$difference" -ne 5000 ]; then
                    # If so, output the line
                    echo "$line" >> "$output_file"
                fi

                # Update the previous line value
                previous=$line
            done < temp.txt

            # Clean up the temporary file
            rm temp.txt

            # Output the result
            echo "Filtered lines are saved in $output_file"
        else
            echo "Directory $dir does not exist."
        fi
    done < missing_runs.txt
else
    echo "missing_runs.txt does not exist."
fi

