#!/bin/bash
# Formula to find the line number: 
# LINE_NUMBER = (((TIME_STEP_OF_INTEREST - INITIAL_TIME_STEP) + 1)/ DUMP_FREQ) * (N_ATOMS + 9)
# Seeing the heading of the time step of interest:
# head -n LINE_NUMBER trj_file | tail -n (N_ATOMS + 9) | head -n 9


# Check if an input file name is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input_file"
    exit 1
fi

# File names
input_file="$1"  # Read from command line argument
output_file="missing_timestep.txt"

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
#rm temp.txt

# Output the result
echo "Filtered lines are saved in $output_file"