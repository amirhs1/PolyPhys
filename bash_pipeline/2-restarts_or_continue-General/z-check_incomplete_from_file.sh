#!/bin/bash
# Path to the file
file="incomplete.txt"

# Check if the file exists
if [[ -f "$file" ]]; then
    # Read each line from the file
    while IFS= read -r folder
    do
        echo "Processing folder: $folder"
        # Add your code here to process each folder
        tail -n 2 ${folder}_res/log.lammps
        tail -n 1 ${folder}_res/slurm*
        filename=${folder}_res/log.lammps
        echo $filename
        awk '/^Performance:/ {count++; if (count==2) {print; exit}}' $filename
        awk '/^Performance:/ {count++; if (count==3) {print; exit}}' $filename
        awk '/^Performance:/ {count++; if (count==4) {print; exit}}' $filename
        filename=${folder}_incomplete/log.lammps
        echo $filename
        awk '/^Performance:/ {count++; if (count==2) {print; exit}}' $filename
        awk '/^Performance:/ {count++; if (count==3) {print; exit}}' $filename
        awk '/^Performance:/ {count++; if (count==4) {print; exit}}' $filename
        echo
    done < "$file"
else
    echo "File not found: $file"
fi