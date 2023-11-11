#!/bin/bash
# Initialize a counter
counter=1
# Loop through files and modify them
cd "data_files"
for file in N2000epsilon*r*; do
    # New file name
    new_file="file_$counter"

    # Get the first line of the file
    first_line=$(head -n 1 "$file")

    # Add the original file name at the end of the first line
    # and then add the rest of the file
    (echo "$first_line $file"; tail -n +2 "$file") > "${new_file}.data"
    
    # Increment the counter
    ((counter++))
done
cd ..