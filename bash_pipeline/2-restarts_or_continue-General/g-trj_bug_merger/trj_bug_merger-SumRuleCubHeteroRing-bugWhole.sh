
#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
dumpStep=2000
echo "Start merging"
#for res in al*ring_res/;do
for res in al*ring_res/;do
    dir=$(echo "$res" | cut -d _ -f 1 )
    echo "$dir"
    mkdir "$dir"/test
    # Extract start and end times
    file=$(find "${dir}_incomplete" -maxdepth 1 -name "slurm*" -type f)
    start_time=$(grep -oP 'Starting run at: \K.*' "$file")
    potential_end_time=$(grep -oP 'Program finished with exit code [0-9]+ at: \K.*|(?<=CANCELLED AT )\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$file" | tail -1)

    if [ "$potential_end_time" != "$start_time" ]; then
        end_time="$potential_end_time"
    else
        "The end time was not found."
    fi
    # Convert times to seconds since epoch
    start_sec=$(date -d "$start_time" '+%s')
    end_sec=$(date -d "$end_time" '+%s')

    # Calculate duration
    duration=$((end_sec - start_sec))
    # Calculate hours, minutes, and seconds
    hours=$((duration / 3600))
    minutes=$(( (duration % 3600) / 60 ))
    seconds=$((duration % 60))

    # Format the duration
    formatted_duration=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
    cd "$dir" || exit
        cp al*bug*lammpstrj ./test/
        cd test || exit
            nMol=$(head -n 4 "$dir".bug.restart.lammpstrj | tail -n 1)
            lastStep=$(head -n 2 "${dir}".bug.restart.lammpstrj | tail -n 1)
            splitStep=$((lastStep + dumpStep))
            echo "timestep to split: $splitStep"
            csplit "${dir}".bug.lammpstrj /$splitStep/
            rm xx01 # this is the rest of incomplete trajectory which is not needed.
            ##sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
            sed -i '$d' xx00 # linux: delete the last line of the xx00 file
            mv "${dir}".bug.lammpstrj "${dir}".bug.incomplete.lammpstrj
            linePerStep=$((nMol + 9))
            # the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
            echo "line per step: $linePerStep"
            #sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx
            sed -i "1,${linePerStep}d" "${dir}".bug.restart.lammpstrj # linux
            cat xx00 "${dir}".bug.restart.lammpstrj > "${dir}".bug.lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
            echo complete
            rm xx00
            tail -n $((nMol + 9)) "${dir}".bug.lammpstrj | head -n 9
            echo
            wc -l "${dir}".bug.lammpstrj
            echo
        cd ..
        sed -i '$d' "${dir}".incomplete.log
        echo "Total wall time: $formatted_duration" >> "${dir}".incomplete.log
        tail -n 2 al*.log
        echo
    cd ..
done
echo "Finished!"
