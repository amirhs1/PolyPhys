
#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
dumpStep=5000
echo "Start merging"
for res in N*_cont/;do
    dir=$(echo "$res" | cut -d _ -f 1 )
    echo "$dir"
    mkdir "$dir"/test
    cd "$dir" || exit
        cp "${dir}".bug.lammpstrj ./test/
        cp "${dir}".bug.continue.lammpstrj ./test/
        cd test || exit
            nMol=$(head -n 4 "$dir".bug.continue.lammpstrj | tail -n 1)
            lastStep=$(head -n 2 "${dir}".bug.continue.lammpstrj | tail -n 1)
            splitStep=$lastStep
            echo "timestep to split: $splitStep"
            csplit "${dir}".bug.lammpstrj /$splitStep/
            rm xx01 # this is the rest of incomplete trajectory which is not needed.
            ##sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
            sed -i '$d' xx00 # linux: delete the last line of the xx00 file
            mv "${dir}".bug.lammpstrj "${dir}".bug.missing_continue.lammpstrj
            linePerStep=$((nMol + 9))
            # the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
            echo "line per step: $linePerStep"
            #sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx
            #sed -i "1,${linePerStep}d" "${dir}".bug.continue.lammpstrj # linux
            cat xx00 "${dir}".bug.continue.lammpstrj > "${dir}".bug.lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
            echo complete
            rm xx00
            tail -n $((nMol + 9)) "${dir}".bug.lammpstrj | head -n 9
            echo
            wc -l "${dir}".bug.lammpstrj
            echo
        cd ..
    cd ..
done
echo "Finished!"
