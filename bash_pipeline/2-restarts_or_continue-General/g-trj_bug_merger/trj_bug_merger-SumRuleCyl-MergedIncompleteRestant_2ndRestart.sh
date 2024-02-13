#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
dumpStep=5000
echo "Start merging"
for res in N*[1-8]/;do
    dir=$(echo "$res" | cut -d / -f 1 )
    echo "$dir"
    mkdir "$dir"/test
    cd "$res" || exit
        cp ${dir}.j02.bug.lammpstrj ./test/
        cp ${dir}.j02.bug.2nd_restart.lammpstrj ./test/
        cd test || exit
            nMol=$(head -n 4 "$dir".j02.bug.2nd_restart.lammpstrj | tail -n 1)
            lastStep=$(head -n 2 "${dir}".j02.bug.2nd_restart.lammpstrj | tail -n 1)
            splitStep=$((lastStep + dumpStep))
            echo "timestep to split: $splitStep"
            csplit "${dir}".j02.bug.lammpstrj /$splitStep/
            rm xx01 # this is the rest of incomplete trajectory which is not needed.
            ##sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
            sed -i '$d' xx00 # linux: delete the last line of the xx00 file
            mv "${dir}".j02.bug.lammpstrj "${dir}".j02.bug.incomplete_and_restart.lammpstrj
            linePerStep=$((nMol + 9))
            # the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
            echo "line per step: $linePerStep"
            #sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx
            sed -i "1,${linePerStep}d" "${dir}".j02.bug.2nd_restart.lammpstrj # linux
            cat xx00 "${dir}".j02.bug.2nd_restart.lammpstrj > "${dir}".j02.bug.lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
            echo complete
            rm xx00
            tail -n $((nMol + 9)) "${dir}".j02.bug.lammpstrj | head -n 9
            echo
            wc -l "${dir}".j02.bug.lammpstrj
            echo
        cd ..
        echo
    cd ..
done
echo "Finished!"