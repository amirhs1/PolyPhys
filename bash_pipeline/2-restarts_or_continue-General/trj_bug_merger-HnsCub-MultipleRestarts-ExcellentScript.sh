#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
echo "Start merging"
for res in N*_res/;do
    dir=$(echo "$res" | cut -d _ -f 1 )
    echo "$dir"
    mkdir "$dir"
    mkdir "$dir"/test
    cp ./"${res}"*.lammpstrj ./"$dir"/
    cp ./"${res}"*.lammpstrj.gz ./"$dir"/
    cp ./"${res}"*.data ./"$dir"/
    cp ./"${res}"log.lammps ./"${dir}"/"${dir}".restart.log
    cp ./"${dir}"_incomplete/*.lammpstrj ./"$dir"/
    cp ./"${dir}"_incomplete/*.lammpstrj.gz ./"$dir"/
    cp ./"${dir}"_incomplete/log.lammps ./"${dir}"/"${dir}".incomplete.log
    day=$(grep -oE "[0-9]{2}-[0-9]{2}:[0-9]{2}" ./"${dir}"_incomplete/submit.sh | cut -d '-' -f 1)
    dayToHrs=$(( day * 24))
    hoursMins=$(grep -oE "[0-9]{2}-[0-9]{2}:[0-9]{2}" ./"${dir}"_incomplete/submit.sh | cut -d '-' -f 2)
    hours=$(echo "$hoursMins" | cut -d ':' -f 1)
    totalHours=$(( hours + dayToHrs ))
    mins=$(echo "$hoursMins" | cut -d ':' -f 2)
    cd "$dir" || exit
        cp N*nucleoid*lammpstrj ./test/
        cd test || exit
            nMol=$(head -n 4 "$dir".nucleoid.restart.lammpstrj | tail -n 1)
            dumpStep=2000
            lastStep=$(head -n 2 "${dir}".nucleoid.restart.lammpstrj | tail -n 1)
            splitStep=$((lastStep + dumpStep))
            echo "timestep to split: $splitStep"
            csplit "${dir}".nucleoid.lammpstrj /$splitStep/
            rm xx01 # this is the rest of incomplete trajectory which is not needed.
            ##sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
            sed -i '$d' xx00 # linux: delete the last line of the xx00 file
            mv "${dir}".nucleoid.lammpstrj "${dir}".nucleoid.incomplete.lammpstrj
            linePerStep=$((nMol + 9))
            # the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
            echo "line per step: $linePerStep"
            #sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx
            sed -i "1,${linePerStep}d" "${dir}".nucleoid.restart.lammpstrj # linux
            cat xx00 "${dir}".nucleoid.restart.lammpstrj > "${dir}".nucleoid.lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
            echo complete
            rm xx00
            tail -n $((nMol + 9)) "${dir}".nucleoid.lammpstrj | head -n 9
            echo
            wc -l "${dir}".nucleoid.lammpstrj
            echo
        cd ..
        sed -i '$d' "${dir}".incomplete.log
        echo "Total wall time: ${totalHours}:${mins}:00" >> "${dir}".incomplete.log
        tail -n 2 N*.log
        echo
    cd ..
done
echo "Finished!"