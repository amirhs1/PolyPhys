#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
for dir in N*[1-8]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j01.bug.lammpstrj"
    resfile="${simname}.j02.bug.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.bug.cont.res.lammpstrj"
    # restart file new name:
    resnew="${simname}.bug.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file:
    mv "$incomfile" "${simname}.bug.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file:
    fullfile="${simname}.j02.bug.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 2009 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done