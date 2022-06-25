#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
ens="4 5 6 7 8"
for i in $ens;do
    dir="N2000epsilon5.0r10.5lz336sig1.0nc90720dt0.005bdump1000adump5000ens${i}/"
    cd "$dir" || exit
    echo 
    incomfile="bug.i$i.lammpstrj"
    echo "$incomfile" "..."
    resfile="bug.i$i.restart.lammpstrj"
    echo "$resfile" "..."
    tstep="$(head -n 2 "$resfile" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file

    #name of incomplete trj (incomfile):
    splitname=$(echo "$incomfile" | awk '{split($0,a,"lammpstrj"); print a[1]}')
    incomplete="${splitname}incomplete.lammpstrj"
    mv "$incomfile" "$incomplete"

    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile=${splitname}lammpstrj
    cat xx00 "$resfile" > "$fullfile"
    echo merging complete
    rm xx00
    tail -n 2009 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    rm "$resfile"
    rm "$incomplete"
    echo "done"
    cd ..
done