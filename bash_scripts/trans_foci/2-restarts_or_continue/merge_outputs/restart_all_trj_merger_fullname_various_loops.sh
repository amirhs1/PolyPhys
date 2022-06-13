#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc6405dt0.005bdump1000adump5000ens[1-8]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j18.all.lammpstrj"
    resfile="${simname}.j19.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j18.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 8414 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc6004dt0.005bdump1000adump5000ens[1-2]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j19.all.lammpstrj"
    resfile="${simname}.j20.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j19.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 8013 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc6004dt0.005bdump1000adump5000ens[3-8]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j18.all.lammpstrj"
    resfile="${simname}.j19.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j18.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 8013 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5604dt0.005bdump1000adump5000ens{5,6,8}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j20.all.lammpstrj"
    resfile="${simname}.j21.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j20.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7613 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5604dt0.005bdump1000adump5000ens{2,3,7}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j21.all.lammpstrj"
    resfile="${simname}.j22.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j21.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7613 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5604dt0.005bdump1000adump5000ens{1,4}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j22.all.lammpstrj"
    resfile="${simname}.j23.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j22.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7613 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5204dt0.005bdump1000adump5000ens8; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j20.all.lammpstrj"
    resfile="${simname}.j21.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j20.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7213 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5204dt0.005bdump1000adump5000ens[4-7]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j21.all.lammpstrj"
    resfile="${simname}.j22.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j21.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7213 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5204dt0.005bdump1000adump5000ens3; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j23.all.lammpstrj"
    resfile="${simname}.j24.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j23.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7213 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc5204dt0.005bdump1000adump5000ens[1-2]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j22.all.lammpstrj"
    resfile="${simname}.j23.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j22.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 7213 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4804dt0.005bdump1000adump5000ens{1,7,8}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j21.all.lammpstrj"
    resfile="${simname}.j22.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j21.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6813 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4804dt0.005bdump1000adump5000ens6; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j22.all.lammpstrj"
    resfile="${simname}.j23.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j22.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6813 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4804dt0.005bdump1000adump5000ens[2-5]; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j20.all.lammpstrj"
    resfile="${simname}.j21.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j21.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6813 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4403dt0.005bdump1000adump5000ens{8,6,3}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j23.all.lammpstrj"
    resfile="${simname}.j24.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j23.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6412 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4403dt0.005bdump1000adump5000ens{7,2,1}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j22.all.lammpstrj"
    resfile="${simname}.j23.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j22.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6412 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4403dt0.005bdump1000adump5000ens4; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j24.all.lammpstrj"
    resfile="${simname}.j25.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j24.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6412 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done

for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4003dt0.005bdump1000adump5000ens{2,4,5,6,7}; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j24.all.lammpstrj"
    resfile="${simname}.j25.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j24.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6012 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done


for dir in N2000epsilon5.0r15.5lz379.5sig4.0nc4003dt0.005bdump1000adump5000ens1; do
    cd "$dir" || exit
    simname=$( echo "$dir" | cut -d / -f 1)
    echo "merging $simname parts..."
    incomfile="${simname}.j23.all.lammpstrj"
    resfile="${simname}.j24.all.lammpstrj"
    # change the name of restart file:
    mv "$resfile" "$simname.last.all.cont.res.lammpstrj"
    # restart file new name
    resnew="${simname}.last.all.cont.res.lammpstrj"
    tstep="$(head -n 2 "$resnew" | tail -n 1)"
    echo "timestep to split: $tstep"
    csplit "$incomfile" /"$tstep"/
    rm xx01
    sed -i '$d' xx00 # delete the last line of the xx00 file
    mv "$incomfile" "${simname}.last.all.cont.incomplete.lammpstrj"
    # cat the restart file to the end of corrected file and copy both to a new file
    fullfile="${simname}.j23.all.lammpstrj"
    cat xx00 "$resnew" > "$fullfile"
    echo "complete"
    rm xx00
    tail -n 6012 "$fullfile" | head -n 9
    echo 
    wc -l "$fullfile"
    cd ..
done