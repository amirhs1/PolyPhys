#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
read -rp "Enter ens number > " i

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
splitname=$(echo ${incomfile} | awk '{split($0,a,"lammpstrj"); print a[1]}')
incomplete="${splitname}incomplete.lammpstrj"
mv "$incomfile" "$incomplete"

# cat the restart file to the end of corrected file and copy both to a new file
fullfile=${splitname}lammpstrj
cat xx00 "$resfile" > "$fullfile"
echo complete
rm xx00
tail -n 2009 "$fullfile" | head -n 9
echo 
wc -l "$fullfile"