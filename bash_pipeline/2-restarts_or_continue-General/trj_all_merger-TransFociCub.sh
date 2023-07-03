#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
echo "If there is a number the same as the timestep of interest, this script does not work."
read -rp "Enter whole name > " whole
read -rp "Trj 1 (Top part) > " trj1
read -rp "Segment number of all Trj 1 > " j1
read -rp "Trj 2 (Bottom part) > " trj2
read -rp "Segment number of Trj 2 > " j2
read -rp "Number of particles > " nMon
dumpStep=5000
lastStep=$(head -n 2 "${trj2}" | tail -n 1)
splitStep=$((lastStep + dumpStep))
echo "timestep to split: $splitStep"
csplit "${trj1}" /$splitStep/
rm xx01 # this is the rest of incomplete trajectory which is not needed.
#sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
sed -i '$d' xx00 # linux: delete the last line of the xx00 file
mv "${trj1}" "${whole}".j${j1}.ring.all.incomplete.lammpstrj
linePerStep=$((nMon + 9))
echo $linePerStep
sed -i "1,${linePerStep}d" "${trj2}" # linux: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
#sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
cat xx00 "${trj2}" > "${whole}".j${j1}.ring.all.lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
echo complete
rm xx00
tail -n $((nMon + 9)) "${whole}".j${j2}.ring.all.restart.lammpstrj | head -n 9
mergedLines=$(wc -l <"${whole}".j${j1}.ring.all.lammpstrj)
echo "Number of lines: $mergedLines"
dividedLines=$((mergedLines / (nMon + 9)))
echo "Number of snapshots/configurations: ${dividedLines}"