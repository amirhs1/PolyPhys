#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
echo "If there is a number the same as the timestep of interest, this script does not work."
read -rp "Enter whole name > " whole
read -rp "Enter dump step > " dumpStep
read -rp "Enter total number of monomers > " nMon
lastStep=$(head -n 2 "$whole".restart.lammpstrj | tail -n 1)
splitStep=$((lastStep + dumpStep))
echo "timestep to split: $splitStep"
csplit "${whole}".lammpstrj /$splitStep/
rm xx01 # this is the rest of incomplete trajectory which is not needed.
#sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
sed -i '$d' xx00 # linux: delete the last line of the xx00 file
mv "${whole}".lammpstrj "${whole}".incomplete.lammpstrj
linePerStep=$((nMon + 9))
echo $linePerStep
sed -i "1,${linePerStep}d" "${whole}.restart.lammpstrj" # linux: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
#sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
cat xx00 "${whole}".restart.lammpstrj > "${whole}".lammpstrj # cat the restart file to the end of corrected file and copy both to a new file
echo complete
rm xx00
tail -n $((nMon + 9)) "${whole}".lammpstrj | head -n 9
echo 
wc -l "${whole}".lammpstrj