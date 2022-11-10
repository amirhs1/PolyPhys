#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
echo "If there is a number the same as the timestep of interest, this script does not work."
read -rp "Enter ens number > " i
read -rp "Enter dump step > " dumpStep
read -rp "Enter total number of monomers > " nMon
incompleteFile=$(echo "al*ens$i.ring.bug.lammpstrj")
restartFile=$(echo "al*ens$i.ring.bug.restart.lammpstrj")
lastStep="$(head -n 2 $restartFile | tail -n 1)"
splitStep=$(("${lastStep}"+"$dumpStep"))
echo "timestep to split: $splitStep"
csplit $incompleteFile /$splitStep/
rm xx01 # this is the rest of incomplete trajectory which is not needed.
#sed -i '' '$d' xx00 # mac osx: delete the last line of the xx00 file
sed -i '$d' xx00 # linux: delete the last line of the xx00 file

splitName=$(echo $incompleteFile | awk '{split($0,a,"lammpstrj"); print a[1]}') #name of incomplete trj (incomfile):
mv $incompleteFile ${splitName}incomplete.lammpstrj
fullFile=${splitName}lammpstrj
echo $fullFile
linePerStep=$(("$nMon"+9))
echo $linePerStep
sed -i "1,${linePerStep}d" ${restartFile} # linux: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
#sed -i '' "1,${linePerStep}d" "${restartFile}" # mac osx: the last of timestep of the split trajectory and the first timestep of the restart trajectory are the same., so delete the first timestep of restart trajectory.
cat xx00 $restartFile > $fullFile # cat the restart file to the end of corrected file and copy both to a new file
echo complete
rm xx00
tail -n $(("$nMon" + 9)) $fullFile | head -n 9
echo 
wc -l "$fullFile"