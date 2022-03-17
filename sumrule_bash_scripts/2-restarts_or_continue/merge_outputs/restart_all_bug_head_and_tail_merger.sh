#simname=$( echo "$dir" | cut -d / -f 1)
ncrowd=5204
nmon=2000
totalj=22
prej=$(( $totalj - 1))
simname="N2000epsilon5.0r15.5lz379.5sig4.0nc5204dt0.005bdump1000adump5000ens6"
echo "merging $simname bug trjs..."
lastcom=${simname}.j02.bug.lammpstrj
extra=${simname}.j03.bug.lammpstrj
lstep=301001000
tstep=$(head -n 2 $extra | tail -n 1)
fstep=$(( ${tstep} + 1000))
echo "timestep to remove extra steps: $lstep"
csplit $extra "/$lstep/"
rm xx01
sed -i '$d' xx00 # delete the last line of the xx00 file
mv xx00 ${simname}.extra.lines.bug.lammpstrj
csplit ${simname}.extra.lines.bug.lammpstrj "/$fstep/"
rm xx00
sed -i '1s/^/ITEM: TIMESTEP\n/' xx01
mv xx01 ${simname}.rest.lines.bug.lammpstrj
# cat the restart file to the end of corrected file and copy both to a new file
fullfile=${simname}.j02.bug.complete.lammpstrj
cat ${lastcom} ${simname}.rest.lines.bug.lammpstrj > "$fullfile"
echo "complete"
nlines=$(( $nmon + 9))
tail -n $nlines $fullfile | head -n 9
echo 
wc -l $fullfile
rm $lastcom
rm $extra
rm ${simname}.extra.lines.bug.lammpstrj
rm ${simname}.rest.lines.bug.lammpstrj
# all trjs
echo "merging $simname all trjs..."
lastcom=${simname}.j${prej}.all.lammpstrj
extra=${simname}.j${totalj}.all.lammpstrj
lstep=301005000
tstep=$(head -n 2 $extra | tail -n 1)
fstep=$(( ${tstep} + 5000))
echo "timestep to remove extra steps: $lstep"
csplit $extra "/$lstep/"
rm xx01
sed -i '$d' xx00 # delete the last line of the xx00 file
mv xx00 ${simname}.extra.lines.all.lammpstrj
csplit ${simname}.extra.lines.all.lammpstrj "/$fstep/"
rm xx00
sed -i '1s/^/ITEM: TIMESTEP\n/' xx01
mv xx01 ${simname}.rest.lines.all.lammpstrj
# cat the restart file to the end of corrected file and copy both to a new file
fullfile=${simname}.j${prej}.all.complete.lammpstrj
cat ${lastcom} ${simname}.rest.lines.all.lammpstrj > "$fullfile"
echo "complete"
nlines=$(( $ncrowd + $nmon + 9))
tail -n $nlines $fullfile | head -n 9
echo 
wc -l $fullfile
rm $lastcom
rm $extra
rm ${simname}.extra.lines.all.lammpstrj
rm ${simname}.rest.lines.all.lammpstrj