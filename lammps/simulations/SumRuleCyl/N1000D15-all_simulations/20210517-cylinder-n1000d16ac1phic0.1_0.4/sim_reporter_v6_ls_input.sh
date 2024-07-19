#!/bin/bash
# This script produce a report of each simulation containing input variables, a list of produced files, and log summaries of each run.
rname=$(pwd | rev | cut -d / -f 1 | rev) #report name
> $rname.txt # name of the report file is the date of siumaltion plus a keyword refering to the main referernce of the simulation.
echo "Simulation name: $rname" >> $rname.txt
echo -e "\n" >> $rname.txt
for directory in N*/; do
    cd ${directory}
    echo $directory | cut -d / -f 1
    #echo $runname
    echo "Run name: "$directory | cut -d / -f 1 >> ../$rname.txt
    echo '1.Inputs:' >> ../$rname.txt
    head -n 7 input.lmp | tail -n 6 >> ../$rname.txt
    echo -e "\n" >> ../$rname.txt
    echo '2. Run summary:' >> ../$rname.txt
    tail -n 45 log.lammps >> ../$rname.txt
    echo -e "\n" >> ../$rname.txt
    echo '3. Created files:' >> ../$rname.txt
    ls -la --block-size=M  >> ../$rname.txt
    echo -e "\n" >> ../$rname.txt
    echo -e "\n" >> ../$rname.txt
    echo -e "#------------------------------------------------------#\n" >> ../$rname.txt
    echo -e "#------------------------------------------------------#" >> ../$rname.txt
    cd ..
done
