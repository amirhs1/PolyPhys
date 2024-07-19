#!/bin/bash
# This script produce a report of each simulation containing input variables, a list of produced files, and log summaries of each run.
name=$(pwd | rev | cut -d / -f 1 | rev) #report name
rname=${name}-summary.txt
> ${rname} # name of the report file is the date of siumaltion plus a keyword refering to the main referernce of the simulation.
echo "Simulation name: $rname" >> $rname.txt
echo -e "\n" >> $rname.txt
for directory in N*[1-8]/; do
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

# move data and trajectory files to the bug_extraction directory
bash trj_extractor_v*.sh
# move Lammps runnning files to run_files directory
rname=run_files-${name}
mkdir ${rname}
mv *.data *.lmp *.sh ${rname}
