#!/bin/bash
# produces a report for each simulation that contains:
# 1. Input variables, 
# 2. a list of generated files; and
# 3. log summaries of each run.
pname=$(pwd | rev | cut -d / -f 1 | rev) #parent name
name=$( echo "$pname" | cut -d - -f 1)
simtype=$( echo "$pname" | cut -d - -f 2 | cut -d _ -f 1) 
# simtype is whether "all" or "cont"
# "all" means that the first time simulations are run.
# "cont" means this is the 2nd time they are run since the system has not yet reached equilibrium.
if [ "$simtype" = "all" ]; then
    rname=${name}-summary # report name
elif [ "$simtype" = "cont" ]; then
    rname=${name}-${simtype}-summary # report name
else
    echo "not a 'all_simulations' or 'cont_simulations' dir"
fi
echo "Simulation name: $rname" > "$rname".txt # name of the report file
echo -e "\n" >> "$rname".txt
for dir in epss*[1-8].ring/; do
    cd "$dir" || exit
    simname=$(echo "$dir" | cut -d / -f 1)
    echo "$simname"
    #echo $runname
    echo "Run name: $simname" >> "../$rname.txt"
    echo "1.Inputs:" >> "../$rname.txt"
    head -n 7 input.lmp | tail -n 6 >> ../"$rname".txt
    echo -e "\n" >> ../"$rname".txt
    echo '2. Run summary:' >> ../"$rname".txt
    tail -n 45 log.lammps >> ../"$rname".txt
    echo -e "\n" >> ../"$rname".txt
    echo '3. Created files:' >> ../"$rname".txt
    ls -la --block-size=M  >> ../"$rname".txt
    echo -e "\n" >> ../"$rname".txt
    echo -e "\n" >> ../"$rname".txt
    echo -e "#------------------------------------------------------#\n" >> ../"$rname".txt
    echo -e "#------------------------------------------------------#" >> ../"$rname".txt
    cd ..
done