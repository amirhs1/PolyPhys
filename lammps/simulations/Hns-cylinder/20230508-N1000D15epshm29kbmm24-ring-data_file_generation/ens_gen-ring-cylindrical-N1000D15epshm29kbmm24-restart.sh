#!/bin/bash
# kbend is in LAMMPS lingo which is half the regualr one, so kbend_mm=2 in LAMMPS means kbend_mm=4 in theory.
#            N		kbend   n_hns   a_crowd n_crowd	r	lz		rseed_group   
P[${#P[@]}]="1000	24		0		2		0		8	628.5	1000000"
P[${#P[@]}]="1000	24		30		2		0		8	628.5	7000000"
P[${#P[@]}]="1000	24		30		2		4243	8	628.5	7000010"
P[${#P[@]}]="1000	24		30		2		21212	8	628.5	7000090"
P[${#P[@]}]="1000	24		30		1		33939	8	628.5	7000100"
P[${#P[@]}]="1000	24		30		1		169695	8	628.5	7000180"

#ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=$(echo "${P[$j]}" | awk '{print $1}') # number of monomers
		kBend=$(echo "${P[$j]}" | awk '{print $2}')
		n_hns=$(echo "${P[$j]}" | awk '{print $3}')
		sig4=$(echo "${P[$j]}" | awk '{print $4}') # sizeo of crowders
		n_crowd=$(echo "${P[$j]}" | awk '{print $5}')
		r=$(echo "${P[$j]}" | awk '{print $6}')
		lz=$(echo "${P[$j]}" | awk '{print $7}') # range of x, y or z from -l to l
		randseed_group=$(echo "${P[$j]}" | awk '{print $8}')
		dirname=N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}.ring
		simname=N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}
		# restart script
		cp restart-input-HnsCyl-MinimizeDna.lmp input.lmp
		echo "variable simname string ${simname}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable kBend11 equal ${kBend}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp restart_minimize_dna.lmp
		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv restart_minimize_dna.lmp "${dirname}/"
		cp submit_restart_minimize_dna.sh "${dirname}/" 
	done
done