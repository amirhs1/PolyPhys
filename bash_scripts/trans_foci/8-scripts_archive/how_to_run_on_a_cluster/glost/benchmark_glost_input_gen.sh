#!/bin/bash

#			 N		l	sig2	n_crowd	epsilon1 n_runs equ_run	
P[${#P[@]}]="50		36	0.3		0	5.0	100000 		10000"
P[${#P[@]}]="50		36	0.3		0	5.0 1000000		100000"
P[${#P[@]}]="50		36	0.3		0	5.0 10000000	1000000"

P[${#P[@]}]="50		36	0.3		450	5.0	100000		10000"
P[${#P[@]}]="50		36	0.3		450	5.0 1000000		100000"
P[${#P[@]}]="50		36	0.3		450	5.0 10000000	1000000"

P[${#P[@]}]="50		36	0.3		950	5.0	100000		10000"
P[${#P[@]}]="50		36	0.3		950	5.0 1000000		100000"
P[${#P[@]}]="50		36	0.3		950	5.0 10000000	1000000"

ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		n_runs=`echo ${P[$j]} | awk '{print $6}'`
		equ_runs=`echo ${P[$j]} | awk '{print $7}'`
		dirname=N${N}ens${i}l${l}sig${sig2}nc${n_crowd}nr${n_runs}
# The command I should use:
		cp cubic.lmp input_${dirname}.lmp
		echo "variable i equal $i" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable l equal $l" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable sig2 equal ${sig2}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable epsilon1 equal ${epsilon1}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp
		echo "variable n_runs equal ${n_runs}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable equ_runs equal ${equ_runs}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo '#Defining inputs' | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp
		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input_${dirname}.lmp ${dirname} 
		cp data.chain.50 ${dirname}
		echo "cd ${dirname} && mkdir archive && lmp_icc_openmpi < input_${dirname}.lmp  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input.txt
#		echo "cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done
