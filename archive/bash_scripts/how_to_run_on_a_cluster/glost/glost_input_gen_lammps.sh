#!/bin/bash

#			 N		l	sig2	n_crowd	epsilon1	
P[${#P[@]}]="50		36	0.3		0		5.0"
P[${#P[@]}]="50		36	0.3		2500	5.0"
P[${#P[@]}]="50		36	0.3		5000	5.0"
P[${#P[@]}]="50		36	0.3		7500	5.0"
P[${#P[@]}]="50		36	0.3		10000	5.0"
P[${#P[@]}]="50		36	0.3		12500	5.0"
P[${#P[@]}]="50		36	0.3		15000	5.0"
P[${#P[@]}]="50		36	0.3		17500	5.0"
P[${#P[@]}]="50		36	0.3		20000	5.0"
P[${#P[@]}]="50		36	0.3		22500	5.0"
P[${#P[@]}]="50		36	0.3		25000	5.0"
P[${#P[@]}]="50		36	0.3		27500	5.0"
P[${#P[@]}]="50		36	0.3		30000	5.0"
P[${#P[@]}]="50		36	0.3		32500	5.0"

P[${#P[@]}]="50		36	0.4		0		5.0"
P[${#P[@]}]="50		36	0.4		2500	5.0"
P[${#P[@]}]="50		36	0.4		5000	5.0"
P[${#P[@]}]="50		36	0.4		7500	5.0"
P[${#P[@]}]="50		36	0.4		10000	5.0"
P[${#P[@]}]="50		36	0.4		12500	5.0"
P[${#P[@]}]="50		36	0.4		15000	5.0"
P[${#P[@]}]="50		36	0.4		17500	5.0"
P[${#P[@]}]="50		36	0.4		20000	5.0"
P[${#P[@]}]="50		36	0.4		22500	5.0"
P[${#P[@]}]="50		36	0.4		25000	5.0"
P[${#P[@]}]="50		36	0.4		30000	5.0"
P[${#P[@]}]="50		36	0.4		32500	5.0"

P[${#P[@]}]="50		36	0.5		0		5.0"
P[${#P[@]}]="50		36	0.5		2500	5.0"
P[${#P[@]}]="50		36	0.5		5000	5.0"
P[${#P[@]}]="50		36	0.5		7500	5.0"
P[${#P[@]}]="50		36	0.5		10000	5.0"
P[${#P[@]}]="50		36	0.5		12500	5.0"
P[${#P[@]}]="50		36	0.5		15000	5.0"
P[${#P[@]}]="50		36	0.5		17500	5.0"
P[${#P[@]}]="50		36	0.5		20000	5.0"
P[${#P[@]}]="50		36	0.5		22500	5.0"
P[${#P[@]}]="50		36	0.5		25000	5.0"

ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		dirname=N${N}ens${i}l${l}sig${sig2}nc${n_crowd}
# The command I should use:
		cp cubic-with_variables.lmp input_${dirname}.lmp
		echo "variable in_filename equal data.chain.50" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable i equal $i" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable l equal $l" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable sig2 equal ${sig2}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo "variable epsilon1 equal ${epsilon1}" | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp	
		echo '#Defining inputs' | cat - input_${dirname}.lmp > temp && mv temp input_${dirname}.lmp
		echo "mkdir ${dirname} && cd ${dirname} && mkdir archive && lmp_icc_openmpi < ${dirname}.in  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
# The command I would like to use
		echo "mkdir ${dirname} && cd ${dirname} && mkdir archive && lmp_icc_openmpi -var in_filename \"../data.chain.${N}\" -var l ${l} -var sig2 ${sig2} -var n_crowd ${n_crowd} -var epsilon1 ${epsilon1} -var i ${i} < ${dirname}.in  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done
