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

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		dirname=N${N}l${r}sig${sig2}nc${n_crowd}
#		echo "mkdir ${dirname} && cd ${dirname} && mkdir archive && lmp_icc_openmpi < ${dirname}.in  > log.out.${dirname}.txt" >> input.txt
		echo "([[ -d ${dirname}]] || mkdir ${dirname}) && cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${dirname}.txt" >> input.txt
	done
done
