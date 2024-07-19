#!/bin/bash

#			 N		l	sig2	n_crowd
P[${#P[@]}]="50		36.0	0.3		990100"

ens='13 14 15 16' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		dirname=N${N}nc${n_crowd}ens${i}l${l}sig${sig2}
# The command I should use:
		cp cubic.lmp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp	
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp			
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp	
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname} 
		cp data.chain.50 ${dirname}
		cp submit_huge.sh submit.sh && mv submit.sh ${dirname}
		cd ${dirname}
		[[ -d archive ]] || mkdir archive
		cd ..
	done
done
