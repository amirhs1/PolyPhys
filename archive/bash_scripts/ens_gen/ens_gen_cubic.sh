#!/bin/bash

#			 N		l	sig2	n_crowd
P[${#P[@]}]="50		36.0	0.3		0"
P[${#P[@]}]="50		36.0	0.3		990100"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		dirname=N${N}ens${i}l${l}r${r}epsilon${espsilon1}sig${sig2}nc${n_crowd}
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
		if [ ${n_crowd} -eq 0 ]; then
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		else
			cp submit_huge.sh submit.sh && mv submit.sh ${dirname}
		fi
		cd ${dirname}
		[[ -d archive ]] || mkdir archive
		cd ..

		
#		echo "cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done