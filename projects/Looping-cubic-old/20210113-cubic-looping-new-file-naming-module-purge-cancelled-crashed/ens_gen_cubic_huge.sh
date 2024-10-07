#!/bin/bash

#			 N		l	sig3	n_crowd randseed_group dt dump sig2 n_large
P[${#P[@]}]="50		15.0	0.3		190986	5020 0.002	2000 4.0 2"
P[${#P[@]}]="50		14.0	0.3		232918	5030 0.002	2000 4.0 2"
#P[${#P[@]}]="50		24.0	0.3		195570	5020"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig3=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		randseed_group=`echo ${P[$j]} | awk '{print $5}'`
		dt=`echo ${P[$j]} | awk '{print $6}'`
		dump=`echo ${P[$j]} | awk '{print $7}'`
		sig2=`echo ${P[$j]} | awk '{print $8}'`
		n_large=`echo ${P[$j]} | awk '{print $9}'`
		dirname=N${N}dl${sig2}nl${n_large}l${l}sig${sig3}nc${n_crowd}ens${i}dt${dt}dump${dump}
# The command I should use:
		cp cubic_looping.lmp input.lmp

		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_bug equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig3 equal ${sig3}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp chain.50_looping.data ${dirname}
		if [ ${n_crowd} -ne 0 ]; then
			if [ ${n_crowd} -le 35000 ]; then
				cp submit_huge_nc_le_35000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 35000  ] && [ ${n_crowd} -le 70000  ]; then
				cp submit_huge_nc_b_35000_le_70000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 70000  ] && [ ${n_crowd} -le 100000  ]; then
				cp submit_huge_nc_b_70000_le_100000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 100000  ] && [ ${n_crowd} -le 200000  ]; then
				cp submit_huge_nc_b_100000_le_200000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_huge_nc_b_200000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
