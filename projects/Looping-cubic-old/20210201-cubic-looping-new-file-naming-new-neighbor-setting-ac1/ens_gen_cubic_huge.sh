#!/bin/bash

#			 N		l	sig3	n_crowd	randseed_group	run_dt bug_dump all_dump sig2 n_large m_large
P[${#P[@]}]="100	30	1.0		82506	5090 0.002	2000	20000	5.0 2	125"
P[${#P[@]}]="100	30	1.0		103133	5090 0.002	2000	20000	5.0 2	125"
P[${#P[@]}]="100	30	1.0		123759	5090 0.002	2000	20000	5.0 2	125"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig3=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		randseed_group=`echo ${P[$j]} | awk '{print $5}'`
		run_dt=`echo ${P[$j]} | awk '{print $6}'`
		bug_dump=`echo ${P[$j]} | awk '{print $7}'`
		all_dump=`echo ${P[$j]} | awk '{print $8}'`
		sig2=`echo ${P[$j]} | awk '{print $9}'`
		n_large=`echo ${P[$j]} | awk '{print $10}'`
		m_large=`echo ${P[$j]} | awk '{print $11}'`
		dirname=N${N}dl${sig2}nl${n_large}l${l}dc${sig3}nc${n_crowd}dt${dt}bdump${bug_dump}adump${all_dump}ens${i}
# The command I should use:
		cp cubic_looping.lmp input.lmp

		
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_bug equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig3 equal ${sig3}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable mass2 equal ${m_large}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable bug_dump equal ${bug_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp chain.100_looping.data ${dirname}
		if [ ${n_crowd} -ne 0 ]; then
			if [ ${n_crowd} -le 35000 ]; then
				cp submit_huge_nc_le_35000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 35000  ] && [ ${n_crowd} -le 70000  ]; then
				cp submit_huge_nc_b_35000_le_70000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 70000  ] && [ ${n_crowd} -le 100000  ]; then
				cp submit_huge_nc_b_70000_le_100000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 100000  ] && [ ${n_crowd} -le 150000  ]; then
				cp submit_huge_nc_b_100000_le_150000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_huge_nc_b_150000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
