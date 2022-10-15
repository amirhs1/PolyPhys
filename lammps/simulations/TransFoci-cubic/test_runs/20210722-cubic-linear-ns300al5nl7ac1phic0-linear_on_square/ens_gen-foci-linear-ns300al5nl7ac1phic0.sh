#!/bin/bash

#			 n_small	l	sig3	n_crowd	randseed_group	run_dt bug_dump all_dump sig2 n_large m_large
P[${#P[@]}]="300	96.5	1.0		0	20000 0.002	1000	5000	5.0 7	125"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		n_small=`echo ${P[$j]} | awk '{print $1}'` # total number of monomers
		l=`echo ${P[$j]} | awk '{print $2}'` # range of x, y or z from -l to l
		sig3=`echo ${P[$j]} | awk '{print $3}'` # crowder size
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		randseed_group=`echo ${P[$j]} | awk '{print $5}'`
		run_dt=`echo ${P[$j]} | awk '{print $6}'`
		bug_dump=`echo ${P[$j]} | awk '{print $7}'`
		all_dump=`echo ${P[$j]} | awk '{print $8}'`
		sig2=`echo ${P[$j]} | awk '{print $9}'` # large monomer size
		n_large=`echo ${P[$j]} | awk '{print $10}'` # number of large monomers
		m_large=`echo ${P[$j]} | awk '{print $11}'`
		dirname=ns${n_small}dl${sig2}nl${n_large}dc${sig3}nc${n_crowd}l${l}dt${dt}bdump${bug_dump}adump${all_dump}ens${i}linear
# The command I should use:
		cp foci-linear-ac_equal-cores_less_8.lmp input.lmp

		echo "variable n_large equal ${n_large}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_small equal ${n_small}" | cat - input.lmp > temp && mv temp input.lmp
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
		cp foci-linear-ns300al5nl7.data ${dirname}

		N=$((${n_large}+${n_small}))
		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu2_le5000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 15000 ]; then
				cp submit_cpu4_gt5000_le15000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 15000 ] && [ $((${n_crowd}+${N})) -le 25000 ]; then
				cp submit_cpu4_gt15000_le25000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 25000 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
				cp submit_cpu8_gt25000_le50000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
				cp submit_cpu16_gt50000_le100000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_cpu32_gt100000_le200000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd_8hrs.sh submit.sh && mv submit.sh ${dirname}
		fi
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
