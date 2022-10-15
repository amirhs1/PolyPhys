#!/bin/bash

#			 n_small	epsilon_small epsilon_big	r	lz	sig3	n_crowd	randseed_group	run_dt bug_dump all_dump sig2 n_big m_big
P[${#P[@]}]="350	5.0	5.0	15.5	56.5	1.0		0	40000	0.005	1000	5000	5.0 7	125"
P[${#P[@]}]="350	5.0	5.0	15.5	56.5	1.0		48867	40010	0.005	1000	5000	5.0 7	125"

ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		n_small=`echo ${P[$j]} | awk '{print $1}'` # total number of monomers
		epsilon_small=`echo ${P[$j]} | awk '{print $2}'`
		epsilon_big=`echo ${P[$j]} | awk '{print $3}'`
		r=`echo ${P[$j]} | awk '{print $4}'`
		lz=`echo ${P[$j]} | awk '{print $5}'` # range of x, y or z from -l to l
		sig3=`echo ${P[$j]} | awk '{print $6}'` # crowder size
		n_crowd=`echo ${P[$j]} | awk '{print $7}'`
		randseed_group=`echo ${P[$j]} | awk '{print $8}'`
		run_dt=`echo ${P[$j]} | awk '{print $9}'`
		bug_dump=`echo ${P[$j]} | awk '{print $10}'`
		all_dump=`echo ${P[$j]} | awk '{print $11}'`
		sig2=`echo ${P[$j]} | awk '{print $12}'` # big monomer size
		n_big=`echo ${P[$j]} | awk '{print $13}'` # number of big monomers
		m_big=`echo ${P[$j]} | awk '{print $14}'`
		dirname=ns${n_small}epss${epsilon_small}epsb${epsilon_big}r${r}dl${sig2}nl${n_big}ml${m_big}dc${sig3}nc${n_crowd}lz${lz}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}ring
# The command I should use:
		if [ $((${n_crowd}+${n_small}+${n_big})) -ge 25000 ]; then
			cp foci-ring-cylinder-ac_equal-cores_equal_more_8.lmp input.lmp
		else
			cp foci-ring-cylinder-ac_equal-cores_less_8.lmp input.lmp
		fi 

		echo "variable n_big equal ${n_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_small equal ${epsilon_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_big equal ${epsilon_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_small equal ${n_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig3 equal ${sig3}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable mass2 equal ${m_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable bug_dump equal ${bug_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp *.data initial_config.data
		mv initial_config.data ${dirname}


		N=$((${n_big}+${n_small}))
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
