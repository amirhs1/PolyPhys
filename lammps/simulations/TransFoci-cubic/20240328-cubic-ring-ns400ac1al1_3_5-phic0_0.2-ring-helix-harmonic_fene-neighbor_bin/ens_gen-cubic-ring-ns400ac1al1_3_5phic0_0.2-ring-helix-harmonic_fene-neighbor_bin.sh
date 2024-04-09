#!/bin/bash
#			n_small	 l	a_crowd	n_crowd		randseed	run_dt 		b_dump 	a_dump 	a_large n_large m_large
P[${#P[@]}]="400	48	1		0			32000		0.005		2000	5000	1 		5		1"
P[${#P[@]}]="400	45	1		139229		32010		0.005		2000	5000	1 		5		1"
P[${#P[@]}]="400	42	1		169798		32020		0.005		2000	5000	1 		5		1"
P[${#P[@]}]="400	42	1		226397		32030		0.005		2000	5000	1 		5		1"
P[${#P[@]}]="400	62	1		0			31000		0.005		2000	5000	5 		5		125"
P[${#P[@]}]="400	45	1		139229		31010		0.005		2000	5000	5 		5		125"
P[${#P[@]}]="400	42	1		169798		31020		0.005		2000	5000	5 		5		125"
P[${#P[@]}]="400	42	1		226397		31030		0.005		2000	5000	5 		5		125"
P[${#P[@]}]="400	55	1		0			30000		0.005		2000	5000	3 		5		27"
P[${#P[@]}]="400	45	1		139229		30010		0.005		2000	5000	3 		5		27"
P[${#P[@]}]="400	42	1		169798		30020		0.005		2000	5000	3 		5		27"
P[${#P[@]}]="400	42	1		226397		30030		0.005		2000	5000	3 		5		27"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
#ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		n_small=$(echo "${P[$j]}" | awk '{print $1}') # number of small monomers
		l=$(echo "${P[$j]}" | awk '{print $2}') # range of x, y or z from -l to l
		sig3=$(echo "${P[$j]}" | awk '{print $3}') # crowder size
		n_crowd=$(echo "${P[$j]}" | awk '{print $4}')
		randseed_group=$(echo "${P[$j]}" | awk '{print $5}')
		run_dt=$(echo "${P[$j]}" | awk '{print $6}')
		bug_dump=$(echo "${P[$j]}" | awk '{print $7}')
		all_dump=$(echo "${P[$j]}" | awk '{print $8}')
		sig2=$(echo "${P[$j]}" | awk '{print $9}') # big monomer size
		n_big=$(echo "${P[$j]}" | awk '{print $10}') # number of big monomers
		m_big=$(echo "${P[$j]}" | awk '{print $11}')
		dirname=al${sig2}nl${n_big}ml${m_big}ns${n_small}ac${sig3}nc${n_crowd}l${l}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}.ring

		# modifying foci-ring-cylinder-init_config_minimize-harmonic.lmp
		cp foci-ring-cubic-init_config_minimize-harmonic.lmp input.lmp

		echo "variable n_big equal ${n_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_small equal ${n_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig3 equal ${sig3}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable mass2 equal ${m_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable bug_dump equal ${bug_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		# rename the input.lmp to min_init_config.lmp
		mv input.lmp minimize_initial_config.lmp
		cp foci-ring-cubic-ac_equal.lmp input.lmp

		echo "variable n_big equal ${n_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_small equal ${n_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig3 equal ${sig3}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable mass2 equal ${m_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable bug_dump equal ${bug_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv minimize_initial_config.lmp "${dirname}"
		mv input.lmp "${dirname}"
		cp ./*.data initial_config.data
		mv initial_config.data "${dirname}"

		N=$((n_big+n_small))
		if [ $((n_crowd+N)) -ne ${N} ]; 
		then
			if [ $((n_crowd+N)) -le 5000 ]; 
			then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh "$dirname"
			elif [ $((n_crowd+N)) -gt 5000 ] && [ $((n_crowd+N)) -le 15000 ]; then
				cp submit_cpu4_gt5000_le15000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 15000 ] && [ $((n_crowd+N)) -le 30000 ]; then
				cp submit_cpu8_gt15000_le30000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 30000 ] && [ $((n_crowd+N)) -le 75000 ]; then
				cp submit_cpu16_gt30000_le75000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 75000 ] && [ $((n_crowd+N)) -le 150000 ]; then
				cp submit_cpu32_gt75000_le150000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 150000 ] && [ $((n_crowd+N)) -le 200000 ]; then
				cp submit_cpu32_gt150000_le200000.sh submit.sh && mv submit.sh "${dirname}"
			else
				cp submit_cpu32_gt200000_le250000.sh submit.sh && mv submit.sh "${dirname}"
			fi
		else
			cp submit_nocrowd_8hrs.sh submit.sh && mv submit.sh "${dirname}"
		fi
		cd "${dirname}" || exit
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done