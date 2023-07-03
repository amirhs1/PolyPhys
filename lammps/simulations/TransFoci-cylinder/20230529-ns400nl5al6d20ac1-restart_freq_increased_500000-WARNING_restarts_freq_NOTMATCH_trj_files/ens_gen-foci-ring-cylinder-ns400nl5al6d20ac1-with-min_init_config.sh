#!/bin/bash
#			n_small	eps_small 	eps_large	r		lz		a_crowd	n_crowd		randseed	run_dt 	b_dump 	a_dump 	a_large n_large m_large
P[${#P[@]}]="400	5			5			10.5	82.5		1		0			16000		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	82.5		1		9900		16010		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	82.5		1		14850		16020		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	82.5		1		19800		16030		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	72			1		19440		16040		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	72			1		21600		16050		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	72			1		23760		16060		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	72			1		25920		16070		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	57			1		22230		16080		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	57			1		23940		16090		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	57			1		25650		16100		0.005	2000	5000	6 		5		216"
P[${#P[@]}]="400	5			5			10.5	57			1		27360		16110		0.005	2000	5000	6 		5		216"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
#ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		n_small=$(echo "${P[$j]}" | awk '{print $1}') # number of small monomers
		epsilon_small=$(echo "${P[$j]}" | awk '{print $2}') #
		epsilon_big=$(echo "${P[$j]}" | awk '{print $3}')
		r=$(echo "${P[$j]}" | awk '{print $4}')
		lz=$(echo "${P[$j]}" | awk '{print $5}') # range of x, y or z from -l to l
		sig3=$(echo "${P[$j]}" | awk '{print $6}') # crowder size
		n_crowd=$(echo "${P[$j]}" | awk '{print $7}')
		randseed_group=$(echo "${P[$j]}" | awk '{print $8}')
		run_dt=$(echo "${P[$j]}" | awk '{print $9}')
		bug_dump=$(echo "${P[$j]}" | awk '{print $10}')
		all_dump=$(echo "${P[$j]}" | awk '{print $11}')
		sig2=$(echo "${P[$j]}" | awk '{print $12}') # big monomer size
		n_big=$(echo "${P[$j]}" | awk '{print $13}') # number of big monomers
		m_big=$(echo "${P[$j]}" | awk '{print $14}')
		dirname=epss${epsilon_small}epsl${epsilon_big}r${r}al${sig2}nl${n_big}ml${m_big}ns${n_small}ac${sig3}nc${n_crowd}lz${lz}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}.ring

		# modifying foci-ring-cylinder-init_config_minimize-harmonic.lmp
		cp foci-ring-cylinder-init_config_minimize-harmonic.lmp input.lmp
		
		echo "variable n_big equal ${n_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_small equal ${epsilon_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_big equal ${epsilon_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
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

		# rename the input.lmp to min_init_config.lmp
		mv input.lmp minimize_initial_config.lmp

		if [ $((n_crowd+n_small+n_big)) -ge 25000 ]; then
			cp foci-ring-cylinder-ac_equal-cores_equal_more_8.lmp input.lmp
		else
			cp foci-ring-cylinder-ac_equal-cores_less_8.lmp input.lmp
		fi

		echo "variable n_big equal ${n_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_small equal ${epsilon_small}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon_big equal ${epsilon_big}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
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

		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv minimize_initial_config.lmp "${dirname}"
		mv input.lmp "${dirname}"
		cp ./*.data initial_config.data
		mv initial_config.data "${dirname}"
		cp submit_minimize_initial_config.sh "${dirname}/" 
		
		N=$((n_big+n_small))
		if [ $((n_crowd+N)) -ne ${N} ]; 
		then
			if [ $((n_crowd+N)) -le 5000 ]; 
			then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh "$dirname"
			elif [ $((n_crowd+N)) -gt 5000 ] && [ $((n_crowd+N)) -le 15000 ]; then
				cp submit_cpu4_gt5000_le15000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 15000 ] && [ $((n_crowd+N)) -le 25000 ]; then
				cp submit_cpu8_gt15000_le25000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N)) -gt 25000 ] && [ $((n_crowd+N)) -le 75000 ]; then
				cp submit_cpu16_gt25000_le75000.sh submit.sh && mv submit.sh "${dirname}"
			else
				cp submit_cpu32_gt75000_le200000.sh submit.sh && mv submit.sh "${dirname}"
			fi
		else
			cp submit_nocrowd_8hrs.sh submit.sh && mv submit.sh "${dirname}"
		fi
		cd "${dirname}" || exit
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
