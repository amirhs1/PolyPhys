#!/bin/bash
#			n_small	eps_small 	eps_large	r		lz		a_crowd	n_crowd	randseed	run_dt 	b_dump 	a_dump 	a_large n_large m_large
#P[${#P[@]}]="400	5			5			10.5	65.5	1		0		18000		0.005	2000	5000	1 		5		1"
#P[${#P[@]}]="400	5			5			10.5	65.5	1		1965	18010		0.005	2000	5000	1 		5		1"
#P[${#P[@]}]="400	5			5			10.5	65.5	1		3930	18020		0.005	2000	5000	1 		5		1"
#P[${#P[@]}]="400	5			5			10.5	65.5	1		5895	18030		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		7860	18040		0.005	2000	5000	1 		5		1"
#P[${#P[@]}]="400	5			5			10.5	65.5	1		9825	18050		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		11790	18060		0.005	2000	5000	1 		5		1"
#P[${#P[@]}]="400	5			5			10.5	65.5	1		13755	18070		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		15720	18080		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		17685	18090		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		19650	18100		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		21615	18110		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		23580	18120		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		25545	18130		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		27510	18140		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		29475	18150		0.005	2000	5000	1 		5		1"
P[${#P[@]}]="400	5			5			10.5	65.5	1		31440	18160		0.005	2000	5000	1 		5		1"

ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
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

		mv input.lmp "${dirname}"
		cp minimized.initial.config.data "${dirname}"

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
