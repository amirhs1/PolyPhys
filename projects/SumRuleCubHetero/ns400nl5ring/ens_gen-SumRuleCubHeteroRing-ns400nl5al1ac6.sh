#!/bin/bash
# Generated on date 20240726
# ns lcube ac nc randseed_group dt b_dump a_dump al nl ml
P[${#P[@]}]="400	45	6	645	61000	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	967	61010	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	1290	61020	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	1451	61030	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	1612	61040	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	1773	61050	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	1934	61060	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	2095	61070	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	2257	61080	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	2418	61090	0.005	2000	5000	1	5	1"
P[${#P[@]}]="400	45	6	2579	61100	0.005	2000	5000	1	5	1"

#ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
ens=1
# ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
for i in ${ens}; do 
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
        simname=al${sig2}nl${n_big}ml${m_big}ns${n_small}ac${sig3}nc${n_crowd}l${l}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}
		
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
        
        cp foci-ring-cubic-ac_equal-neigh_multi.lmp input.lmp

		echo "variable simname string ${simname}" | cat - input.lmp > temp && mv temp input.lmp
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
		cp ./submit_init_config_minimize.sh "${dirname}"

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
			elif [ $((n_crowd+N)) -gt 200000 ] && [ $((n_crowd+N)) -le 250000 ]; then
				cp submit_cpu32_gt200000_le250000.sh submit.sh && mv submit.sh "${dirname}"
			else
				cp submit_cpu32_gt250000_le300000.sh submit.sh && mv submit.sh "${dirname}"
			fi
		else
			cp submit_nocrowd_8hrs.sh submit.sh && mv submit.sh "${dirname}"
		fi
		cd "${dirname}" || exit
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done

		