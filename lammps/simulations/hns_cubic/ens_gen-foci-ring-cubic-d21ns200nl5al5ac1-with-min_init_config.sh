#!/bin/bash
#			 N		epsilon_hns_mon	n_hns	a_crowd	n_crowd		m_crowd	l	randseed_group	run_dt 	n_dump 	a_dump 	
P[${#P[@]}]="200	29				36		2		5969		8		25	1000		0.005	2000 	5000"

#ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
ens='2' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=$(echo "${P[$j]}" | awk '{print $1}') # number of monomers
		eps_hm=$(echo "${P[$j]}" | awk '{print $2}')
		n_hns=$(echo "${P[$j]}" | awk '{print $3}')
		sig4=$(echo "${P[$j]}" | awk '{print $4}') # sizeo of crowders
		n_crowd=$(echo "${P[$j]}" | awk '{print $5}')
		m_crowd=$(echo "${P[$j]}" | awk '{print $6}')
		l=$(echo "${P[$j]}" | awk '{print $7}') # range of x, y or z from -l to l
		randseed_group=$(echo "${P[$j]}" | awk '{print $8}')
		run_dt=$(echo "${P[$j]}" | awk '{print $9}') # integration time delta 
		nucleoid_dump=$(echo "${P[$j]}" | awk '{print $10}')
		all_dump=$(echo "${P[$j]}" | awk '{print $11}')
		dirname=N${N}epshm${eps_hm}nh${n_hns}ac${sig4}nc${n_crowd}mc${m_crowd}l${l}dt${run_dt}ndump${nucleoid_dump}adump${all_dump}ens${i}
		# modifying foci-ring-cylinder-init_config_minimize-harmonic.lmp
		cp hns-cubic_v4.lmp input.lmp
		echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable nucleoid_dump equal ${nucleoid_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable m_crowd equal ${m_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig4 equal ${sig4}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_hns equal ${n_hns}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable eps_hm equal ${eps_hm}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable N equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv input.lmp "${dirname}"
	
		cp ./nucleoid_crowders.data "./${dirname}/"

		N_no_crowd=$((N+n_hns))
		if [ $((n_crowd+N_no_crowd)) -ne ${N_no_crowd} ]; 
		then
			if [ $((n_crowd+N)) -le 5000 ]; 
			then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh "$dirname"
			elif [ $((n_crowd+N_no_crowd)) -gt 5000 ] && [ $((n_crowd+N_no_crowd)) -le 15000 ]; then
				cp submit_cpu4_gt5000_le15000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N_no_crowd)) -gt 15000 ] && [ $((n_crowd+N_no_crowd)) -le 30000 ]; then
				cp submit_cpu8_gt15000_le30000.sh submit.sh && mv submit.sh "${dirname}"
			elif [ $((n_crowd+N_no_crowd)) -gt 30000 ] && [ $((n_crowd+N_no_crowd)) -le 75000 ]; then
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
