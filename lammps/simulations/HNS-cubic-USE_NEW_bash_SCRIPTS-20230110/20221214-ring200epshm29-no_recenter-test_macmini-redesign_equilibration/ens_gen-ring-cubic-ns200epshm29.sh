#!/bin/bash
#			 N		epsilon_hns_mon	kbend	n_hns	a_crowd	n_crowd		l	rseed_group	run_dt 	ndump 	adump 	
#P[${#P[@]}]="200	29				2		0		2		0			25	1000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				2		12		2		0			25	10000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				2		24		2		0			25	13000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				2		36		2		0			25	2000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				2		48		2		0			25	3000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				10		0		2		0			50	4000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				10		12		2		0			50	11000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				10		24		2		0			50	14000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				10		36		2		0			50	5000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				10		48		2		0			50	6000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				20		0		2		0			60	7000		0.005	2000 	5000"
P[${#P[@]}]="200	29				20		12		2		1000		60	12000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				20		24		2		0			60	15000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				20		36		2		0			60	8000		0.005	2000 	5000"
#P[${#P[@]}]="200	29				20		48		2		0			60	9000		0.005	2000 	5000"

#ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=$(echo "${P[$j]}" | awk '{print $1}') # number of monomers
		eps_hm=$(echo "${P[$j]}" | awk '{print $2}')
		kBend=$(echo "${P[$j]}" | awk '{print $3}')
		n_hns=$(echo "${P[$j]}" | awk '{print $4}')
		sig4=$(echo "${P[$j]}" | awk '{print $5}') # sizeo of crowders
		n_crowd=$(echo "${P[$j]}" | awk '{print $6}')
		l=$(echo "${P[$j]}" | awk '{print $7}') # range of x, y or z from -l to l
		randseed_group=$(echo "${P[$j]}" | awk '{print $8}')
		run_dt=$(echo "${P[$j]}" | awk '{print $9}') # integration time delta 
		nucleoid_dump=$(echo "${P[$j]}" | awk '{print $10}')
		all_dump=$(echo "${P[$j]}" | awk '{print $11}')
		dirname=N${N}epshm${eps_hm}kbmm${kBend}nh${n_hns}ac${sig4}nc${n_crowd}l${l}dt${run_dt}ndump${nucleoid_dump}adump${all_dump}ens${i}.ring
		# minimization and equilibration script
		cp hns-cubic-init_config_minimize.lmp input.lmp
		echo "variable adump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable ndump equal ${nucleoid_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig4 equal ${sig4}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_hns equal ${n_hns}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable kBend11 equal ${kBend}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable eps_hm equal ${eps_hm}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable N equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp minimize_initial_config.lmp
		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv minimize_initial_config.lmp "${dirname}"
		# production script
		if [ $((n_crowd)) -eq 0 ];
		then
			cp hns-cubic_v8-ac_equal_a.lmp input.lmp
		else
			if [ $((sig4)) -ne 1 ];
			then
				cp hns-cubic_v8.lmp input.lmp	
			else
				cp hns-cubic_v8-ac_equal_a.lmp input.lmp
			fi
		fi
		echo "variable adump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable ndump equal ${nucleoid_dump}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig4 equal ${sig4}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_hns equal ${n_hns}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable kBend11 equal ${kBend}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable eps_hm equal ${eps_hm}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable N equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv input.lmp "${dirname}"
	
		cp ./dna.data ./"${dirname}/"
		cp ./hns.mol ./"${dirname}/"

		N_no_crowd=$((N+n_hns))
		if [ $((n_crowd+N_no_crowd)) -ne ${N_no_crowd} ];
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
				cp submit_cpu32_gt150000_le200000.sh submit.sh && mv submit.sh "${dirname}"
			else
				cp submit_cpu32_gt250000_le300000.sh submit.sh && mv submit.sh "${dirname}"
			fi
		else
			cp submit_nocrowd_20hrs.sh submit.sh && mv submit.sh "${dirname}"
		fi
		cd "${dirname}" || exit
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
