#!/bin/bash
# kbend is in LAMMPS lingo which is half the regualr one, so kbend_mm=2 in LAMMPS means kbend_mm=4 in theory.
#            N		kbend   n_hns   a_crowd n_crowd	r	lz		rseed_group   
P[${#P[@]}]="1000	24		0		2		0		8	628.5	100000"
P[${#P[@]}]="1000	24		30		2		0		8	628.5	700000"
P[${#P[@]}]="1000	24		30		2		4243	8	628.5	700010"
P[${#P[@]}]="1000	24		30		2		21212	8	628.5	700090"
P[${#P[@]}]="1000	24		30		1		33939	8	628.5	700100"
P[${#P[@]}]="1000	24		30		1		169695	8	628.5	700180"

#ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
ens='1 2' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

# Equilibrated data file generated in "N1000kbmm24r8nh0ac2lz628.5nc0ens1.j03.ring.dnaMinimize" simulation are used:
data_dir="equilibrated_data_files-N1000kbmm24r8nh0ac2lz628.5nc0ens1.j03.ring.dnaMinimize"
first_step=250200000
delta_step=200000
counter=0  # Initialize a counter variable
for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=$(echo "${P[$j]}" | awk '{print $1}') # number of monomers
		kBend=$(echo "${P[$j]}" | awk '{print $2}')
		n_hns=$(echo "${P[$j]}" | awk '{print $3}')
		sig4=$(echo "${P[$j]}" | awk '{print $4}') # sizeo of crowders
		n_crowd=$(echo "${P[$j]}" | awk '{print $5}')
		r=$(echo "${P[$j]}" | awk '{print $6}')
		lz=$(echo "${P[$j]}" | awk '{print $7}') # range of x, y or z from -l to l
		randseed_group=$(echo "${P[$j]}" | awk '{print $8}')
		dir_name=N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}.ring
		sim_name=N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}		
		# production script
		[[ -d "${dir_name}" ]] || mkdir "${dir_name}"
		file_number=$((first_step + counter*delta_step))  # Calculate the file number
        file_name="time_step-${file_number}.data"  # Generate the file name
		cp "./${data_dir}/${file_name}" "./${dir_name}/minimized.dna.data"
		counter=$((counter + 1))
		cp hns.mol "${dir_name}/"
		if [ $((n_crowd)) -eq 0 ];
		then
			cp hns-cylindrical_v1-ac_equal_a.lmp input.lmp
		else
			if [ $((sig4)) -ne 1 ];
			then
				cp hns-cylindrical_v1.lmp input.lmp	
			else
				cp hns-cylindrical_v1-ac_equal_a.lmp input.lmp
			fi
		fi
		echo "variable simname string ${sim_name}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig4 equal ${sig4}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_hns equal ${n_hns}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable kBend11 equal ${kBend}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp "${dir_name}/"
		N_no_crowd=$((N+n_hns))
		if [ $((n_crowd+N_no_crowd)) -ne ${N_no_crowd} ];
		then
			if [ $((n_crowd+N)) -le 5000 ]; 
			then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh "$dir_name"
			elif [ $((n_crowd+N)) -gt 5000 ] && [ $((n_crowd+N)) -le 15000 ]; then
				cp submit_cpu4_gt5000_le15000.sh submit.sh && mv submit.sh "${dir_name}"
			elif [ $((n_crowd+N)) -gt 15000 ] && [ $((n_crowd+N)) -le 30000 ]; then
				cp submit_cpu8_gt15000_le30000.sh submit.sh && mv submit.sh "${dir_name}"
			elif [ $((n_crowd+N)) -gt 30000 ] && [ $((n_crowd+N)) -le 75000 ]; then
				cp submit_cpu16_gt30000_le75000.sh submit.sh && mv submit.sh "${dir_name}"
			elif [ $((n_crowd+N)) -gt 75000 ] && [ $((n_crowd+N)) -le 150000 ]; then
				cp submit_cpu32_gt75000_le150000.sh submit.sh && mv submit.sh "${dir_name}"
			elif [ $((n_crowd+N)) -gt 150000 ] && [ $((n_crowd+N)) -le 200000 ]; then
				cp submit_cpu32_gt150000_le200000.sh submit.sh && mv submit.sh "${dir_name}"
			elif [ $((n_crowd+N)) -gt 200000 ] && [ $((n_crowd+N)) -le 250000 ]; then
				cp submit_cpu32_gt150000_le200000.sh submit.sh && mv submit.sh "${dir_name}"
			else
				cp submit_cpu32_gt250000_le300000.sh submit.sh && mv submit.sh "${dir_name}"
			fi
		else
			cp submit_nocrowd_20hrs.sh submit.sh && mv submit.sh "${dir_name}"
		fi
		cd "${dir_name}" || exit
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done