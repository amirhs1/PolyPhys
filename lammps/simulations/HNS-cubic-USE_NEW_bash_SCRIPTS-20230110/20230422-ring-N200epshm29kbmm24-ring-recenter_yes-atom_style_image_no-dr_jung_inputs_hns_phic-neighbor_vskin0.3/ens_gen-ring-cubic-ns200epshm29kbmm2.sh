#!/bin/bash
# kbend is in LAMMPS lingo which is half the regualr one, so kbend_mm=2 in LAMMPS means kbend_mm=4 in theory.
#            N		kbend   n_hns   a_crowd n_crowd     l	rseed_group 
P[${#P[@]}]="200 	2		0		1		0			25	50000"
P[${#P[@]}]="200 	2		0		1		19099		25	50010"
P[${#P[@]}]="200 	2		0		1		28648		25	50020"
P[${#P[@]}]="200 	2		0		1		38198		25	50030"
P[${#P[@]}]="200 	2		0		1		47747		25	50040"
P[${#P[@]}]="200 	2		0		1		57296		25	50050"
P[${#P[@]}]="200 	2		0		1		66846		25	50060"
P[${#P[@]}]="200 	2		0		1		76395		25	50070"
P[${#P[@]}]="200 	2		0		1		85944		25	50080"
P[${#P[@]}]="200 	2		0		1		95493		25	50090"
P[${#P[@]}]="200 	2		4		1		0			25	51000"
P[${#P[@]}]="200 	2		4		1		19099		25	51010"
P[${#P[@]}]="200 	2		4		1		28648		25	51020"
P[${#P[@]}]="200 	2		4		1		38198		25	51030"
P[${#P[@]}]="200 	2		4		1		47747		25	51040"
P[${#P[@]}]="200 	2		4		1		57296		25	51050"
P[${#P[@]}]="200 	2		4		1		66846		25	51060"
P[${#P[@]}]="200 	2		4		1		76395		25	51070"
P[${#P[@]}]="200 	2		4		1		85944		25	51080"
P[${#P[@]}]="200 	2		4		1		95493		25	51090"
P[${#P[@]}]="200 	2		8		1		0			25	52000"
P[${#P[@]}]="200 	2		8		1		19099		25	52010"
P[${#P[@]}]="200 	2		8		1		28648		25	52020"
P[${#P[@]}]="200 	2		8		1		38198		25	52030"
P[${#P[@]}]="200 	2		8		1		47747		25	52040"
P[${#P[@]}]="200 	2		8		1		57296		25	52050"
P[${#P[@]}]="200 	2		8		1		66846		25	52060"
P[${#P[@]}]="200 	2		8		1		76395		25	52070"
P[${#P[@]}]="200 	2		8		1		85944		25	52080"
P[${#P[@]}]="200 	2		8		1		95493		25	52090"
P[${#P[@]}]="200 	2		12		1		0			25	53000"
P[${#P[@]}]="200 	2		12		1		19099		25	53010"
P[${#P[@]}]="200 	2		12		1		28648		25	53020"
P[${#P[@]}]="200 	2		12		1		38198		25	53030"
P[${#P[@]}]="200 	2		12		1		47747		25	53040"
P[${#P[@]}]="200 	2		12		1		57296		25	53050"
P[${#P[@]}]="200 	2		12		1		66846		25	53060"
P[${#P[@]}]="200 	2		12		1		76395		25	53070"
P[${#P[@]}]="200 	2		12		1		85944		25	53080"
P[${#P[@]}]="200 	2		12		1		95493		25	53090"
P[${#P[@]}]="200 	2		16		1		0			25	54000"
P[${#P[@]}]="200 	2		16		1		19099		25	54010"
P[${#P[@]}]="200 	2		16		1		28648		25	54020"
P[${#P[@]}]="200 	2		16		1		38198		25	54030"
P[${#P[@]}]="200 	2		16		1		47747		25	54040"
P[${#P[@]}]="200 	2		16		1		57296		25	54050"
P[${#P[@]}]="200 	2		16		1		66846		25	54060"
P[${#P[@]}]="200 	2		16		1		76395		25	54070"
P[${#P[@]}]="200 	2		16		1		85944		25	54080"
P[${#P[@]}]="200 	2		16		1		95493		25	54090"
P[${#P[@]}]="200 	2		20		1		0			25	55000"
P[${#P[@]}]="200 	2		20		1		19099		25	55010"
P[${#P[@]}]="200 	2		20		1		28648		25	55020"
P[${#P[@]}]="200 	2		20		1		38198		25	55030"
P[${#P[@]}]="200 	2		20		1		47747		25	55040"
P[${#P[@]}]="200 	2		20		1		57296		25	55050"
P[${#P[@]}]="200 	2		20		1		66846		25	55060"
P[${#P[@]}]="200 	2		20		1		76395		25	55070"
P[${#P[@]}]="200 	2		20		1		85944		25	55080"
P[${#P[@]}]="200 	2		20		1		95493		25	55090"
P[${#P[@]}]="200 	2		30		1		0			25	56000"
P[${#P[@]}]="200 	2		30		1		19099		25	56010"
P[${#P[@]}]="200 	2		30		1		28648		25	56020"
P[${#P[@]}]="200 	2		30		1		38198		25	56030"
P[${#P[@]}]="200 	2		30		1		47747		25	56040"
P[${#P[@]}]="200 	2		30		1		57296		25	56050"
P[${#P[@]}]="200 	2		30		1		66846		25	56060"
P[${#P[@]}]="200 	2		30		1		76395		25	56070"
P[${#P[@]}]="200 	2		30		1		85944		25	56080"
P[${#P[@]}]="200 	2		30		1		95493		25	56090"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

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
		dirname=N${N}kbmm${kBend}nh${n_hns}ac${sig4}l${l}nc${n_crowd}ens${i}.ring
		simname=N${N}kbmm${kBend}nh${n_hns}ac${sig4}l${l}nc${n_crowd}ens${i}
		# minimization and equilibration script
		cp hns-cubic-minimize_dna.lmp input.lmp
		echo "variable randseed equal $((i+randseed_group))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable l equal ${l}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable kBend11 equal ${kBend}" | cat - input.lmp > temp && mv temp input.lmp
		echo '# Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp minimize_dna.lmp
		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv minimize_dna.lmp "${dirname}/"
		cp submit_minimize_dna.sh "${dirname}/" 
		# production script
		if [ $((n_crowd)) -eq 0 ];
		then
			cp hns-cubic_v9-ac_equal_a.lmp input.lmp
		else
			if [ $((sig4)) -ne 1 ];
			then
				cp hns-cubic_v9.lmp input.lmp	
			else
				cp hns-cubic_v9-ac_equal_a.lmp input.lmp
			fi
		fi
		echo "variable simname string ${simname}" | cat - input.lmp > temp && mv temp input.lmp
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

		mv input.lmp "${dirname}/"
		cp dna.data "${dirname}/"
		cp hns.mol "${dirname}/"

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