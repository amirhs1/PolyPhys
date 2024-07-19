#!/bin/bash

#			 N	    eps1 	r	    lz	    sig2	n_c		r_seed	run_dt  bdump   adump
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    1263	300000	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    1895	300020	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    2526	300040	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    2842	300060	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    3157	300080	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    3473	300100	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    3789	300120	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    4104	300140	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    4420	300160	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    4736	300180	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    4.0	    5051	300200	0.005	5000	10000"

P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    0	    300220	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    10102	300440	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    15153	300240	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    20204	300260	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    22729	300280	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    25254	300300	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    27780	300320	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    30305	300340	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    32831	300360	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    35358	300380	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    37881	300400	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    2.0	    40407	300420	0.005	5000	10000"

P[${#P[@]}]="2000	5.0	    13.0	431	    1.0	    80813	300460	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	431	    1.0	    121219	300480	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	300	    1.0	    112500	300500	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	300	    1.0	    126563	300520	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	285.0	1.0	    133594	300540	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	259.0	1.0	    133547	300560	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	259.0	1.0	    145688	300580	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	265.0	1.0	    161485	300600	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	265.0	1.0	    173903	300620	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	265.0	1.0	    186329	300640	0.005	5000	10000"
P[${#P[@]}]="2000	5.0	    13.0	265.0	1.0	    198750	300660	0.005	5000	10000"

ens='1 2 3 4 5 6 7 8' 
for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		epsilon1=`echo ${P[$j]} | awk '{print $2}'`
		r=`echo ${P[$j]} | awk '{print $3}'`
		lz=`echo ${P[$j]} | awk '{print $4}'`
		sig2=`echo ${P[$j]} | awk '{print $5}'`
		n_crowd=`echo ${P[$j]} | awk '{print $6}'`
		randseed_group=`echo ${P[$j]} | awk '{print $7}'`
		run_dt=`echo ${P[$j]} | awk '{print $8}'`
		bug_dump=`echo ${P[$j]} | awk '{print $9}'`
		all_dump=`echo ${P[$j]} | awk '{print $10}'`
		dirname=N${N}epsilon${epsilon1}r${r}lz${lz}sig${sig2}nc${n_crowd}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}
		simname=N${N}epsilon${epsilon1}r${r}lz${lz}sig${sig2}nc${n_crowd}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}

		# production script
		cp continue.lmp input.lmp
		echo "variable simname string ${simname}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_bug equal $N" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
        echo "read_restart restart.after.71000000" | cat - input.lmp > temp && mv temp input.lmp 
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}

		# Loop through each folder and move the corresponding file
		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu2_le5000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 9000 ]; then
				cp submit_cpu4_gt5000_le9000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 9000 ] && [ $((${n_crowd}+${N})) -le 12500 ]; then
				cp submit_cpu4_gt9000_le12500.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 12500 ] && [ $((${n_crowd}+${N})) -le 22500 ]; then
				cp submit_cpu8_gt12500_le22500.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 22500 ] && [ $((${n_crowd}+${N})) -le 32500 ]; then
				cp submit_cpu16_gt22500_le32500.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 32500 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
				cp submit_cpu16_gt32500_le50000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 75000 ]; then
				cp submit_cpu32_gt50000_le75000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 75000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
				cp submit_cpu32_gt75000_le100000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 100000 ] && [ $((${n_crowd}+${N})) -le 130000 ]; then
				cp submit_cpu32_gt100000_le130000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 130000 ] && [ $((${n_crowd}+${N})) -le 150000 ]; then
				cp submit_cpu32_gt130000_le150000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 150000 ] && [ $((${n_crowd}+${N})) -le 200000 ]; then
				cp submit_cpu32_gt150000_le200000.sh submit.sh && mv submit.sh ${dirname}		
			else
				cp submit_cpu32_gt200000_le230000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_cpu2_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done