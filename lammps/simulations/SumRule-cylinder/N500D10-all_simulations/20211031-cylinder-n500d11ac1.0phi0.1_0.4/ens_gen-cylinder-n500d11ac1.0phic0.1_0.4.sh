#!/bin/bash
# randseed_group is the mother of the seeds for an ensemble with same input parameters

#			 N	epsilon1	r	lz	sig2	n_crowd		randseed_group	run_dt bug_dump all_dump

P[${#P[@]}]="500	5.0	5.5	205.5	1.0	6150	6230	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	9225	6240	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	12300	6250	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	13838	6260	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	15375	6270	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	16913	6280	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	18450	6290	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	19988	6300	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	21525	6310	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	23063	6320	0.005	1000	5000"
P[${#P[@]}]="500	5.0	5.5	205.5	1.0	24600	6330	0.005	1000	5000"
ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

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
# The command I should use:
		
		if echo "${sig2} < 1.0" | bc -l | grep -q 1; then
			if [ $((${n_crowd}+${N})) -le 25000 ]; then
				cp cylinder_v8_ac_smaller.lmp input.lmp
			else
				cp cylinder_v8_ac_smaller_cpu8.lmp input.lmp # processors 2 2 *
			fi
		else
			if [ $((${n_crowd}+${N})) -le 25000 ]; then
				cp cylinder_v8_ac_larger.lmp input.lmp
			else
				cp cylinder_v8_ac_larger_cpu8.lmp input.lmp # processors 2 2 *
			fi
		fi
		

		#if [ ${n_crowd} -eq 0  -a ${sig2} == 0.3 ]; then
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		#elif [ ${n_crowd} -eq 20400  -a ${sig2} ==  0.3 ]; then
		 #       echo 'variable randseed equal ${i}+1120' | cat - input.lmp > temp && mv temp input.lmp
		#elif [ ${n_crowd} -eq 61200  -a ${sig2} == 0.3 ]; then
		#        echo 'variable randseed equal ${i}+1130' | cat - input.lmp > temp && mv temp input.lmp
		#elif [ ${n_crowd} -eq 0  -a ${sig2} == 0.2 ]; then
		#        echo 'variable randseed equal ${i}+1140' | cat - input.lmp > temp && mv temp input.lmp
		#else [ ${n_crowd} -eq 68850  -a ${sig2} == 0.2 ]
		#        echo 'variable randseed equal ${i}+1150' | cat - input.lmp > temp && mv temp input.lmp
		#fi

		echo "variable n_bug equal $N" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp chain.$N.data ${dirname}

		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu2_le5000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 25000 ]; then
				cp submit_cpu4_gt5000_le25000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 25000 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
				cp submit_cpu8_gt25000_le50000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
				cp submit_cpu16_gt50000_le100000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_cpu32_gt100000_le200000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
	