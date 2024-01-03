#!/bin/bash

#			 N	    eps1 	r	    lz	    sig2	n_c		r_seed	run_dt  bdump   adump
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    0	    400000	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    20824	400020	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    31236	400040	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    41648	400060	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    46854	400080	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    52060	400100	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    57266	400120	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    62472	400140	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    67678	400160	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    72884	400180	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    78090	400200	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	1.0	    83295	400220	0.005	5000	10000"

P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    2603	400240	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    3905	400260	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    5206	400280	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    5857	400300	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    6508	400320	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    7159	400340	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    7809	400360	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    8460	400380	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    9111	400400	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    9762	400420	0.005	5000	10000"
P[${#P[@]}]="1000	5.0	    8.0	    308.5	2.0	    10412	400440	0.005	5000	10000"

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
			elif [ $((${n_crowd}+${N})) -gt 22500 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
				cp submit_cpu16_gt22500_le50000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
				cp submit_cpu32_gt50000_le100000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 100000 ] && [ $((${n_crowd}+${N})) -le 130000 ]; then
				cp submit_cpu32_gt100000_le130000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 130000 ] && [ $((${n_crowd}+${N})) -le 150000 ]; then
				cp submit_cpu32_gt130000_le150000.sh submit.sh && mv submit.sh ${dirname}	
			else
				cp submit_cpu32_gt150000_le201000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
            if [ $((${n_crowd}+${N})) -le 1000 ]; then
                cp submit_cpu1_nocrowd.sh submit.sh && mv submit.sh ${dirname}
            else
                cp submit_cpu2_nocrowd.sh submit.sh && mv submit.sh ${dirname}
            fi
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done