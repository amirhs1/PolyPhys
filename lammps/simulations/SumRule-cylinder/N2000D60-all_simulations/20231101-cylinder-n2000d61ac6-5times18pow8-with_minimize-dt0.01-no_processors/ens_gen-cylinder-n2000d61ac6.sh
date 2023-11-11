#!/bin/bash
# randseed_group is the mother of the seeds for an ensemble with same input parameters

#			 N		epsilon1	r		lz		sig2	n_crowd		randseed_group	run_dt 	bug_dump 	all_dump
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		0 			100				0.005	5000		10000"

ens='1' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
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
		simname=N${N}epsilon${epsilon1}r${r}lz${lz}sig${sig2}nc${n_crowd}dt${run_dt}bdump${bug_dump}adump${all_dump}ens${i}
		# minimization and equilibration script
		cp cylinder-minimize_chain.lmp input.lmp
		echo "variable n_bug equal $N" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp minimize_chain.lmp
		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv minimize_chain.lmp "${dirname}/"
		cp submit_minimize_chain.sh "${dirname}/"
		
		# production script
		cp cylinder_v8_ac_larger_bin0.3.lmp input.lmp
		echo "variable simname string ${simname}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
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
		cp chain.2000.data ${dirname}

		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_cpu8_gt5000_le20000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
	