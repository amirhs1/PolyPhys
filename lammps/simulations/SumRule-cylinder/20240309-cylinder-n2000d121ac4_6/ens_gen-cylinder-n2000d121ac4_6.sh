#!/bin/bash
# The data file generated at the end of simulation designed below
# 20231101-cylinder-n2000d61ac6-5times18pow8-with_minimize-dt0.01-no_processors
# data file original name:
# N2000epsilon5.0r60.5lz144sig4.0nc0dt0.01bdump5000adump10000ens1.bug

#			 N		epsilon1	r		lz		sig2	n_crowd		randseed_group	run_dt 	bug_dump 	all_dump
P[${#P[@]}]="2000	5.0			60.5	144	4.0		0 			7000			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		9083 		7010			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		13624 		7020			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		18166 		7030			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		20436 		7040			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		22707 		7050			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		24978 		7060			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		27248 		7070			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		29519		7080			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		31790 		7090			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		34061 		7100			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	4.0		36331 		7110			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		2599 		7120			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		3899 		7130			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		5198 		7140			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		5848 		7150			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		6498 		7160			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		7148 		7170			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		7798 		7180			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		8447 		7190			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		9097 		7200			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		9747 		7210			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			60.5	144	6.0		10397 		7220			0.005	5000		10000"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

# Reset the counter for moving files to folders

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
		cp cylinder-system_minimize.lmp input.lmp
		echo "variable n_bug equal $N" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable randseed equal $(($i+${randseed_group}))" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp
		mv input.lmp system_minimize.lmp
		[[ -d "${dirname}" ]] || mkdir "${dirname}"
		mv system_minimize.lmp "${dirname}/"
		cp submit_system_minimize.sh "${dirname}/"
		cp chain.2000.data "${dirname}/chain.2000.data"
		
		# production script
		cp cylinder_v9_ac_larger_bin0.3.lmp input.lmp
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

		# Loop through each folder and move the corresponding file
		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu2_le5000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 10000 ]; then
				cp submit_cpu4_gt5000_le10000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 10000 ] && [ $((${n_crowd}+${N})) -le 20000 ]; then
				cp submit_cpu8_gt10000_le20000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 20000 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
				cp submit_cpu16_gt20000_le50000.sh submit.sh && mv submit.sh ${dirname}	
			elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
				cp submit_cpu32_gt50000_le100000.sh submit.sh && mv submit.sh ${dirname}	
			else
				cp submit_cpu32_gt100000_le160000.sh submit.sh && mv submit.sh ${dirname}	
			fi
		else
			cp submit_cpu1_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
	

