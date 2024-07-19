#!/bin/bash
# Here, we use 8*23 equilbrated from a long long of a single chain in the
# cylindrical confinement and use them as input for simulation. The data file
# generated in MDAanalysis. The simulations crashes due to missed atom and fene
# bond too long. I cannot use this strategy.
#
# randseed_group is the mother of the seeds for an ensemble with same input parameters

#			 N		epsilon1	r		lz		sig2	n_crowd		randseed_group	run_dt 	bug_dump 	all_dump
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		0 			4000			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		3941 		4010			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		5911 		4020			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		7881 		4030			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		8866 		4040			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		9851 		4050			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		10836 		4060			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		11821 		4070			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		12807 		4080			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		13792 		4090			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		14777 		4100			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	4.0		15762 		4110			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		1168 		4120			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		1752 		4130			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		2335 		4140			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		2627 		4150			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		2919 		4160			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		3211 		4170			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		3503 		4180			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		3795 		4190			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		4087		4200			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		4379 		4210			0.005	5000		10000"
P[${#P[@]}]="2000	5.0			30.5	233.5	6.0		4670 		4220			0.005	5000		10000"

ens='1 2 3 4 5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

# Reset the counter for moving files to folders
counter=1

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

		# Loop through each folder and move the corresponding file
		if [ -d "$dirname" ]; then
			# Construct the file name
			file="file_${counter}.data"
			
			# Check if the file exists
			if [ -f "./data_files/$file" ]; then
				# Move the file to the folder
				cp "./data_files/$file" "chain.2000.data"
				mv "chain.2000.data" "${dirname}/"
			else
				echo "File $file does not exist"
			fi

			# Increment the counter
			((counter++))
    	fi


		if [ $((${n_crowd}+${N})) -ne ${N} ]; then
			if [ $((${n_crowd}+${N})) -le 5000 ]; then
				cp submit_cpu4_le5000.sh submit.sh && mv submit.sh ${dirname}
			elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 10000 ]; then
				cp submit_cpu4_gt5000_le10000.sh submit.sh && mv submit.sh ${dirname}	
			else
				cp submit_cpu8_gt10000_le20000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_cpu2_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..
	done
done
	

