runs_stat_20200616.txt#!/bin/bash
# randseed_group is the mother of the seeds for an ensemble with same input parameters
#			 randseed_group
P[${#P[@]}]="1010"
P[${#P[@]}]="1010"
P[${#P[@]}]="1010"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		randseed_group=`echo ${P[$j]} | awk '{print $1}'`
		dirname=ens${i}
# The command I should use:
		cp restart.lmp input.lmp

		#if [ ${n_crowd} -eq 0  -a ${sig2} == 0.3 ]; then
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
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

		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp data.chain.80 ${dirname}
		cp submit_huge_nc_b_35000_le_70000.sh submit.sh && mv submit.sh ${dirname}
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..


#		echo "cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done
