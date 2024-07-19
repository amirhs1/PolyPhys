#!/bin/bash
# randseed_group is the mother of the seeds for an ensemble with same input parameters
#			 N	epsilon1	r	lz	sig2	n_crowd	randseed_group
P[${#P[@]}]="80	5.0	3.5	45.0	0.3	0	1090"
P[${#P[@]}]="80	5.0	3.5	45.0	0.3	13307	1160"
P[${#P[@]}]="80	5.0	3.5	45.0	0.3	26412	1170"
P[${#P[@]}]="80	5.0	3.5	38.5	0.3	33981	1180"
P[${#P[@]}]="80	5.0	3.5	36.6	0.3	40162	1190"
P[${#P[@]}]="80	5.0	3.5	27.6	0.3	34960	1200"
P[${#P[@]}]="80	5.0	3.5	18.0	0.3	25404	1210"

P[${#P[@]}]="80	5.0	3.5	45.0	0.2	0	1220"
P[${#P[@]}]="80	5.0	3.5	45.0	0.2	51840	1230"
P[${#P[@]}]="80	5.0	3.5	40.5	0.2	77760	1240"
P[${#P[@]}]="80	5.0	3.5	29.25	0.2	78624	1250"
P[${#P[@]}]="80	5.0	3.5	13.5	0.2	41472	1260"

ens='1 2 3 4' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
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
		dirname=N${N}epsilon${epsilon1}r${r}lz${lz}sig${sig2}nc${n_crowd}ens${i}
# The command I should use:
		cp cylinder_unscale_dumps_rand_seed_changed.lmp input.lmp

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

		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp data.chain.80 ${dirname}
		if [ ${n_crowd} -ne 0 ]; then
			if [ ${n_crowd} -le 35000 ]; then
				cp submit_huge_nc_le_35000.sh submit.sh && mv submit.sh ${dirname}
			elif [ ${n_crowd} -gt 35000  ] && [ ${n_crowd} -le 70000  ]; then
				cp submit_huge_nc_b_35000_le_70000.sh submit.sh && mv submit.sh ${dirname}
			else
				cp submit_huge_nc_b_70000.sh submit.sh && mv submit.sh ${dirname}
			fi
		else
			cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
		fi
		cd ${dirname}
		[[ -d restarts ]] || mkdir restarts
		cd ..


#		echo "cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done
