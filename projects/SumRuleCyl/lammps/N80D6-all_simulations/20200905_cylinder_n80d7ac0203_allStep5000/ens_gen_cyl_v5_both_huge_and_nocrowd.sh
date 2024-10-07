#!/bin/bash
# randseed_group is the mother of the seeds for an ensemble with same input parameters
#			 N	epsilon1	r	lz	sig2	n_crowd	randseed_group
P[${#P[@]}]="80	5.0	3.5 45.0	0.2	24300	1600"
P[${#P[@]}]="80	5.0	3.5	40.5	0.2	54675	1610"
P[${#P[@]}]="80	5.0	3.5	40.5	0.2	65610	1620"
P[${#P[@]}]="80	5.0	3.5	16.5	0.2	49005	1650"
P[${#P[@]}]="80	5.0	3.5	16.5	0.2	57915	1660"
P[${#P[@]}]="80	5.0	3.5	16.5	0.2	62370	1670"
P[${#P[@]}]="80	5.0	3.5 16.5	0.2	66825	1680"

P[${#P[@]}]="80	5.0	3.5 45.0	0.3	18000	1690"
P[${#P[@]}]="80	5.0	3.5	45.0	0.3	32400	1700"
P[${#P[@]}]="80	5.0	3.5	38.5	0.3	36000	1710"
P[${#P[@]}]="80	5.0	3.5	27.0	0.3	40040	1720"
P[${#P[@]}]="80	5.0	3.5	27.0	0.3	32400	1730"
P[${#P[@]}]="80	5.0	3.5	27.0	0.3	36720	1740"
P[${#P[@]}]="80	5.0	3.5	27.0	0.3	41040	1750"
P[${#P[@]}]="80	5.0	3.5	27.0	0.3	43200	1760"



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
		cp cylinder_v5_ac_smaller_unscale_dumps_rand_seed_changed_write_data.lmp input.lmp

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

		echo "variable n_bug equal ${N}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
		echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
		echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

		[[ -d ${dirname} ]] || mkdir ${dirname}
		mv input.lmp ${dirname}
		cp chain.80.data ${dirname}
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
