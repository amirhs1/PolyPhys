#!/bin/bash

#			 N		lz	sig2	n_crowd	epsilon1	r
#P[${#P[@]}]="80		51.0	0.3		0	5.0	3.0"
P[${#P[@]}]="80		51.0	0.3		61200	5.0	3.0"

ens='5 6 7 8' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		lz=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		r=`echo ${P[$j]} | awk '{print $6}'`
		dirname=N${N}nc${n_crowd}ens${i}lz${lz}r${r}epsilon${espsilon1}sig${sig2}
# The command I should use:
		cp cylinder.lmp input.lmp
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
		cp submit_huge.sh submit.sh && mv submit.sh ${dirname}
		cd ${dirname}
		[[ -d archive ]] || mkdir archive
		cd ..
#		echo "cd ${dirname} && mkdir archive && echo Hello World!  > log.out.${CURRENT_DATE}.${dirname}.txt" >> glost_input_${CURRENT_DATE}.txt
	done
done
