# require more sig2
# r=4.0 (D=7) is quite big for N=80	 chain
#			 N		r	sig2	n_crowd	epsilon1	lz
P[${#P[@]}]="80		4.0	0.3		0		5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		2500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		5000	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		7500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		10000	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		12500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		15000	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		17500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		20000	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		22500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		25000	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		27500	5.0			14.0"
P[${#P[@]}]="80		4.0	0.3		30000	5.0			14.0"

P[${#P[@]}]="80		3.0	0.3		0		5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		1250	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		2500	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		3750	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		5000	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		6250	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		7500	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		8750	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		10000	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		11250	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		12500	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		13750	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		15000	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		16250	5.0			18.0"
P[${#P[@]}]="80		3.0	0.3		17500	5.0			18.0"

ens='1 2 3 4 5 6 7 8'
exe='./lammps'

for i in ${ens}; do
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		r=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		lz=`echo ${P[$j]} | awk '{print $6}'`

		dirname=N${N}r${r}sig${sig2}nc${n_crowd}
		[[ -d ${dirname} ]] || mkdir ${dirname}
		cd ${dirname}
		[[ -d archive ]] || mkdir archive

		# Make script file and submit
		# This part should be modified: maybe no sqsub in kisti HPC
		echo 	sqsub -q mpi -n 4 -r 7d --mpp 3G -o out${i}.txt ${exe} \
					-i '../in.chain1' \
					-var in_filename "../data.chain.${N}" \
					-var r ${r} \
					-var lz ${lz} \
					-var sig2 ${sig2} \
					-var n_crowd ${n_crowd} \
					-var epsilon1 ${epsilon1} \
					-var i ${i} > in${i}.sh

		bash "in${i}.sh"

		#jobID=`bash in${i}.sh 2>&1`
		#echo ${jobID}
		#jobID=`echo ${jobID} | cut -d ' ' -f4`
		#echo ${CLUSTER} ${jobID} > in${i}.txt
		cd ..
	done
done
