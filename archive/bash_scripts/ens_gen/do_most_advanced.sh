#!/bin/bash

# Load the modules
#module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps-omp/20170811
#echo "Starting run at: 'date'"

#			 N		l	sig2	n_crowd	epsilon1	
P[${#P[@]}]="50		36	0.3		0		5.0"

#P[${#P[@]}]="80		45	0.3		0		5.0"
#P[${#P[@]}]="80		45	0.3		2500	5.0"
#P[${#P[@]}]="80		45	0.3		5000	5.0"
#P[${#P[@]}]="80		45	0.3		7500	5.0"
#P[${#P[@]}]="80		45	0.3		10000	5.0"
#P[${#P[@]}]="80		45	0.3		12500	5.0"
#P[${#P[@]}]="80		45	0.3		15000	5.0"
#P[${#P[@]}]="80		45	0.3		17500	5.0"
#P[${#P[@]}]="80		45	0.3		20000	5.0"
#P[${#P[@]}]="80		45	0.3		22500	5.0"
#P[${#P[@]}]="80		45	0.3		25000	5.0"
#P[${#P[@]}]="80		45	0.3		27500	5.0"

#P[${#P[@]}]="2000	315	0.3		0		5.0"

ens='1 2' # Here, we define a global variable called "ens" which is the total ensemble we use in our simulation.

# Next 3 lines modified for Graham on SharCNet/ComputeCanada
lmp_exe='lmp_icc_openmpi' # This is also another global variable which contains the name of module/routine/program performing simulation. In my case, I will use the LAMMPS molecular dynamics (MD) package for simulations
CURRENT_DATE="$(date +%Y%m%d)"

for i in ${ens}; do # ${NAME} uses the value saved in variable called NAME; it is like the LAMMPS itself.
	for ((j=0; j<${#P[@]}; j++ )); do
		N=`echo ${P[$j]} | awk '{print $1}'`
		l=`echo ${P[$j]} | awk '{print $2}'`
		sig2=`echo ${P[$j]} | awk '{print $3}'`
		n_crowd=`echo ${P[$j]} | awk '{print $4}'`
		epsilon1=`echo ${P[$j]} | awk '{print $5}'`
		dirname=N${N}l${r}sig${sig2}nc${n_crowd}
		[[ -d ${dirname} ]] || mkdir ${dirname}
		cd ${dirname}
		[[ -d archive ]] || mkdir archive
#		[[ -d image.${i} ]] || mkdir image.${i}
		lmp_output="out${i}.txt"
		# Make script file and submit
		# This part should be modified: maybe no sqsub in kisti HPC
		echo "Starting run ensemble ${i} at: $(date)" 
		echo 'module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps-omp/20170811' >> in${i}.sh
		echo '#SBATCH --account=def-byha'>> in${i}.sh
		echo '#SBATCH --ntasks=1' >> in${i}.sh
		echo '#SBATCH --mem-per-cpu=2500M' >> in${i}.sh
		echo '#SBATCH --time=0-00:30' >> in${i}.sh
		echo ${lmp_exe} \
					-var in_filename "../data.chain.${N}" \
					-var run_date ${CURRENT_DATE} \
					-var l ${l} \
					-var sig2 ${sig2} \
					-var n_crowd ${n_crowd} \
					-var epsilon1 ${epsilon1} \
					-var i ${i} \
					-i '../20190229_chain_cube_sharcnet.lmp'\
					"> ${lmp_output}"\
					>> in${i}.sh
		bash "in${i}.sh"
#		jobID=`bash in${i}.sh 2>&1`
#		echo ${jobID}
#		jobID=`echo ${jobID} | cut -d ' ' -f4`
#		echo ${CLUSTER} ${jobID} > in${i}.txt
		# Load the modules
		echo "Ensemble ${i} finished with exit code $? at: $(date)"
		cd ..
	done
done
