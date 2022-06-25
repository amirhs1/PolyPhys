 #!/bin/bash

 #SBATCH --account=def-byha
 #SBATCH --ntasks=1
 #SBATCH --mem-per-cpu=2500M      # memory; default unit is megabytes.
 #SBATCH --time=0-00:30           # time (DD-HH:MM).
  
 # Load the module: 
 module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps-omp/20170811
 
 echo "Starting run at: `date`"

 lmp_exec=lmp_icc_openmpi
 lmp_input="modi.lmp"
 lmp_output="output.txt" 
 
 
 ${lmp_exec} -var in_filename data.chain.80  -var run_date 20190303 -var l 36 -var sig2 0.3 -var n_crowd 0 -var epsilon1 5.0 -var i 1 < ${lmp_input} > ${lmp_output}
 
 echo "Program finished with exit code $? at: `date`" 
