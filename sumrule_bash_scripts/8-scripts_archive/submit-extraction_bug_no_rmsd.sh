#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-07:30   
#SBATCH --account=rrg-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com  
#SBATCH --mail-type=ALL     

# record environment to unclutter gnu parallel run
parallel --record-env

# Load python and generate your venv
module load StdEnv/2020 python/3.8
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index matplotlib
pip install --no-index DateTime
pip install --no-index numpy
pip install --no-index scipy
pip install --no-index pandas
pip install --no-index seaborn
pip install --no-index sympy
pip install --no-index MDAnalysis

# Create a function to execute your job
exe(){
dir=${1}
file=$(echo $dir | cut -d / -f 1)
(cd ${file} && python extraction_bug_no_rmsd.py > ${file}-gnuarallel_out-bug_no_rmsd.txt)
#echo (cd ${file} && python analyze_run.py) this for testing.
}


echo "Starting run at: `date`"

# Run GNU Parallel
#export the function
export -f exe

# run the loop in parallel
parallel --will-cite --env _ exe {}-gnuparallel-bug_no_rmsd.txt ::: N*/

echo "Program glost_launch finished with exit code $? at: `date`"
