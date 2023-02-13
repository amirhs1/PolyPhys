#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=00-30:00   
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com  
#SBATCH --mail-type=ALL     

# record environment to unclutter gnu parallel run
parallel --record-env

# Load python and generate your venv
module load StdEnv/2020 python/3.9
virtualenv --no-download "${SLURM_TMPDIR}"/env
# trunk-ignore(shellcheck/SC1091)
source "$SLURM_TMPDIR"/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index matplotlib
pip install --no-index DateTime
pip install --no-index numpy
pip install --no-index scipy
pip install --no-index pandas
pip install --no-index seaborn
pip install --no-index sympy
pip install --no-index statsmodels
pip install --no-index MDAnalysis==2.3.0

# Create a function to execute your job
exe(){
dir=${1}
file=$(echo "$dir" | cut -d / -f 1)
parent=$(pwd "${dir}")
cp -vR "${parent}"/"${file}" "${SLURM_TMPDIR}"/
ls "${SLURM_TMPDIR}"/"${file}"/
cd "${SLURM_TMPDIR}"/"${file}"/ || exit
python probe-hns_cub-all_in_one-nucleoid_wholes.py > "${SLURM_TMPDIR}/${file}"-probe-hns_cub-all_in_one-nucleoid_wholes.txt
cp "${SLURM_TMPDIR}"/"${file}"/*.npy "${parent}"/"${file}"/
cp "${SLURM_TMPDIR}"/"${file}"/*.csv "${parent}"/"${file}"/
cp "${SLURM_TMPDIR}"/"${file}"/*.txt "${parent}"/"${file}"/
# trunk-ignore(shellcheck/SC2115)
rm -r "${SLURM_TMPDIR}"/"${file}"
}

echo "Starting run at: $(date)"

# Run GNU Parallel
#export the function
export -f exe

# run the loop in parallel
parallel --will-cite --ungroup  --env _ exe {}-gnuparallel_out-probe-hns_cub-all_in_one-nucleoid_wholes.txt ::: N*/

echo "Program glost_launch finished with exit code $? at: $(date)"