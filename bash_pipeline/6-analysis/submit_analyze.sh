#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-06:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL



# Load python and generate your venv
module load StdEnv/2020 python/3.9
virtualenv --no-download "${SLURM_TMPDIR}"/env
# trunk-ignore(shellcheck/SC1091)
source "${SLURM_TMPDIR}"/env/bin/activate
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

# python file
echo "Starting run at: $(date)"
python analyze.py > analyze.txt
echo "Program finished with exit code $? at: $(date)"
