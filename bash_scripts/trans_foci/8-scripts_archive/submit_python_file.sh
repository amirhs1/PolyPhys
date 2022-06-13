#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=0-05:30   
#SBATCH --account=rrg-byha
#SBATCH --mail-user=54184877+amirhs1@users.noreply.github.com  
#SBATCH --mail-type=ALL     

# Load python and generate your venv
module load StdEnv/2020 python/3.8
source $HOME/daskEnv/bin/activate

echo "Starting run at: `date`"
python analyze_segments_hists.py
echo "Program glost_launch finished with exit code $? at: `date`"