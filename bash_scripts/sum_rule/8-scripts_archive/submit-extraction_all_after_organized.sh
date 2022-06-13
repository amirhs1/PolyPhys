#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-20:00   
#SBATCH --account=rrg-byha
#SBATCH --mail-user=54184877+amirhs1@users.noreply.github.com  
#SBATCH --mail-type=ALL     

python extraction_bug_after_organized.py
echo "Program glost_launch finished with exit code $? at: `date`"