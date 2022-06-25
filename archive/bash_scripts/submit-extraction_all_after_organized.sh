#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-20:00   
#SBATCH --account=rrg-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com  
#SBATCH --mail-type=ALL     

python extraction_bug_after_organized.py
echo "Program glost_launch finished with exit code $? at: `date`"