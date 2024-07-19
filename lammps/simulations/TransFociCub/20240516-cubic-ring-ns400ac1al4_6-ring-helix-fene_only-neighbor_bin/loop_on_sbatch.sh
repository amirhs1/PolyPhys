#!/bin/bash
echo "Start to submit"
for directory in al*/; do # TransFociCub
#for directory in N*[1-4]/; do # TransFociCub
#for directory in eps*/; do # TransFociCyl
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
echo "Finished!"
