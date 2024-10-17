#!/bin/bash
echo "Start to submit"
for dir in am*ens[1-4]/; do # TransFociCub
   cd $dir
   sbatch submit.sh
   sleep 5
   cd ..
done
echo "Finished!"
