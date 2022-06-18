while IFS= read -r line; do
    echo "$line"
    cd $line
    sbatch submit.sh
    sleep 5
    cd ..
done < unfinished_run.txt