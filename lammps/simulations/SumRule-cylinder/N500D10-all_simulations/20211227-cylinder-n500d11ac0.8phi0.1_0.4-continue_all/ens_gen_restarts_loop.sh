firstseed=100206 # first random seed for this simulation group
TEMPFILE="./rand.txt"
echo $firstseed > $TEMPFILE
cat $TEMPFILE
for file in N*.restart; do
    echo $file
    dirname=$(echo ${file} | cut -d . -f -6)
    mkdir $dirname
    randseed=$[$(cat $TEMPFILE)] #current randseed
    arr=($(grep -Eo '[[:alpha:]]+|[^a-zA-Z]+' <<< "$file"))
    N=${arr[1]}
    epsilon1=${arr[3]}
    r=${arr[5]}
    lz=${arr[7]}
    sig2=${arr[9]}
    n_crowd=${arr[11]}
    run_dt=${arr[13]}
    bug_dump=${arr[15]}
    all_dump=${arr[17]}
    ens=${arr[19]}
    i=${ens:0:1}
    cp restart_input_cylinder_sumrule_loop.lmp input.lmp
    echo "variable randseed equal $randseed" | cat - input.lmp > temp && mv temp input.lmp
    echo "variable n_bug equal $N" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable r equal $r" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable i equal $i" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable lz equal ${lz}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable n_crowd equal ${n_crowd}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable sig2 equal ${sig2}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable epsilon1 equal ${epsilon1}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable run_dt equal ${run_dt}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable bug_dump equal ${bug_dump}" | cat - input.lmp > temp && mv temp input.lmp
	echo "variable all_dump equal ${all_dump}" | cat - input.lmp > temp && mv temp input.lmp
	echo "read_restart  ${file}" | cat - input.lmp > temp && mv temp input.lmp
	echo '#Defining input parameters:' | cat - input.lmp > temp && mv temp input.lmp

    mv input.lmp ${dirname}
    if [ $((${n_crowd}+${N})) -ne ${N} ]; then
		if [ $((${n_crowd}+${N})) -le 5000 ]; then
			cp submit_cpu2_le5000.sh submit.sh && mv submit.sh ${dirname}
		elif [ $((${n_crowd}+${N})) -gt 5000 ] && [ $((${n_crowd}+${N})) -le 10000 ]; then
			cp submit_cpu4_gt5000_le10000.sh submit.sh && mv submit.sh ${dirname}
		elif [ $((${n_crowd}+${N})) -gt 10000 ] && [ $((${n_crowd}+${N})) -le 25000 ]; then
			cp submit_cpu8_gt10000_le25000.sh submit.sh && mv submit.sh ${dirname}
		elif [ $((${n_crowd}+${N})) -gt 25000 ] && [ $((${n_crowd}+${N})) -le 50000 ]; then
			cp submit_cpu8_gt25000_le50000.sh submit.sh && mv submit.sh ${dirname}
		elif [ $((${n_crowd}+${N})) -gt 50000 ] && [ $((${n_crowd}+${N})) -le 100000 ]; then
			cp submit_cpu16_gt50000_le100000.sh submit.sh && mv submit.sh ${dirname}
		else
			cp submit_cpu32_gt100000_le200000.sh submit.sh && mv submit.sh ${dirname}
		fi
	else
		cp submit_nocrowd.sh submit.sh && mv submit.sh ${dirname}
	fi
    
    cp $file $dirname
    cd ${dirname}
	[[ -d restarts ]] || mkdir restarts
	cd ..

    COUNTER=$[$(cat $TEMPFILE) + 1]
    echo $COUNTER
    echo $COUNTER > $TEMPFILE
done

