for tar in N*.tar;do
    name=$(echo ${tar} | cut -d . -f -6)
    arr=($(grep -Eo '[[:alpha:]]+|[^a-zA-Z]+' <<< "$name"))
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

    mother=${name}-zip
    #restart_after.i1.gz
    tar -xvf $tar ${mother}/restart_after.i${i}.gz
    mv ${mother}/restart_after.i${i}.gz ./${name}.restart.gz
done