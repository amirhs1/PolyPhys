hns='0 4 8 12 16 20 30'
ac='1 2'
D=8.0
for n in $hns; do
    for a in $ac;do
        simdir=N200D${D}nh${n}ac${a}-all_simulations
        mkdir ${simdir}
        mv N200kbmm2r4.5nh${n}ac${a}*ring  ${simdir}/
        cp -R run_files run_files-N200D${D}nh${n}ac${a}
        mv run_files-N200D${D}nh${n}ac${a} ${simdir}/
        cp check*.sh loop_*.sh trj_ex*.sh sim*.sh ${simdir}/
    done
done