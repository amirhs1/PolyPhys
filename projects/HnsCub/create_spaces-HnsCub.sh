hns='0 4 8 12 16 20 30'
epshc='0 1'
ac = '2'
for n in $hns; do
    for a in $ac: do
        for e in $epshc;do
            # N#kbmm#nh#ac#epshc#
            simdir=N200kbmm2.0nh${n}ac${a}epshc${e}-all_simulations
            mkdir ${simdir}
            # N${N}kbmm${kBend}nh${n_hns}ac${sig4}l${l}epshc${eps_hc}nc${n_crowd}ens${i}.ring
            mv N200kbmm2nh${n}ac${a}*epshc${e}*ring  ${simdir}/
            cp -R run_files run_files-N200kbmm2.0nh${n}ac${a}epshc${e}
            mv run_files-N200kbmm2.0nh${n}ac${a}epshc${e}} ${simdir}/
            cp check*.sh loop_*.sh trj_ex*.sh sim*.sh ${simdir}/
        done
    done
done