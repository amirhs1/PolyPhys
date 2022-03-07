ens="4 5 6 7 8"
for i in $ens;do
    incdir="N2000epsilon5.0r10.5lz336sig1.0nc90720dt0.005bdump1000adump5000ens${i}_incomplete"
    resdir="N2000epsilon5.0r10.5lz336sig1.0nc90720dt0.005bdump1000adump5000ens${i}_res"
    cp ${incdir}/restarts/restart.${i}.66000000 ${resdir}/
done