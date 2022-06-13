# this only works for padding with 0 upt to 10:
for j in $(seq -f "%02g" 2 10);do
    jold=$((10#$j-1)) # 10# means base 10
    in="N500epsilon5.0r5.5lz205.5sig0.8nc12012dt0.002bdump1000adump5000ens2.j0${jold}.all.lammpstrj"
    out="N500epsilon5.0r5.5lz205.5sig0.8nc12012dt0.002bdump1000adump5000ens2.j${j}.all.lammpstrj"
    tail -n 12521 $in | cat - $out > temp && mv temp $out
done
