for tar in N*.tar;do
    mother=${tar::-4}-zip
    tar -xvf $tar ${mother}/log.lammps.gz
done