for dir in N*-restarts/; do
    name=${dir::-1}
    echo $name
    nohup rsync -axvH --no-g --no-p  $HOME/amirhsi_rrg/cylinder_simulations/$name $HOME/scratch/ > ${name}.log &
done