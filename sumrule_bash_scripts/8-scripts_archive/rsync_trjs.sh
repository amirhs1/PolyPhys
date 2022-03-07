for dir in N2*ac1*-trjs_all/; do
    name=${dir::-1}
    echo ${dir::-1}
    nohup ls $dir > ${dir::-1}-new.txt &
    echo nohup rsync -axvH --no-g --no-p  $HOME/amirhsi_projects/cylinder_simulations/$name $HOME/scratch/ &
done
