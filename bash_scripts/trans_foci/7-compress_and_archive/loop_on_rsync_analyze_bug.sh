for file in N*-analyze_bug.tar;do
    rsync -axvH --no-g --no-p  $HOME/scratch/analyze_bugs/$file $HOME/user_projects/cylinder_simulations/
done