salloc --ntasks=1 --cpus-per-task=4 --mem=8G --time=1:0:0 --account=def-someuser
source daskEnv/bin/activate
jupyter notebook --no-browser --ip $(hostname -f)
ssh -L [local_port]:[hostname]:[host_port] [user]@[cluster] # or
sshuttle --dns -Nr <username>@<cluster>.computecanada.ca
http://localhost:local_port
