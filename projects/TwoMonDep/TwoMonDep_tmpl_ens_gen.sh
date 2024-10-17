# Initialize ensemble values
ens='1 2 3 4' # Global variable defining total ensembles for simulation.

# Loop through ensemble values
for i in ${ens}; do 
    # Loop through parameter sets
    for ((j=0; j < ${#P[@]}; j++ )); do
        # Extract parameters from array P
        read -r sig1 n1 sig2 n2 hl sd randseed_group <<< "${P[$j]}"

        # Define directory and simulation names
        dirname="am${sig1}nm${n1}ac${sig2}nc${n2}hl${hl}sd${sd}dt${dt}bdump${bdump}adump${adump}tdump${adump}ens${i}"
        simname="am${sig1}nm${n1}ac${sig2}nc${n2}hl${hl}sd${sd}dt${dt}bdump${bdump}adump${adump}tdump${adump}ens${i}"

        # Create simulation input file from template
        cp TwoMonDep_tmpl_input.lmp input.lmp
        {
            echo "# ---------- Input Parameters ----------"
            echo "variable        i equal ${i} # Simulation counter"
            echo "variable        sig1 equal ${sig1} # Diameter of monomers"
            echo "variable        n1 equal ${n1} # Number of monomers"
            echo "variable        sig2 equal ${sig2} # Diameter of crowders"
            echo "variable        n2 equal ${n2} # Number of crowders"
            echo "variable        sd equal ${sd} # Surface-to-surface distance between monomers"
            echo "variable        hl equal ${hl} # Box size half-length"
            echo "variable        randseed equal $((i+randseed_group)) # Random seed"
            echo "variable        simname string ${simname} # Simulation/run name"
        } | cat - input.lmp > temp && mv temp input.lmp

        # Create directory if it does not exist and move input file
        mkdir -p "${dirname}"
        mv input.lmp "${dirname}"
        
        # Select appropriate submission script based on total number of particles
        total_particles=$((n2 + n1))
        if [ "${total_particles}" -ne "${n1}" ]; then
            if [ "${total_particles}" -le 5000 ]; then
                cp submit_cpu4_le5000.sh submit.sh
            elif [ "${total_particles}" -le 15000 ]; then
                cp submit_cpu4_gt5000_le15000.sh submit.sh
            elif [ "${total_particles}" -le 30000 ]; then
                cp submit_cpu8_gt15000_le30000.sh submit.sh
            elif [ "${total_particles}" -le 75000 ]; then
                cp submit_cpu16_gt30000_le75000.sh submit.sh
            elif [ "${total_particles}" -le 150000 ]; then
                cp submit_cpu32_gt75000_le150000.sh submit.sh
            elif [ "${total_particles}" -le 200000 ]; then
                cp submit_cpu32_gt150000_le200000.sh submit.sh
            elif [ "${total_particles}" -le 250000 ]; then
                cp submit_cpu32_gt200000_le250000.sh submit.sh
            else
                cp submit_cpu32_gt250000_le300000.sh submit.sh
            fi
        else
            cp submit_nocrowd_8hrs.sh submit.sh
        fi
        
        # Move the submission script to the directory
        mv submit.sh "${dirname}"
        
        # Enter the directory and create necessary sub-directory
        cd "${dirname}" || exit
        mkdir -p restarts
        cd ..
    done
done