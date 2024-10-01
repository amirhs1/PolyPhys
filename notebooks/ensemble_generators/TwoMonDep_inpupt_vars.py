from datetime import date
import numpy as np

# Constants
PROJECT_NAME = 'TwoMonDep'
PHI_CS = np.array([0.3, 0.4])                # Crowder volume fractions
RAND_SEED_PARENT = {0.3: 10000, 0.4: 20000}  # Initial seeds for randomness
SIG1 = 5.0                                   # Diameter of monomers
SIG2 = 1.0                                   # Diameter of crowders
N1 = 2                                       # Number of monomers
BOX_SIZE_FACTOR = 5                          # Box size correction factor
DT = 0.0005                                  # Integration time step
B_DUMP_FREQ = 2000                           # Bug dumping frequency
A_DUMP_FREQ = 5000                           # All dumping frequency
T_DUMP_FREQ = 2000                           # Thermo dumping frequency

def generate_ensemble_scripts():
    # Generate date string for script generation info
    today = date.today()
    formatted_date = today.strftime("%Y%m%d")

    bash_tmpl_path = f"./{PROJECT_NAME}_tmpl_ens_gen.sh"

    # Read the template bash file
    with open(bash_tmpl_path) as f:
        bash_template_lines = f.readlines()

    # Compute additional parameters
    wcl_rcut1 = SIG1 * round(2 ** (1 / 6), 6)  # WCA cutoff for monomers
    side_length = BOX_SIZE_FACTOR * SIG1       # Side length of the simulation box
    half_length = side_length / 2              # Half of the simulation box

    # Define surface-to-surface distances
    s_ds_head = np.arange(0, 3 * SIG2, 0.1, dtype=float)
    s_ds_tail = np.arange(3 * SIG2, 5.5 * SIG2, 0.5 * SIG2)
    s_ds = np.round(np.concatenate((s_ds_head, s_ds_tail)), 2)

    # Loop over each crowder volume fraction
    for phi_c in PHI_CS:
        params_list = []
        rand_seeds = np.arange(
            RAND_SEED_PARENT[phi_c],
            RAND_SEED_PARENT[phi_c] + (len(s_ds) + 1) * 10, 10, dtype=int)

        # Loop over each surface-to-surface distance
        for s_d, rand_seed in zip(s_ds, rand_seeds):
            # Calculate number of crowders based on volume fraction
            n2 = int(np.ceil((phi_c * side_length ** 3) / (np.pi * SIG2 ** 3 / 6)))

            # Append formatted parameter string
            params_list.append(f'P[${{#P[@]}}]=\"{SIG1}\t{N1}\t{SIG2}\t{n2}\t{half_length}\t{s_d}\t{rand_seed}\"')

        # Create the script's upper part with all parameter sets
        upper_part = "\n".join(params_list)

        # Combine the upper part with the template and additional variables
        bash_script = f'''#!/bin/bash
# Generated on {formatted_date}
# Parameters: sig1, n1, sig2, n2, hl, s_d, rand_seed
{upper_part}

dt={DT}
adump={A_DUMP_FREQ}
bdump={B_DUMP_FREQ}
tdump={T_DUMP_FREQ}

{''.join(bash_template_lines)}
'''

        # Generate a filename and write the script to a file
        filename = f"ens_gen-{PROJECT_NAME}-am{SIG1}nm{N1}ac{SIG2}nc{n2}.sh"
        with open(filename, 'w') as f:
            f.write(bash_script)
        print(f"Bash script generated: '{filename}'")


if __name__ == "__main__":
    generate_ensemble_scripts()
