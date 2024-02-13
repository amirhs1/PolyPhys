import numpy as np
from glob import glob
from pathlib import Path

group = 'bug'
incompletes = glob('./N*/N*-'+group+'-*T*.npy')
incompletes = sorted(incompletes)
n_frames = 46000

for f in incompletes:  # Changed 'files' to 'incompletes'
    file_path = Path(f)
    directory = str(file_path.parent)  # Convert Path object to string
    a = np.load(f)
    n_missing = n_frames - a.shape[0]
    idxs = np.arange(0, a.shape[0], dtype=int)
    choices = np.random.choice(idxs[-100:], size=n_missing, replace=False)
    random_elements = a[choices]
    b = np.concatenate((a, random_elements), axis=0)

    # Corrected file path concatenation using pathlib
    original_file_path = file_path.parent / (file_path.stem + '-original.npy')
    new_file_path = file_path.parent / file_path.name

    np.save(original_file_path, a)
    np.save(new_file_path, b)
