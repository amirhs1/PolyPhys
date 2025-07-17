"""\
generate_test_data
==================

Utility functions for generating and preparing LAMMPS test data files.

This module contains tools to transform large raw simulation data files
(e.g., LAMMPS trajectory files) into smaller, compressed test data files
used internally in the PolyPhys test suite. The original source data files
are not part of the repository, but the functions here document and 
automate the generation process for transparency and reproducibility.

Functions
---------
- read_n_atoms(trj_path)
    Read number of atoms from line 4 of a LAMMPS trajectory file.

- split_trj_gz(trj_path, n_frames, n_files, output_prefix)
    Extract the last N frames from a LAMMPS trajectory file and write them
    as compressed `.trj.gz` test files, optionally overlapping frames between
    files.

Notes
-----
This module is intended for development use only and should not be imported
as part of the core PolyPhys package.
"""

from pathlib import Path


def read_n_atoms(trj_path: Path) -> int:
    with trj_path.open() as f:
        for i, line in enumerate(f):
            if i == 3:
                return int(line.strip())
    raise ValueError("File too short to contain number of atoms at line 4.")


def split_trj(trj_path: Path, n_frames: int = 199, n_files: int = 2, output_prefix="split"):
    n_atoms = read_n_atoms(trj_path)
    lines_per_frame = n_atoms + 9

    with trj_path.open() as f:
        all_lines = f.readlines()

    total_lines = len(all_lines)
    total_frames = total_lines // lines_per_frame
    if total_frames < n_frames:
        raise ValueError(f"File only contains {total_frames} frames, less than {n_frames}.")

    # Extract only the last `n_frames`
    start_line = (total_frames - n_frames) * lines_per_frame
    recent_lines = all_lines[start_line:]

    frames = [
        recent_lines[i * lines_per_frame: (i + 1) * lines_per_frame]
        for i in range(n_frames)
    ]

    # Determine split points with 1 overlapping frame
    chunk_size = (n_frames + 1) // n_files
    for i in range(n_files):
        start = i * (chunk_size - 1)
        end = start + chunk_size
        chunk = frames[start:end]
        out_path = trj_path.parent / f"{output_prefix}_{i+1}.trj"
        with out_path.open("w") as f:
            for frame in chunk:
                f.writelines(frame)
        print(f"Wrote {len(chunk)} frames to {out_path}")

