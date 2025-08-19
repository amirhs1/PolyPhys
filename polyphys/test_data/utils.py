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

from typing import Optional
from pathlib import Path


def read_n_atoms_from_lammps_trj(trj_path) -> int:
    trj_path = Path(trj_path)
    with trj_path.open() as f:
        for i, line in enumerate(f):
            if i == 3:
                return int(line.strip())
    raise ValueError("File too short to contain number of atoms at line 4.")


def split_lammps_trj(
    filename: str,
    n_wholes: int = 1,
    n_segments: int = 2,
    frames_per_segment: int = 50,
    prefix: Optional[str] = None,
    suffix: Optional[str] = None
) -> None:
    trj_path = Path(filename)
    n_atoms = read_n_atoms_from_lammps_trj(trj_path)
    lines_per_frame = n_atoms + 9
    n_files = n_wholes * n_segments
    n_frames = n_files * frames_per_segment
    n_frames = (n_files - 1) * (frames_per_segment - 1) + frames_per_segment
    total_lines = n_frames * lines_per_frame

    with trj_path.open() as f:
        all_lines = [next(f) for _ in range(total_lines)]

    frames = [
        all_lines[i * lines_per_frame: (i + 1) * lines_per_frame]
        for i in range(n_frames)
    ]
    if prefix is not None:
        if not prefix.endswith('.'):
            prefix += '.'
        else:
            pass
    for i in range(n_wholes):
        for j in range(n_segments):
            file_index = i * n_segments + j
            start = file_index * (frames_per_segment - 1)
            end = start + frames_per_segment
            chunk = frames[start:end]
            out_path = Path(f"{prefix}ens{i+1}.j{j+1}.{suffix}.lammpstrj")
            out_path = trj_path.parent.parent / out_path
            with out_path.open("w") as f:
                for frame in chunk:
                    f.writelines(frame)
            print(f"Wrote {len(chunk)} frames to {out_path}")
