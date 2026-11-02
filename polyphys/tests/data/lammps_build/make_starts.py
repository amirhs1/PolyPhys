"""Write the start configurations that ``in.kg`` checks.

Run from the project environment::

    python make_starts.py OUTDIR

It writes, with :mod:`polyphys.simulate.build`, ``helix.data`` (a linear
helix of 200 beads, radius 2, pitch 2.5, atom style ``bond``),
``ring.data`` (a plectoneme ring of 200 beads, the same radius and pitch,
atom style ``bond``), and ``ring_angle.data`` (the same ring with angles,
atom style ``angle``), each in a box that is periodic along z only and wide
enough for a cylinder of diameter 5.2 about the z axis, and prints the
:func:`~polyphys.simulate.build.check_conformation` report of each. The
same configurations are built by ``test_build.py``; this script is for a
run by hand with ``check_starts.sh``.
"""

import sys
from pathlib import Path

import numpy as np

from polyphys.simulate.build import (
    chain_data,
    check_conformation,
    helix,
    plectoneme_ring,
)
from polyphys.simulate.core import Box
from polyphys.simulate.data import write_data

N_BEADS = 200
RADIUS = 2.0
PITCH = 2.5
HALF_WIDTH = 5.0        # x and y half-width of the box; the wall at 3.1 fits
MARGIN = 2.0            # z margin on each side, beyond the ghost cutoff
GHOST_CUTOFF = 1.122462 + 0.3
PERIODIC = (False, False, True)


def start_box(positions: np.ndarray) -> Box:
    """Return a box of the configuration's height plus the margins."""
    half_height = 0.5 * float(np.ptp(positions[:, 2])) + MARGIN
    return Box(lo=np.array([-HALF_WIDTH, -HALF_WIDTH, -half_height]),
               hi=np.array([HALF_WIDTH, HALF_WIDTH, half_height]))


def main(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    line = helix(N_BEADS, radius=RADIUS, pitch=PITCH)
    ring = plectoneme_ring(N_BEADS, radius=RADIUS, pitch=PITCH)
    files = {
        "helix.data": chain_data(line, start_box(line)),
        "ring.data": chain_data(ring, start_box(ring), ring=True),
        "ring_angle.data": chain_data(ring, start_box(ring), ring=True,
                                      angles=True),
    }
    for name, data in files.items():
        cylinder = {"dcyl": 5.2, "dmon": 1.0} if name == "helix.data" else {}
        report = check_conformation(data, r0=1.5, ghost_cutoff=GHOST_CUTOFF,
                                    periodic=PERIODIC, **cylinder)
        path = write_data(data, out_dir / name,
                          title=f"{name}: polyphys.simulate.build start")
        print(f"{path.name}: {report}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python make_starts.py OUTDIR")
    main(Path(sys.argv[1]))
