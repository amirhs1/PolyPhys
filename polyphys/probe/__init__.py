"""The extract submodule contains python scripts that directly calculate \
physical properties from the LAMMPS trjaectory file, using MDAnalysis \
package."""
from polyphys.probe import prober

__all__ = ["prober"]
