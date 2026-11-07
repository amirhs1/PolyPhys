"""
==========================================================
:mod:`polyphys.probe.pairs`
==========================================================

The :mod:`polyphys.probe.pairs` module measures, in every frame of a
trajectory, the distances between the particles of one species, the pairs
in direct contact, and the bonds and clusters those contacts form: the
foci of the *TransFoci* studies, the large monomers of a ring among small
ones. :class:`PairDistances` is an MDAnalysis
:class:`~MDAnalysis.analysis.base.AnalysisBase`, so it runs over a slice of
frames in serial or, split across workers, in parallel, and stores one
array per quantity in its ``results``.

The research code measured the same four quantities in one loop over the
frames of each *TransFoci* study's probe function, without the minimum
image, and wrote the LAMMPS atom IDs of the foci on the diagonal of the
distance matrix. Here the diagonal carries the foci's 0-based backbone
indices (see :class:`PairDistances`, "Notes"), and the minimum image along
the periodic directions is optional.

Contents
========

- :class:`PairDistances`, the per-frame pair analysis.

Dependencies
============
- Python standard library: `math`, `pathlib`, `typing`.
- `NumPy <https://numpy.org/>`_.
- `MDAnalysis <https://docs.mdanalysis.org/>`_ 2.8 or later, through the
  ``mda`` extra: :class:`~MDAnalysis.analysis.base.AnalysisBase` and
  :class:`~MDAnalysis.analysis.results.ResultsGroup`.
- :mod:`polyphys.analyze.clusters`:
  :func:`~polyphys.analyze.clusters.pair_distances`,
  :func:`~polyphys.analyze.clusters.find_direct_contacts`,
  :func:`~polyphys.analyze.clusters.count_bonds`, and
  :func:`~polyphys.analyze.clusters.count_clusters`.
- :mod:`polyphys.spec`:
  :class:`~polyphys.spec.ProbeSpec`.
- :mod:`polyphys.manage.parser`:
  :class:`~polyphys.manage.parser.ParserBase`.

Notation
========
:math:`n`, :math:`F`
    The number of particles of the species and the frames of a run.
:math:`d_{ij}`, :math:`r_c`
    The distance between particles :math:`i < j` in a frame, in the dump's
    length unit, and the contact cutoff: :math:`i` and :math:`j` are in
    direct contact when :math:`d_{ij} \\le r_c`.
:math:`g_i`
    The backbone index of particle :math:`i`, 0-based, from which the
    pipeline's :class:`~polyphys.analyze.pipeline.FociPairDistances` takes
    the genomic distance of each pair.

Examples
========
Three foci at the corners of a 3-4-5 triangle, one frame, with a cutoff of
3.5: the first two are in contact, the third is alone.

>>> import MDAnalysis as mda
>>> import numpy as np
>>> universe = mda.Universe.empty(3, trajectory=True)
>>> universe.atoms.positions = np.array(
...     [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0]])
>>> analysis = PairDistances(universe.atoms, 3.5, [2, 5, 8]).run()
>>> analysis.results['distMatTFoci'][0]
array([[2., 3., 5.],
       [0., 5., 4.],
       [0., 0., 8.]])
>>> analysis.results['bondsHistTFoci'], analysis.results['clustersHistTFoci']
(array([[1, 2, 0]]), array([[0, 1, 1, 0]]))

.. versionadded:: 0.7.0
"""
import math
from pathlib import Path
from typing import Any, Self, Sequence

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe

from ..analyze import clusters
from ..manage.parser import ParserBase
from ..spec import ProbeSpec
from ..manage.types import EntityT, GeometryT, Pathish

# The periodic directions of each geometry, for the minimum image: a
# cylinder is periodic along its axis, z, a cube in every direction.
_PERIODIC: dict[GeometryT, tuple[bool, bool, bool]] = {
    'cubic': (True, True, True),
    'cylindrical': (False, False, True),
}

# The quantities, in the order `results` and `save` list them, as the
# research code named them before the species label.
_QUANTITIES = ('distMatT', 'directContactsMatT', 'bondsHistT',
               'clustersHistT')

_BACKENDS = ('serial', 'multiprocessing', 'dask')


class PairDistances(AnalysisBase):
    """
    The pair distances, direct contacts, bonds, and clusters of one
    species, measured in every frame.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    atomgroup : :class:`~MDAnalysis.core.groups.AtomGroup`
        The particles, at least two, in atom-ID order, which
        ``Universe.select_atoms`` gives; they index the matrices.
    cutoff : float
        The contact cutoff :math:`r_c`, finite and positive, in the dump's
        length unit.
    genomic_positions : sequence of int or numpy.ndarray
        The backbone index :math:`g_i` of each particle, 0-based,
        non-negative, and distinct, one per particle in the group's order;
        written on the diagonal of every frame's distance matrix.
    species : :py:data:`EntityT`, default 'Foci'
        The species label, which names ``results`` and the files
        :meth:`save` writes, ``<name>-<quantity><species>.npy``.
    periodic : tuple of three bool or None, default None
        The directions along which the distances use the minimum image,
        with the frame's box lengths; ``None`` for none, as the research
        code measured.
    inclusive : bool, default True
        Whether a pair at exactly the cutoff is in contact,
        :math:`d_{ij} \\le r_c`, as the research code counted, rather than
        :math:`d_{ij} < r_c`.
    verbose : bool, default False
        Whether to show a progress bar during :meth:`run`.

    Attributes
    ----------
    atomgroup : :class:`~MDAnalysis.core.groups.AtomGroup`
        The particles.
    cutoff : float
        The contact cutoff.
    genomic_positions : numpy.ndarray
        The backbone indices, int64, shape ``(n,)``.
    species : :py:data:`EntityT`
        The species label.
    periodic : tuple of three bool or None
        The minimum-image directions.
    inclusive : bool
        Whether the cutoff itself counts as contact.
    results : :class:`~MDAnalysis.analysis.results.Results`
        After :meth:`run`, with ``S`` the species label: ``distMatTS``,
        float64 ``(n_frames, n, n)``, each frame's :math:`d_{ij}` for
        :math:`i < j` above the diagonal, :math:`g_i` on it, and zeros
        below; ``directContactsMatTS``, bool ``(n_frames, n, n)``, the
        symmetric contact matrix with ``True`` on the diagonal, as
        :func:`~polyphys.analyze.clusters.find_direct_contacts` returns it;
        ``bondsHistTS``, int64 ``(n_frames, n)``, ``[t, d]`` the number of
        particles with ``d`` contacts in frame ``t``; ``clustersHistTS``,
        int64 ``(n_frames, n + 1)``, ``[t, s]`` the number of clusters of
        size ``s`` in frame ``t``, with ``[t, 0]`` always 0.
    n_frames : int
        After :meth:`run`, the number of frames analyzed, from
        :class:`~MDAnalysis.analysis.base.AnalysisBase`.
    frames : numpy.ndarray
        After :meth:`run`, the indices of the frames analyzed, shape
        ``(n_frames,)``.
    times : numpy.ndarray
        After :meth:`run`, the times of the frames analyzed, shape
        ``(n_frames,)``, in the trajectory's time unit.

    Raises
    ------
    ValueError
        If `atomgroup` has fewer than two particles; if `cutoff` is not a
        finite positive number; if `genomic_positions` has not one
        non-negative integer per particle or repeats one; if `species` is
        empty; or if `periodic` is neither ``None`` nor three booleans.

    Notes
    -----
    Per frame, :func:`~polyphys.analyze.clusters.pair_distances` gives the
    full distance matrix, with the minimum image along each direction of
    `periodic` using that frame's box length;
    :func:`~polyphys.analyze.clusters.find_direct_contacts` the contacts
    at the cutoff; :func:`~polyphys.analyze.clusters.count_bonds` the
    bonds histogram; and :func:`~polyphys.analyze.clusters.count_clusters`
    the clusters histogram, from the connected components of the contact
    matrix. The distance matrix is stored as the research code wrote it:
    its upper triangle, with the diagonal replaced.

    **Indexing.** LAMMPS numbers atoms from 1 (``read_data``: "a unique
    value from 1 to Natoms for each atom"), and MDAnalysis keeps those
    numbers as ``AtomGroup.ids`` while numbering its own atoms from 0 as
    ``AtomGroup.indices``. The research code wrote ``foci.atoms.ids`` on
    the diagonal, so the pair tags of its analysis, ``loc<i>loc<j>``,
    carry 1-based IDs. Here the diagonal carries the 0-based position of
    each particle along its chain, which
    :func:`~polyphys.analyze.clusters.genomic_distance` requires and
    :meth:`from_spec` computes; a tag of this package is therefore one
    less than the research tag of the same pair. The genomic distance
    itself is the same under either convention.

    The MDAnalysis contract: :meth:`_prepare` allocates one array per
    quantity for the frames of the run, :meth:`_single_frame` writes the
    current frame's entry, and :meth:`_conclude` does nothing. The class is
    marked parallelizable, with ``serial``, ``multiprocessing``, and
    ``dask`` backends; each worker fills the arrays of its own frames,
    which :meth:`_get_aggregator` concatenates along the frame axis
    (``ResultsGroup.ndarray_vstack``).

    Cost per frame: :math:`O(n^2)` time and space for the distance and
    contact matrices and the connected components. Over a run,
    :math:`8 F n^2` bytes for the distance matrices, :math:`F n^2` for the
    contact matrices, and :math:`8 F (2 n + 1)` for the two histograms:
    for the five foci of the *TransFoci* studies, about 230 bytes per
    frame. Positions come from MDAnalysis as float32, so each distance
    carries a relative round-off of about :math:`10^{-7}`.

    Examples
    --------
    Two particles across the periodic boundary of a box of length 20
    along z, at z = 0.5 and 19.5, are 1 apart with the minimum image and
    19 without:

    >>> import MDAnalysis as mda
    >>> import numpy as np
    >>> universe = mda.Universe.empty(2, trajectory=True)
    >>> universe.load_new(
    ...     np.array([[[0.0, 0.0, 0.5], [0.0, 0.0, 19.5]]]), format='MEMORY',
    ...     dimensions=np.array([[20.0, 20.0, 20.0, 90.0, 90.0, 90.0]]))
    <Universe with 2 atoms>
    >>> wrapped = PairDistances(universe.atoms, 1.5, [0, 9],
    ...                         periodic=(False, False, True)).run()
    >>> wrapped.results['distMatTFoci'][0, 0, 1]
    np.float64(1.0)
    >>> wrapped.results['directContactsMatTFoci'][0]
    array([[ True,  True],
           [ True,  True]])
    >>> PairDistances(universe.atoms, 1.5, [0, 9]).run().results[
    ...     'distMatTFoci'][0, 0, 1]
    np.float64(19.0)
    """

    _analysis_algorithm_is_parallelizable = True

    def __init__(
        self,
        atomgroup: AtomGroup,
        cutoff: float,
        genomic_positions: Sequence[int] | np.ndarray,
        *,
        species: EntityT = 'Foci',
        periodic: tuple[bool, bool, bool] | None = None,
        inclusive: bool = True,
        verbose: bool = False,
    ) -> None:
        super().__init__(atomgroup.universe.trajectory, verbose=verbose)
        n_atoms = atomgroup.n_atoms
        if n_atoms < 2:
            raise ValueError(
                f"'atomgroup' must hold at least two particles, not "
                f"{n_atoms}.")
        if not (isinstance(cutoff, (int, float))
                and not isinstance(cutoff, bool)
                and math.isfinite(cutoff) and cutoff > 0):
            raise ValueError(
                f"'cutoff' must be a finite positive number, not "
                f"{cutoff!r}.")
        positions = np.asarray(genomic_positions)
        if (positions.shape != (n_atoms,)
                or not np.issubdtype(positions.dtype, np.integer)
                or np.any(positions < 0)):
            raise ValueError(
                f"'genomic_positions' must hold one non-negative integer "
                f"per particle, {n_atoms} in all, not "
                f"{list(np.ravel(genomic_positions))!r}.")
        if np.unique(positions).size != n_atoms:
            raise ValueError(
                f"'genomic_positions' repeats an index: "
                f"{positions.tolist()}.")
        if not (isinstance(species, str) and species):
            raise ValueError(
                f"'species' must be a non-empty label, not {species!r}.")
        if periodic is not None and not (
                isinstance(periodic, (tuple, list)) and len(periodic) == 3
                and all(isinstance(flag, bool) for flag in periodic)):
            raise ValueError(
                f"'periodic' must be None or three booleans, not "
                f"{periodic!r}.")
        self.atomgroup = atomgroup
        self.cutoff = float(cutoff)
        self.genomic_positions = positions.astype(np.int64)
        self.species = species
        self.periodic: tuple[bool, bool, bool] | None = (
            None if periodic is None
            else (bool(periodic[0]), bool(periodic[1]), bool(periodic[2])))
        self.inclusive = inclusive
        self._upper = np.triu_indices(n_atoms, 1)
        self._diagonal = np.diag_indices(n_atoms)

    @classmethod
    def from_spec(
        cls,
        universe: Universe,
        spec: ProbeSpec,
        parser: ParserBase,
        *,
        verbose: bool = False,
    ) -> Self:
        """
        Build the analysis of a whole's pair species from its study's probe
        spec.

        The particles are ``spec.pair_species`` selected with its string in
        ``spec.selections``; the cutoff is ``spec.contact_cutoff`` resolved
        for this whole; the backbone indices are the particles' positions
        among the atoms of ``spec.chain``, in atom-ID order; and the
        minimum image follows the geometry when ``spec.minimum_image`` is
        set.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        universe : :class:`~MDAnalysis.core.universe.Universe`
            The trajectory, as :func:`~polyphys.probe.core.make_universe`
            opens it.
        spec : :class:`~polyphys.spec.ProbeSpec`
            The study's probe spec, with ``pair_species`` and
            ``contact_cutoff`` set.
        parser : :class:`~polyphys.manage.parser.ParserBase`
            The whole's parser, for the cutoff.
        verbose : bool, default False
            Whether to show a progress bar during :meth:`run`.

        Returns
        -------
        PairDistances
            The analysis, not yet run.

        Raises
        ------
        ValueError
            If the spec names no pair species, or if a selected particle is
            not among the chain's atoms.

        Examples
        --------
        A ring of six monomers whose atoms 2, 4, and 6 (LAMMPS IDs) are the
        large ones: their backbone indices are 1, 3, and 5.

        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> from polyphys.manage.parser import TransFociCyl
        >>> from polyphys.spec import ProbeSpec, ProjectSpec
        >>> project = ProjectSpec(
        ...     parser=TransFociCyl, group='bug', divisor=0.025,
        ...     physical_attrs=('space', 'nmon_small', 'nmon_large'),
        ...     species_sizes={'Mon': 'dmon_small', 'Foci': 'dmon_large',
        ...                    'Crd': 'dcrowd'},
        ...     tseries_properties=('gyrTMon',))
        >>> spec = ProbeSpec(
        ...     project, {'Mon': 'resid 1', 'Foci': 'type 2'}, 'Mon',
        ...     ('gyrT',), unwrap=False, pair_species='Foci',
        ...     contact_cutoff=lambda info: info.dmon_large + info.dcrowd)
        >>> whole = TransFociCyl(
        ...     'epss5epsl5r5al2nl3ml8ns3ac1nc0lz10dt0.005bdump1adump1ens1',
        ...     'whole', 'bug')
        >>> universe = mda.Universe.empty(
        ...     6, n_residues=1, atom_resindex=[0] * 6, trajectory=True)
        >>> universe.add_TopologyAttr('resid', [1])
        >>> universe.add_TopologyAttr('ids', [1, 2, 3, 4, 5, 6])
        >>> universe.add_TopologyAttr('types', ['1', '2'] * 3)
        >>> universe.load_new(np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        ...     [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 1.0, 0.0],
        ...     [3.0, 2.0, 0.0]]]), format='MEMORY', dimensions=np.array(
        ...     [[10.0, 10.0, 20.0, 90.0, 90.0, 90.0]]))
        <Universe with 6 atoms>
        >>> analysis = PairDistances.from_spec(universe, spec, whole)
        >>> analysis.genomic_positions, analysis.cutoff, analysis.periodic
        (array([1, 3, 5]), 3.0, (False, False, True))
        >>> analysis.run().results['distMatTFoci'][0]
        array([[1.        , 2.        , 2.82842712],
               [0.        , 3.        , 2.        ],
               [0.        , 0.        , 5.        ]])
        """
        if spec.pair_species is None or spec.contact_cutoff is None:
            raise ValueError(
                f"{spec.name}: the probe spec names no pair species.")
        group = universe.select_atoms(spec.selections[spec.pair_species])
        chain = universe.select_atoms(spec.selections[spec.chain])
        in_chain = np.isin(group.ids, chain.ids)
        if not np.all(in_chain):
            raise ValueError(
                f"{spec.name}: the {spec.pair_species} atoms with IDs "
                f"{group.ids[~in_chain].tolist()} are not among the "
                f"{spec.chain} atoms.")
        genomic = np.flatnonzero(np.isin(chain.ids, group.ids))
        periodic = _PERIODIC[spec.geometry] if spec.minimum_image else None
        return cls(group, ProbeSpec.resolve(spec.contact_cutoff, parser),
                   genomic, species=spec.pair_species, periodic=periodic,
                   verbose=verbose)

    @classmethod
    def get_supported_backends(cls) -> tuple[str, ...]:
        """
        The backends :meth:`run` accepts: ``'serial'``,
        ``'multiprocessing'``, and ``'dask'``.

        Examples
        --------
        >>> PairDistances.get_supported_backends()
        ('serial', 'multiprocessing', 'dask')
        """
        return _BACKENDS

    def _key(self, quantity: str) -> str:
        return f"{quantity}{self.species}"

    def _get_aggregator(self) -> ResultsGroup:
        """
        Concatenate each worker's arrays along the frame axis.
        """
        return ResultsGroup(lookup={
            self._key(quantity): ResultsGroup.ndarray_vstack
            for quantity in _QUANTITIES})

    def _prepare(self) -> None:
        """
        Allocate one array per quantity for the frames of this run or
        worker.
        """
        n_atoms = self.atomgroup.n_atoms
        self.results[self._key('distMatT')] = np.zeros(
            (self.n_frames, n_atoms, n_atoms), dtype=np.float64)
        self.results[self._key('directContactsMatT')] = np.zeros(
            (self.n_frames, n_atoms, n_atoms), dtype=np.bool_)
        self.results[self._key('bondsHistT')] = np.zeros(
            (self.n_frames, n_atoms), dtype=np.int64)
        self.results[self._key('clustersHistT')] = np.zeros(
            (self.n_frames, n_atoms + 1), dtype=np.int64)

    def _box(self) -> Any:
        """The box lengths of the periodic directions, 0 elsewhere."""
        if self.periodic is None:
            return None
        lengths = np.asarray(self._ts.dimensions[:3], dtype=np.float64)
        return np.where(self.periodic, lengths, 0.0)

    def _single_frame(self) -> None:
        """
        Measure the four quantities on the current frame and store them in
        this frame's entries.
        """
        distances = clusters.pair_distances(
            self.atomgroup.positions.astype(np.float64), self._box())
        contacts = clusters.find_direct_contacts(
            distances, self.cutoff, inclusive=self.inclusive)
        matrix = np.zeros_like(distances)
        matrix[self._upper] = distances[self._upper]
        matrix[self._diagonal] = self.genomic_positions
        index = self._frame_index
        self.results[self._key('distMatT')][index] = matrix
        self.results[self._key('directContactsMatT')][index] = contacts
        self.results[self._key('bondsHistT')][index] = (
            clusters.count_bonds(contacts))
        self.results[self._key('clustersHistT')][index] = (
            clusters.count_clusters(contacts))

    def save(self, save_to: Pathish, name: str) -> list[Path]:
        """
        Save each quantity's array as ``<name>-<quantity><species>.npy``.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        save_to : str or pathlib.Path
            An existing directory.
        name : str
            The artifact's name, ``<parser.name>-<parser.group>`` in the
            studies, such as ``'epss5...ens1-bug'``.

        Returns
        -------
        list of pathlib.Path
            The files written: the distance matrices, the contact matrices,
            the bonds histograms, and the clusters histograms.

        Raises
        ------
        RuntimeError
            If :meth:`run` has not been called.

        Examples
        --------
        >>> import tempfile
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> universe = mda.Universe.empty(2, trajectory=True)
        >>> universe.atoms.positions = np.array([[0.0, 0.0, 0.0],
        ...                                      [0.0, 0.0, 3.0]])
        >>> analysis = PairDistances(universe.atoms, 1.0, [0, 4]).run()
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     paths = analysis.save(tmp, 'N2ens1-bug')
        ...     len(paths), paths[0].name, np.load(paths[2]).tolist()
        (4, 'N2ens1-bug-distMatTFoci.npy', [[2, 0]])
        """
        if any(self._key(quantity) not in self.results
               for quantity in _QUANTITIES):
            raise RuntimeError("Call run() before save().")
        paths = []
        for quantity in _QUANTITIES:
            path = Path(save_to) / f"{name}-{self._key(quantity)}.npy"
            np.save(path, self.results[self._key(quantity)])
            paths.append(path)
        return paths
