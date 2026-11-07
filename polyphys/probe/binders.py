"""
==========================================================
:mod:`polyphys.probe.binders`
==========================================================

The :mod:`polyphys.probe.binders` module counts, in every frame of a
trajectory, how the two-patch binders of a study bind its polymer chain:
how many patches touch a monomer, how many binders are free, dangling by
one patch, or bridging two monomers, how many bridged pairs are cis or
trans, and the histogram of the bridged pairs' genomic distances, the loop
lengths. :class:`BinderStatistics` is an MDAnalysis
:class:`~MDAnalysis.analysis.base.AnalysisBase`, so it runs over a slice of
frames in serial or, split across workers, in parallel, and stores one
array per count in its ``results``.

The research code measured the same eight counts in the H-NS probe of each
study: the distances between every monomer and every patch, with no
minimum image, the contacts below a cutoff, and the counts of its
``hns_binding``. Here the distances take the minimum image along the
directions the study's geometry makes periodic, the cutoff and the
comparison are the study's :class:`~polyphys.spec.BinderSpec`,
and the counting is :func:`~polyphys.analyze.clusters.binding_statistics`.

Contents
========

- :class:`BinderStatistics`, the per-frame binding analysis.

Dependencies
============
- Python standard library: `math`, `pathlib`, `typing`.
- `NumPy <https://numpy.org/>`_.
- `MDAnalysis <https://docs.mdanalysis.org/>`_ 2.8 or later, through the
  ``mda`` extra: :class:`~MDAnalysis.analysis.base.AnalysisBase` and
  :class:`~MDAnalysis.analysis.results.ResultsGroup`.
- :mod:`polyphys.analyze.clusters`:
  :func:`~polyphys.analyze.clusters.find_direct_contacts`,
  :func:`~polyphys.analyze.clusters.binding_statistics`, and
  :class:`~polyphys.analyze.clusters.BindingStats`.
- :mod:`polyphys.spec`:
  :data:`~polyphys.spec.BINDER_SERIES` and
  :class:`~polyphys.spec.ProbeSpec`.

Notation
========
:math:`N_m`, :math:`N_p`, :math:`N_b = N_p / 2`
    The numbers of monomers, patches, and binders; the two patches of
    binder :math:`k` are the patches :math:`2k` and :math:`2k + 1` in
    atom-ID order.
:math:`r_{ij}`
    The distance between monomer :math:`i` and patch :math:`j`, in the
    dump's length unit, the monomer diameter in the studies; along a
    periodic direction of box length :math:`L`, the component is reduced
    to :math:`[-L/2, L/2)` first (the minimum image).
:math:`r_\\mathrm{thres}`
    The binding threshold, `cutoff`: monomer :math:`i` and patch :math:`j`
    are in contact when :math:`r_{ij} \\le r_\\mathrm{thres}` (`inclusive`)
    or :math:`r_{ij} < r_\\mathrm{thres}`.
:math:`d`
    The genomic distance of a bridged pair, in bonds along the backbone,
    the shorter way round on a ring; a bridge is cis when
    :math:`d \\le` `cis_threshold` and trans otherwise.

Examples
========
A linear chain of four beads one unit apart and one binder whose two
patches sit 0.3 above beads 0 and 2, in a one-frame universe:

>>> import MDAnalysis as mda
>>> import numpy as np
>>> universe = mda.Universe.empty(6, trajectory=True)
>>> universe.atoms.positions = np.array(
...     [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0],
...      [0.0, 0.3, 0.0], [2.0, 0.3, 0.0]])
>>> analysis = BinderStatistics(
...     universe.atoms[:4], universe.atoms[4:], cutoff=0.5,
...     topology='linear', cis_threshold=1).run()
>>> analysis.results['n_binders_bridging'], analysis.results['n_pairs_trans']
(array([1]), array([1]))
>>> analysis.results['loop_length_histogram']
array([[0, 0, 1, 0]])

.. versionadded:: 0.7.0
"""
import math
from pathlib import Path
from typing import Mapping, Self, cast, get_args

import numpy as np
import numpy.typing as npt
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe

from ..analyze.clusters import (
    BindingStats,
    Topology,
    binding_statistics,
    find_direct_contacts,
)
from ..spec import BINDER_SERIES, ProbeSpec
from ..manage.types import EntityT, GeometryT, Pathish, TopologyT

#: The result key of the per-frame loop-length histogram.
LOOP_KEY = 'loop_length_histogram'

#: The result key of the per-frame monomer-to-patch distance matrices, when
#: they are stored.
DISTANCES_KEY = 'distances'

# The directions that each geometry makes periodic: every one in a periodic
# cube, the axis only in a cylinder (LAMMPS boundary ``ff ff pp``).
_PERIODIC: Mapping[GeometryT, tuple[bool, bool, bool]] = {
    'cubic': (True, True, True),
    'cylindrical': (False, False, True),
}

_BACKENDS = ('serial', 'multiprocessing', 'dask')

# The number of patches per binder that binding_statistics counts.
_TWO_PATCHES = 2


def _pair_distances(
    reference: npt.NDArray[np.float64],
    configuration: npt.NDArray[np.float64],
    box_lengths: npt.NDArray[np.float64],
    periodic: tuple[bool, bool, bool],
) -> npt.NDArray[np.float64]:
    """
    The distance between every row of `reference` and every row of
    `configuration`, shape ``(n_reference, n_configuration)``, with the
    minimum image along each periodic direction of `box_lengths`.
    """
    diff = reference[:, None, :] - configuration[None, :, :]
    for axis in np.flatnonzero(periodic):
        length = box_lengths[axis]
        diff[..., axis] -= length * np.round(diff[..., axis] / length)
    return np.asarray(np.sqrt(np.einsum('ijk,ijk->ij', diff, diff)),
                      dtype=np.float64)


class BinderStatistics(AnalysisBase):
    """
    The binding of two-patch binders to one polymer chain, counted in every
    frame.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    chain : :class:`~MDAnalysis.core.groups.AtomGroup`
        The chain, at least two monomers, in atom-ID order: the order along
        the backbone, which ``Universe.select_atoms`` gives.
    patches : :class:`~MDAnalysis.core.groups.AtomGroup`
        The binders' patches, in atom-ID order, from the same universe as
        `chain`: the two patches of binder :math:`k` at positions
        :math:`2k` and :math:`2k + 1`, so their number is a positive
        multiple of `patches_per_binder`. The binder cores are not needed.
    cutoff : float
        The binding threshold :math:`r_\\mathrm{thres}`, in the dump's
        length unit; finite and positive.
    topology : :py:data:`TopologyT`
        ``'linear'`` or ``'ring'``, whether the chain is open or closed,
        which sets the genomic distance of a bridged pair and the length of
        the loop-length histogram.
    inclusive : bool, default True
        Whether a pair at exactly `cutoff` is a contact
        (:math:`r \\le r_\\mathrm{thres}`) or not
        (:math:`r < r_\\mathrm{thres}`).
    cis_threshold : int, default 4
        The largest genomic distance, in bonds, of a cis bridge;
        non-negative.
    patches_per_binder : int, default 2
        The number of patches per binder; only ``2`` is accepted, which
        :func:`~polyphys.analyze.clusters.binding_statistics` counts.
    periodic : tuple of three bool, default (False, False, False)
        Whether the minimum image is taken along ``x``, ``y``, and ``z``,
        with the box lengths of each frame. The research probe took none;
        :meth:`from_spec` sets every direction in a cube and ``z`` only in
        a cylinder.
    store_distances : bool, default False
        Whether the monomer-to-patch distance matrix of every frame is
        kept in ``results`` and written by :meth:`save`.
    chain_label, patch_label, binder_label : :py:data:`EntityT`
        The species labels that name the files :meth:`save` writes;
        ``'Mon'``, ``'HnsPatch'``, and ``'HnsCore'`` by default, the H-NS
        studies' names.
    verbose : bool, default False
        Whether to show a progress bar during :meth:`run`.

    Attributes
    ----------
    chain, patches : :class:`~MDAnalysis.core.groups.AtomGroup`
        The chain and the patches.
    cutoff : float
        The binding threshold.
    topology : :py:data:`TopologyT`
        The chain's topology.
    inclusive : bool
        Whether a pair at exactly `cutoff` is a contact.
    cis_threshold : int
        The largest genomic distance of a cis bridge, in bonds.
    patches_per_binder : int
        The number of patches per binder, ``2``.
    periodic : tuple of three bool
        The directions along which the minimum image is taken.
    store_distances : bool
        Whether the distance matrices are kept.
    chain_label, patch_label, binder_label : :py:data:`EntityT`
        The species labels of the saved files.
    n_lengths : int
        The length of a loop-length histogram row: ``n_monomers // 2 + 1``
        on a ring, ``n_monomers`` on a linear chain, as
        :class:`~polyphys.analyze.clusters.BindingStats` sizes it.
    results : :class:`~MDAnalysis.analysis.results.Results`
        After :meth:`run`: one int64 array of shape ``(n_frames,)`` per
        attribute of :class:`~polyphys.analyze.clusters.FrameBinding` in
        :data:`~polyphys.spec.BINDER_SERIES`, keyed by that
        name; ``'loop_length_histogram'``, int64 of shape
        ``(n_frames, n_lengths)``, whose row ``[t, d]`` is the number of
        distinct bridged pairs at :math:`d` bonds in frame ``t``, with
        column 0 always 0; and, with `store_distances`, ``'distances'``,
        float64 of shape ``(n_frames, n_monomers, n_patches)``.
    n_frames : int
        After :meth:`run`, the number of frames analyzed, from
        :class:`~MDAnalysis.analysis.base.AnalysisBase`.
    frames, times : numpy.ndarray
        After :meth:`run`, the indices and times of the frames analyzed,
        shape ``(n_frames,)``.

    Raises
    ------
    ValueError
        If `chain` has fewer than two atoms; if `patches` is empty, is not
        from the universe of `chain`, or has a number of atoms that is not
        a multiple of `patches_per_binder`; if `cutoff` is not a finite
        positive number; if `topology` is not ``'linear'`` or ``'ring'``;
        if `cis_threshold` is negative; if `patches_per_binder` is not
        ``2``; if `periodic` is not three booleans; or if a label is empty.
        During :meth:`run`, if a periodic direction is asked for and the
        frame carries no box.

    Notes
    -----
    Per frame, the analysis forms :math:`r_{ij}` for every monomer and
    patch, with the minimum image along the periodic directions, marks the
    contacts with :func:`~polyphys.analyze.clusters.find_direct_contacts`,
    and counts them with
    :func:`~polyphys.analyze.clusters.binding_statistics`, whose rules hold
    here: the patch and binder counts use the contacts as marked, so a
    patch within the cutoff of two adjacent monomers is one engaged patch;
    the bridged pairs use one monomer per patch, the lowest-index one, so
    each bridging binder joins exactly two monomers; and cis and trans
    count distinct monomer pairs. The histogram row is the bincount of the
    bridged pairs' genomic distances.

    The binder model is the bivalent patchy protein of Brackley et al.
    [1]_, with a patch at each pole of a repulsive core, that the H-NS
    studies follow [2]_; the dangling, cis, and trans modes and the cis
    threshold of two intervening monomers are those of Joyeux and Vreede
    [3]_, whose cis bridge joins sites "separated by at most two DNA
    beads", a genomic distance of at most 3 bonds, which `cis_threshold`
    of 4 reproduces on the studies' chains. A patch that sits between two
    consecutive monomers can touch both, which Brackley [4]_ notes as the
    one way a single patch reaches two beads; the single-patch rule
    resolves it.

    The MDAnalysis contract: :meth:`_prepare` allocates the arrays for the
    frames of the run, :meth:`_single_frame` writes the current frame's
    row, and :meth:`_conclude` does nothing. The class is marked
    parallelizable, with ``serial``, ``multiprocessing``, and ``dask``
    backends; each worker fills the arrays of its own frames, which
    :meth:`_get_aggregator` concatenates along the frame axis
    (``ResultsGroup.ndarray_hstack`` for the eight series,
    ``ResultsGroup.ndarray_vstack`` for the histogram and the distances).

    Cost per frame: :math:`O(N_m N_p)` time for the distances and the
    contacts, plus :math:`O(N_b \\log N_b)` to make the bridged pairs
    distinct, and :math:`O(N_m N_p)` space for the distance and contact
    matrices. Over a run of :math:`F` frames, :math:`O(F)` space per
    series, :math:`O(F N_m / 2)` for the histogram, and, with
    `store_distances`, :math:`8 F N_m N_p` bytes for the float64
    distances, held in memory until :meth:`save`, as the research code
    held its list of matrices. Positions come from MDAnalysis as float32,
    so each distance carries a relative round-off of about :math:`10^{-7}`
    before it is compared with the cutoff.

    References
    ----------
    .. [1] C. A. Brackley, S. Taylor, A. Papantonis, P. R. Cook, and D.
       Marenduzzo, "Nonspecific bridging-induced attraction drives
       clustering of DNA-binding proteins and genome organization", Proc.
       Natl. Acad. Sci. U.S.A. 110, E3605-E3611 (2013).
       https://doi.org/10.1073/pnas.1302950110
    .. [2] Y. Jung, A. Sadeghi, and B.-Y. Ha, "Modeling the compaction of
       bacterial chromosomes by biomolecular crowding and the cross-linking
       protein H-NS", Sci. Rep. 14, 139 (2024).
       https://doi.org/10.1038/s41598-023-50355-2
    .. [3] M. Joyeux and J. Vreede, "A Model of H-NS Mediated Compaction of
       Bacterial DNA", Biophys. J. 104, 1615-1622 (2013).
       https://doi.org/10.1016/j.bpj.2013.02.043
    .. [4] C. A. Brackley, "Polymer compaction and bridging-induced
       clustering of protein-inspired patchy particles", J. Phys.: Condens.
       Matter 32, 314002 (2020). https://doi.org/10.1088/1361-648X/ab7f6c

    Examples
    --------
    A ring of six beads and two binders over one frame: binder 0 bridges
    beads 0 and 2, two bonds apart; binder 1 dangles from bead 4 with its
    other patch far away.

    >>> import MDAnalysis as mda
    >>> import numpy as np
    >>> chain = np.column_stack([np.arange(6.0), np.zeros(6), np.zeros(6)])
    >>> patches = np.array([[0.0, 0.3, 0.0], [2.0, 0.3, 0.0],
    ...                     [4.0, 0.3, 0.0], [9.0, 9.0, 9.0]])
    >>> universe = mda.Universe.empty(10, trajectory=True)
    >>> universe.atoms.positions = np.vstack([chain, patches])
    >>> analysis = BinderStatistics(
    ...     universe.atoms[:6], universe.atoms[6:], cutoff=0.5,
    ...     topology='ring', cis_threshold=4).run()
    >>> [int(analysis.results[key][0]) for key in BINDER_SERIES]
    [3, 1, 3, 0, 1, 1, 1, 0]
    >>> analysis.results['loop_length_histogram']
    array([[0, 0, 1, 0]])

    The eight values are, in the order of
    :data:`~polyphys.spec.BINDER_SERIES`: three patch-monomer
    contacts, one free and three engaged patches, no free binder, one
    bridging and one dangling binder, one cis pair and no trans pair.

    With the minimum image along ``x`` in a box of length 10, a patch at
    ``x = 9.8`` touches bead 0:

    >>> universe.trajectory.ts.dimensions = [10.0, 10.0, 10.0, 90, 90, 90]
    >>> universe.atoms[9].position = [9.8, 0.0, 0.0]
    >>> wrapped = BinderStatistics(
    ...     universe.atoms[:6], universe.atoms[6:], cutoff=0.5,
    ...     topology='ring', periodic=(True, False, False)).run()
    >>> int(wrapped.results['n_binders_bridging'][0])
    2
    """

    _analysis_algorithm_is_parallelizable = True

    def __init__(
        self,
        chain: AtomGroup,
        patches: AtomGroup,
        *,
        cutoff: float,
        topology: TopologyT,
        inclusive: bool = True,
        cis_threshold: int = 4,
        patches_per_binder: int = _TWO_PATCHES,
        periodic: tuple[bool, bool, bool] = (False, False, False),
        store_distances: bool = False,
        chain_label: EntityT = 'Mon',
        patch_label: EntityT = 'HnsPatch',
        binder_label: EntityT = 'HnsCore',
        verbose: bool = False,
    ) -> None:
        super().__init__(chain.universe.trajectory, verbose=verbose)
        if chain.n_atoms < 2:
            raise ValueError(
                f"'chain' must hold at least two monomers, not "
                f"{chain.n_atoms}.")
        if patches.universe is not chain.universe:
            raise ValueError(
                "'patches' must come from the universe of 'chain'.")
        if patches_per_binder != _TWO_PATCHES:
            raise ValueError(
                f"'patches_per_binder' must be {_TWO_PATCHES}, not "
                f"{patches_per_binder!r}: binding_statistics counts "
                "two-patch binders only.")
        if patches.n_atoms < 1 or patches.n_atoms % patches_per_binder:
            raise ValueError(
                f"'patches' must hold a positive multiple of "
                f"{patches_per_binder} atoms, not {patches.n_atoms}.")
        if not (isinstance(cutoff, (int, float))
                and not isinstance(cutoff, bool)
                and math.isfinite(cutoff) and cutoff > 0):
            raise ValueError(
                f"'cutoff' must be a finite positive number, not "
                f"{cutoff!r}.")
        if topology not in get_args(Topology):
            raise ValueError(
                f"'topology' must be one of {list(get_args(Topology))}, "
                f"not {topology!r}.")
        if cis_threshold < 0:
            raise ValueError(
                f"'cis_threshold' must be non-negative, not "
                f"{cis_threshold!r}.")
        periodic = tuple(periodic)  # type: ignore[assignment]
        if len(periodic) != 3 \
                or any(not isinstance(flag, (bool, np.bool_))
                       for flag in periodic):
            raise ValueError(
                f"'periodic' must be three booleans, not {periodic!r}.")
        for name, label in (('chain_label', chain_label),
                            ('patch_label', patch_label),
                            ('binder_label', binder_label)):
            if not (isinstance(label, str) and label):
                raise ValueError(
                    f"'{name}' must be a non-empty species label, not "
                    f"{label!r}.")
        self.chain = chain
        self.patches = patches
        self.cutoff = float(cutoff)
        self.topology = topology
        self.inclusive = bool(inclusive)
        self.cis_threshold = int(cis_threshold)
        self.patches_per_binder = patches_per_binder
        self.periodic = cast(tuple[bool, bool, bool],
                             tuple(bool(flag) for flag in periodic))
        self.store_distances = bool(store_distances)
        self.chain_label = chain_label
        self.patch_label = patch_label
        self.binder_label = binder_label
        self.n_lengths = BindingStats(
            chain.n_atoms, cast(Topology, topology),
            self.cis_threshold).loop_length_histogram.size

    @classmethod
    def from_spec(
        cls,
        universe: Universe,
        spec: ProbeSpec,
        *,
        verbose: bool = False,
    ) -> Self:
        """
        Build the analysis of a study's binders from its probe spec.

        The chain is ``universe.select_atoms(spec.selections[spec.chain])``
        and the patches
        ``universe.select_atoms(spec.selections[spec.binders.patches])``;
        the cutoff, the comparison, the cis threshold, the patches per
        binder, whether the distances are stored, and the labels are the
        spec's :class:`~polyphys.spec.BinderSpec`; the topology
        is the project's; and the minimum image is taken along every
        direction in a cube and along ``z`` only in a cylinder.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        universe : :class:`~MDAnalysis.core.universe.Universe`
            The trajectory, as :func:`~polyphys.probe.core.make_universe`
            opens it.
        spec : :class:`~polyphys.spec.ProbeSpec`
            The study's probe spec, with `binders` set.
        verbose : bool, default False
            Whether to show a progress bar during :meth:`run`.

        Returns
        -------
        BinderStatistics
            The analysis, not yet run.

        Raises
        ------
        ValueError
            If the spec has no `binders`, or if the selections pick fewer
            than two monomers or a number of patches that is not a positive
            multiple of the patches per binder.

        Examples
        --------
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> from polyphys.manage.parser import HnsCub
        >>> from polyphys.spec import (
        ...     BinderSpec, ProbeSpec, ProjectSpec)
        >>> project = ProjectSpec(
        ...     parser=HnsCub, group='nucleoid', divisor=0.01,
        ...     physical_attrs=('space', 'nmon', 'nhns'),
        ...     species_sizes={'Mon': 'dmon', 'Hns': 'dhns'},
        ...     tseries_properties=(
        ...         'nBoundTHnsPatch', 'nFreeTHnsPatch', 'nEngagedTHnsPatch',
        ...         'nFreeTHnsCore', 'nBridgeTHnsCore', 'nDangleTHnsCore',
        ...         'nCisTHnsCore', 'nTransTHnsCore', 'gyrTMon'))
        >>> spec = ProbeSpec(
        ...     project, {'Mon': 'type 1', 'HnsPatch': 'type 2'}, 'Mon',
        ...     ('gyrT',), unwrap=True,
        ...     binders=BinderSpec('HnsPatch', 'HnsCore', cutoff=0.685))
        >>> universe = mda.Universe.empty(6, trajectory=True)
        >>> universe.add_TopologyAttr('types', ['1'] * 4 + ['2'] * 2)
        >>> analysis = BinderStatistics.from_spec(universe, spec)
        >>> analysis.chain.n_atoms, analysis.patches.n_atoms
        (4, 2)
        >>> analysis.topology, analysis.periodic, analysis.cutoff
        ('ring', (True, True, True), 0.685)
        """
        if spec.binders is None:
            raise ValueError(
                f"{spec.name}: the probe spec has no 'binders'.")
        binders = spec.binders
        chain = universe.select_atoms(spec.selections[spec.chain])
        patches = universe.select_atoms(spec.selections[binders.patches])
        return cls(chain, patches, cutoff=binders.cutoff,
                   topology=spec.project.topology,
                   inclusive=binders.inclusive,
                   cis_threshold=binders.cis_threshold,
                   patches_per_binder=binders.patches_per_binder,
                   periodic=_PERIODIC[spec.geometry],
                   store_distances=binders.store_distances,
                   chain_label=spec.chain, patch_label=binders.patches,
                   binder_label=binders.binders, verbose=verbose)

    @classmethod
    def get_supported_backends(cls) -> tuple[str, ...]:
        """
        The backends :meth:`run` accepts: ``'serial'``,
        ``'multiprocessing'``, and ``'dask'``.

        Examples
        --------
        >>> BinderStatistics.get_supported_backends()
        ('serial', 'multiprocessing', 'dask')
        """
        return _BACKENDS

    @property
    def keys(self) -> tuple[str, ...]:
        """
        The keys of ``results``, in the order :meth:`save` writes them:
        the eight series of
        :data:`~polyphys.spec.BINDER_SERIES`, the loop-length
        histogram, and, with `store_distances`, the distances.

        Examples
        --------
        >>> import MDAnalysis as mda
        >>> universe = mda.Universe.empty(4, trajectory=True)
        >>> analysis = BinderStatistics(
        ...     universe.atoms[:2], universe.atoms[2:], cutoff=0.5,
        ...     topology='linear')
        >>> analysis.keys[-2:]
        ('n_pairs_trans', 'loop_length_histogram')
        """
        keys = (*BINDER_SERIES, LOOP_KEY)
        if self.store_distances:
            keys = (*keys, DISTANCES_KEY)
        return keys

    def _get_aggregator(self) -> ResultsGroup:
        """
        Concatenate each worker's arrays along the frame axis.
        """
        return ResultsGroup(lookup={
            key: (ResultsGroup.ndarray_hstack if key in BINDER_SERIES
                  else ResultsGroup.ndarray_vstack)
            for key in self.keys})

    def _prepare(self) -> None:
        """
        Allocate the arrays for the frames of this run or worker: int64
        series and histogram rows, float64 distances.
        """
        for attr in BINDER_SERIES:
            self.results[attr] = np.empty(self.n_frames, dtype=np.int64)
        self.results[LOOP_KEY] = np.zeros(
            (self.n_frames, self.n_lengths), dtype=np.int64)
        if self.store_distances:
            self.results[DISTANCES_KEY] = np.empty(
                (self.n_frames, self.chain.n_atoms, self.patches.n_atoms),
                dtype=np.float64)

    def _single_frame(self) -> None:
        """
        Form the distances of the current frame, mark the contacts, count
        them, and store the counts and the histogram row for this frame.
        """
        if any(self.periodic):
            dimensions = self._ts.dimensions
            if dimensions is None:
                raise ValueError(
                    f"frame {self._ts.frame} carries no box, which the "
                    f"periodic directions {self.periodic} need.")
            box_lengths = np.asarray(dimensions[:3], dtype=np.float64)
        else:
            box_lengths = np.zeros(3)
        distances = _pair_distances(
            np.asarray(self.chain.positions, dtype=np.float64),
            np.asarray(self.patches.positions, dtype=np.float64),
            box_lengths, self.periodic)
        contacts = find_direct_contacts(distances, self.cutoff,
                                        inclusive=self.inclusive)
        frame = binding_statistics(contacts, cast(Topology, self.topology),
                                   self.cis_threshold)
        index = self._frame_index
        for attr in BINDER_SERIES:
            self.results[attr][index] = getattr(frame, attr)
        self.results[LOOP_KEY][index] = np.bincount(
            frame.bridged_pairs[:, 2], minlength=self.n_lengths)
        if self.store_distances:
            self.results[DISTANCES_KEY][index] = distances

    def save(self, save_to: Pathish, name: str) -> list[Path]:
        """
        Save each array as ``<name>-<series>.npy``, under the research
        code's names.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        save_to : str or pathlib.Path
            An existing directory.
        name : str
            The artifact's name, ``<parser.name>-<parser.group>`` in the
            studies, such as ``'N200...ens1-nucleoid'``.

        Returns
        -------
        list of pathlib.Path
            The files written, in the order of :attr:`keys`: the eight
            series as ``<stem><patch_label>`` or ``<stem><binder_label>``
            (``nBoundTHnsPatch``, ``nFreeTHnsPatch``, ``nEngagedTHnsPatch``,
            ``nFreeTHnsCore``, ``nBridgeTHnsCore``, ``nDangleTHnsCore``,
            ``nCisTHnsCore``, ``nTransTHnsCore`` with the default labels),
            the histogram as ``loopLengthHistT<chain_label>``, and, with
            `store_distances`, the distances as
            ``distMatT<chain_label><patch_label>``.

        Raises
        ------
        RuntimeError
            If :meth:`run` has not been called.

        Notes
        -----
        The research code summed the loop-length histogram over the frames
        into ``loopLengthHist<chain_label>``; the per-frame rows saved here
        sum to it along the frame axis.

        Examples
        --------
        >>> import tempfile
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> universe = mda.Universe.empty(4, trajectory=True)
        >>> universe.atoms.positions = np.array(
        ...     [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        ...      [0.0, 0.3, 0.0], [5.0, 5.0, 5.0]])
        >>> analysis = BinderStatistics(
        ...     universe.atoms[:2], universe.atoms[2:], cutoff=0.5,
        ...     topology='linear').run()
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     paths = analysis.save(tmp, 'N2ens1-nucleoid')
        ...     names = [path.name for path in paths]
        ...     dangling = np.load(paths[5]).tolist()
        >>> names[0]
        'N2ens1-nucleoid-nBoundTHnsPatch.npy'
        >>> names[5], dangling
        ('N2ens1-nucleoid-nDangleTHnsCore.npy', [1])
        >>> names[8]
        'N2ens1-nucleoid-loopLengthHistTMon.npy'
        """
        if any(key not in self.results for key in self.keys):
            raise RuntimeError("Call run() before save().")
        labels = {'patches': self.patch_label, 'binders': self.binder_label}
        stems = {attr: f"{stem}{labels[which]}"
                 for attr, (stem, which) in BINDER_SERIES.items()}
        stems[LOOP_KEY] = f"loopLengthHistT{self.chain_label}"
        stems[DISTANCES_KEY] = f"distMatT{self.chain_label}{self.patch_label}"
        paths = []
        for key in self.keys:
            path = Path(save_to) / f"{name}-{stems[key]}.npy"
            np.save(path, self.results[key])
            paths.append(path)
        return paths
