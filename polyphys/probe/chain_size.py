"""
==========================================================
:mod:`polyphys.probe.chain_size`
==========================================================

The :mod:`polyphys.probe.chain_size` module measures the size and shape of
one polymer chain in every frame of a trajectory: its end-to-end distance,
its span along the cylinder axis, its transverse size, its radius of
gyration, and the asphericity, shape parameter, and principal axes of its
gyration tensor. :class:`ChainSize` is an MDAnalysis
:class:`~MDAnalysis.analysis.base.AnalysisBase`, so it runs over a slice of
frames in serial or, split across workers, in parallel, and stores one
array per observable in its ``results``.

The research code measured the same seven quantities in one loop over the
frames of each study's probe function, after shifting every position by
half the box. Every quantity here is invariant under such a shift, because
each one first recenters the chain on its center of geometry or of mass,
so the shift is dropped.

Contents
========

- :class:`ChainSize`, the per-frame chain-size analysis.

Dependencies
============
- Python standard library: `pathlib`, `typing`.
- `NumPy <https://numpy.org/>`_.
- `MDAnalysis <https://docs.mdanalysis.org/>`_ 2.8 or later, through the
  ``mda`` extra: :class:`~MDAnalysis.analysis.base.AnalysisBase`,
  :class:`~MDAnalysis.analysis.results.ResultsGroup`, and the shape methods
  of :class:`~MDAnalysis.core.groups.AtomGroup`.
- :mod:`polyphys.analyze.measurer`: :func:`~polyphys.analyze.measurer.fsd`,
  :func:`~polyphys.analyze.measurer.end_to_end`, and
  :func:`~polyphys.analyze.measurer.transverse_size`.
- :mod:`polyphys.spec`:
  :data:`~polyphys.spec.CHAIN_SIZE_OBSERVABLES` and
  :class:`~polyphys.spec.ProbeSpec`.

Notation
========
:math:`N`
    The number of beads in the chain, ``atomgroup.n_atoms``.
:math:`\\mathbf{r}_i`, :math:`m_i`
    The position and mass of bead :math:`i`, in the dump's length and mass
    units; the studies use the monomer diameter and mass.
:math:`\\mathbf{r}_\\mathrm{cg}`, :math:`\\mathbf{r}_\\mathrm{cm}`
    The center of geometry :math:`\\sum_i \\mathbf{r}_i / N` and the center
    of mass :math:`\\sum_i m_i \\mathbf{r}_i / \\sum_i m_i`.
:math:`\\mathbf{S}`, :math:`\\lambda_1 \\ge \\lambda_2 \\ge \\lambda_3`
    The gyration tensor
    :math:`S_{ab} = \\sum_i m_i (r_{i,a} - r_{\\mathrm{cm},a})
    (r_{i,b} - r_{\\mathrm{cm},b}) / \\sum_i m_i`, and its eigenvalues,
    which MDAnalysis's ``gyration_moments`` returns.

Examples
========
A straight chain of four beads one unit apart, in a one-frame universe:

>>> import MDAnalysis as mda
>>> import numpy as np
>>> universe = mda.Universe.empty(4, trajectory=True)
>>> universe.add_TopologyAttr('masses', np.ones(4))
>>> universe.atoms.positions = np.array(
...     [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
>>> analysis = ChainSize(universe.atoms, ('rfloryT', 'asphericityT')).run()
>>> analysis.results['rfloryT'], analysis.results['asphericityT']
(array([3.]), array([1.]))

.. versionadded:: 0.7.0
"""
from pathlib import Path
from typing import Any, Callable, Mapping, Self, Sequence

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.core.universe import Universe

from ..analyze import measurer
from ..spec import CHAIN_SIZE_OBSERVABLES, ProbeSpec
from ..manage.types import AxisT, EntityT, Pathish, PropertyT

# One frame's value of an observable: a function of the chain, the cylinder
# axis, and whether to unwrap the chain.
_Measure = Callable[[AtomGroup, AxisT, bool], Any]


def _trans_size(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return measurer.transverse_size(group.positions, axis)


def _fsd(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return measurer.fsd(group.positions, axis)


def _gyr(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return group.radius_of_gyration()


def _rflory(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return measurer.end_to_end(group.positions)


def _asphericity(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return group.asphericity(wrap=False, unwrap=unwrap)


def _shape(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return group.shape_parameter(wrap=False)


def _principal(group: AtomGroup, axis: AxisT, unwrap: bool) -> Any:
    return group.principal_axes(wrap=False)


# The measure of each observable, keyed as CHAIN_SIZE_OBSERVABLES names them.
_MEASURES: Mapping[PropertyT, _Measure] = {
    'transSizeT': _trans_size,
    'fsdT': _fsd,
    'gyrT': _gyr,
    'rfloryT': _rflory,
    'asphericityT': _asphericity,
    'shapeT': _shape,
    'principalT': _principal,
}

# The shape of one frame's value, after the frame axis; a scalar otherwise.
_FRAME_SHAPES: Mapping[PropertyT, tuple[int, ...]] = {'principalT': (3, 3)}

_BACKENDS = ('serial', 'multiprocessing', 'dask')


class ChainSize(AnalysisBase):
    """
    The size and shape of one polymer chain, measured in every frame.

    .. versionadded:: 0.7.0

    Parameters
    ----------
    atomgroup : :class:`~MDAnalysis.core.groups.AtomGroup`
        The chain, at least two beads, in atom-ID order: the order along
        the backbone, which ``Universe.select_atoms`` gives.
    observables : sequence of :py:data:`PropertyT`
        The observables to measure, each a name in
        :data:`~polyphys.spec.CHAIN_SIZE_OBSERVABLES`, in the
        order ``results`` and :meth:`save` list them.
    chain : :py:data:`EntityT`, default 'Mon'
        The species label of the chain, which names the files
        :meth:`save` writes, ``<name>-<observable><chain>.npy``.
    unwrap : bool, default False
        Whether the chain is unwrapped across the periodic boundaries
        before its gyration tensor is formed, for ``asphericityT``. The
        research code unwrapped in a periodic cube and not in a cylinder.
    axis : :py:data:`AxisT`, default 2
        The cylinder axis, ``0``, ``1``, or ``2`` for ``x``, ``y``, or
        ``z``, along which ``fsdT`` spans and across which ``transSizeT``
        spreads.
    verbose : bool, default False
        Whether to show a progress bar during :meth:`run`.

    Attributes
    ----------
    atomgroup : :class:`~MDAnalysis.core.groups.AtomGroup`
        The chain.
    observables : tuple of :py:data:`PropertyT`
        The observables, in order.
    chain : :py:data:`EntityT`
        The species label of the chain.
    unwrap : bool
        Whether the chain is unwrapped for ``asphericityT``.
    axis : :py:data:`AxisT`
        The cylinder axis.
    results : :class:`~MDAnalysis.analysis.results.Results`
        After :meth:`run`, one float64 array per observable, keyed by its
        name: shape ``(n_frames,)`` for every scalar, ``(n_frames, 3, 3)``
        for ``principalT``.
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
        If `observables` is empty, repeats a name, or names an observable
        that is not in
        :data:`~polyphys.spec.CHAIN_SIZE_OBSERVABLES`; if
        `chain` is empty; if `axis` is not ``0``, ``1``, or ``2``; or if
        `atomgroup` has fewer than two beads.

    Notes
    -----
    The observables, with the species label ``Mon`` as the studies name
    them:

    ``rfloryTMon``
        :math:`|\\mathbf{r}_N - \\mathbf{r}_1|`, the end-to-end distance
        (:func:`~polyphys.analyze.measurer.end_to_end`).
    ``fsdTMon``
        :math:`\\max_i z_i - \\min_i z_i` along `axis`, the span or Feret
        diameter (:func:`~polyphys.analyze.measurer.fsd`).
    ``transSizeTMon``
        :math:`2 \\langle |\\mathbf{r}_{i,\\perp} -
        \\mathbf{r}_{\\mathrm{cg},\\perp}| \\rangle_i`, twice the mean
        distance from the center of geometry in the plane across `axis`
        (:func:`~polyphys.analyze.measurer.transverse_size`).
    ``gyrTMon``
        :math:`R_g = \\sqrt{\\lambda_1 + \\lambda_2 + \\lambda_3}`, the
        mass-weighted radius of gyration
        (``AtomGroup.radius_of_gyration``).
    ``asphericityTMon``
        :math:`\\frac{3}{2} \\sum_k (\\lambda_k - \\bar\\lambda)^2 /
        (\\sum_k \\lambda_k)^2`, from 0 for a spherically symmetric chain
        to 1 for a straight one
        (``AtomGroup.asphericity(wrap=False, unwrap=unwrap)``).
    ``shapeTMon``
        :math:`27 \\prod_k (\\lambda_k - \\bar\\lambda) /
        (\\sum_k \\lambda_k)^3`, from 0 for a spherically symmetric chain
        to 2 for a straight one (``AtomGroup.shape_parameter(wrap=False)``).
    ``principalTMon``
        The :math:`3 \\times 3` matrix of principal axes of the moment of
        inertia, one eigenvector per row, by decreasing eigenvalue
        (``AtomGroup.principal_axes(wrap=False)``).

    The asphericity and shape parameter follow Dima and Thirumalai [1]_,
    as MDAnalysis implements them. ``shapeT`` and ``principalT`` read the
    positions as the dump wraps them, whatever `unwrap` is, as the research
    code did; `unwrap` reaches ``asphericityT`` only.

    No per-frame shift of the positions is applied: each observable
    recenters the chain itself, so the research code's shift of every
    position by half the box changed none of these values.

    The MDAnalysis contract: :meth:`_prepare` allocates one float64 array
    per observable for the frames of the run, :meth:`_single_frame` writes
    the current frame's row, and :meth:`_conclude` does nothing. The class
    is marked parallelizable, with ``serial``, ``multiprocessing``, and
    ``dask`` backends; each worker fills the arrays of its own frames,
    which :meth:`_get_aggregator` concatenates along the frame axis
    (``ResultsGroup.ndarray_hstack`` for a scalar series,
    ``ResultsGroup.ndarray_vstack`` for ``principalT``).

    Cost per frame: :math:`O(N)` time for every observable, since each is a
    sum over the beads followed by at most a :math:`3 \\times 3`
    eigenproblem, plus the :math:`O(N + N_\\mathrm{bonds})` unwrap of the
    chain when `unwrap` is set; :math:`O(N)` extra space for the recentered
    positions. Over a run, :math:`O(F)` space per scalar series and
    :math:`O(9 F)` for ``principalT``, for :math:`F` frames. Positions come
    from MDAnalysis as float32, so each value carries a relative round-off
    of about :math:`10^{-7}`.

    References
    ----------
    .. [1] R. I. Dima and D. Thirumalai, "Asymmetry in the Shapes of Folded
       and Denatured States of Proteins", J. Phys. Chem. B 108, 6564-6570
       (2004). https://doi.org/10.1021/jp037128y

    Examples
    --------
    Two frames of a three-bead chain, the second frame stretched along z:

    >>> import MDAnalysis as mda
    >>> import numpy as np
    >>> frames = np.array([
    ...     [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
    ...     [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]]])
    >>> universe = mda.Universe.empty(3, trajectory=True)
    >>> universe.add_TopologyAttr('masses', np.ones(3))
    >>> universe.load_new(frames, format='MEMORY')
    <Universe with 3 atoms>
    >>> analysis = ChainSize(universe.atoms, ('rfloryT', 'fsdT', 'gyrT'),
    ...                      axis=2)
    >>> analysis = analysis.run()
    >>> np.round(analysis.results['rfloryT'], 6)
    array([1.414214, 2.      ])
    >>> analysis.results['fsdT']
    array([0., 2.])

    The frames to analyze are chosen as in any MDAnalysis analysis:

    >>> ChainSize(universe.atoms, ('fsdT',)).run(start=1).results['fsdT']
    array([2.])
    """

    _analysis_algorithm_is_parallelizable = True

    def __init__(
        self,
        atomgroup: AtomGroup,
        observables: Sequence[PropertyT],
        *,
        chain: EntityT = 'Mon',
        unwrap: bool = False,
        axis: AxisT = 2,
        verbose: bool = False,
    ) -> None:
        super().__init__(atomgroup.universe.trajectory, verbose=verbose)
        observables = tuple(observables)
        if not observables:
            raise ValueError("'observables' must not be empty.")
        unknown = [obs for obs in observables if obs not in _MEASURES]
        if unknown:
            raise ValueError(
                f"'observables' names {unknown}, which are not in "
                f"{list(CHAIN_SIZE_OBSERVABLES)}.")
        if len(set(observables)) != len(observables):
            raise ValueError(
                f"'observables' repeats a name: {list(observables)}.")
        if not (isinstance(chain, str) and chain):
            raise ValueError(
                f"'chain' must be a non-empty species label, not {chain!r}.")
        if axis not in (0, 1, 2):
            raise ValueError(f"'axis' must be 0, 1, or 2, not {axis!r}.")
        if atomgroup.n_atoms < 2:
            raise ValueError(
                f"'atomgroup' must hold at least two beads, not "
                f"{atomgroup.n_atoms}.")
        self.atomgroup = atomgroup
        self.observables = observables
        self.chain = chain
        self.unwrap = unwrap
        self.axis = axis

    @classmethod
    def from_spec(
        cls,
        universe: Universe,
        spec: ProbeSpec,
        *,
        axis: AxisT = 2,
        verbose: bool = False,
    ) -> Self:
        """
        Build the analysis of a study's chain from its probe spec.

        The chain is ``universe.select_atoms(spec.selections[spec.chain])``;
        the observables, the chain label, and `unwrap` are the spec's.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        universe : :class:`~MDAnalysis.core.universe.Universe`
            The trajectory, as :func:`~polyphys.probe.core.make_universe`
            opens it.
        spec : :class:`~polyphys.spec.ProbeSpec`
            The study's probe spec.
        axis : :py:data:`AxisT`, default 2
            The cylinder axis.
        verbose : bool, default False
            Whether to show a progress bar during :meth:`run`.

        Returns
        -------
        ChainSize
            The analysis, not yet run.

        Raises
        ------
        ValueError
            If the selection picks fewer than two beads.

        Examples
        --------
        >>> import MDAnalysis as mda
        >>> import numpy as np
        >>> from polyphys.manage.parser import SumRuleCyl
        >>> from polyphys.spec import ProbeSpec, ProjectSpec
        >>> project = ProjectSpec(
        ...     parser=SumRuleCyl, group='bug', divisor=0.025,
        ...     physical_attrs=('space', 'nmon'),
        ...     species_sizes={'Mon': 'dmon', 'Crd': 'dcrowd'},
        ...     tseries_properties=('gyrTMon', 'rfloryTMon'))
        >>> spec = ProbeSpec(project, {'Mon': 'resid 1'}, 'Mon',
        ...                  ('rfloryT',), unwrap=False)
        >>> universe = mda.Universe.empty(3, n_residues=1,
        ...                               atom_resindex=[0, 0, 0],
        ...                               trajectory=True)
        >>> universe.add_TopologyAttr('resid', [1])
        >>> universe.add_TopologyAttr('masses', np.ones(3))
        >>> universe.atoms.positions = np.array(
        ...     [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        >>> analysis = ChainSize.from_spec(universe, spec)
        >>> analysis.chain, analysis.atomgroup.n_atoms
        ('Mon', 3)
        >>> analysis.run().results['rfloryT']
        array([2.])
        """
        chain = universe.select_atoms(spec.selections[spec.chain])
        return cls(chain, spec.observables, chain=spec.chain,
                   unwrap=spec.unwrap, axis=axis, verbose=verbose)

    @classmethod
    def get_supported_backends(cls) -> tuple[str, ...]:
        """
        The backends :meth:`run` accepts: ``'serial'``,
        ``'multiprocessing'``, and ``'dask'``.

        Examples
        --------
        >>> ChainSize.get_supported_backends()
        ('serial', 'multiprocessing', 'dask')
        """
        return _BACKENDS

    def _get_aggregator(self) -> ResultsGroup:
        """
        Concatenate each worker's arrays along the frame axis.
        """
        return ResultsGroup(lookup={
            obs: (ResultsGroup.ndarray_vstack if obs in _FRAME_SHAPES
                  else ResultsGroup.ndarray_hstack)
            for obs in self.observables})

    def _prepare(self) -> None:
        """
        Allocate one float64 array per observable for the frames of this
        run or worker.
        """
        for obs in self.observables:
            self.results[obs] = np.empty(
                (self.n_frames, *_FRAME_SHAPES.get(obs, ())),
                dtype=np.float64)

    def _single_frame(self) -> None:
        """
        Measure every observable on the current frame and store it in its
        array's row for this frame.
        """
        for obs in self.observables:
            self.results[obs][self._frame_index] = _MEASURES[obs](
                self.atomgroup, self.axis, self.unwrap)

    def save(self, save_to: Pathish, name: str) -> list[Path]:
        """
        Save each observable's array as ``<name>-<observable><chain>.npy``.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        save_to : str or pathlib.Path
            An existing directory.
        name : str
            The artifact's name, ``<parser.name>-<parser.group>`` in the
            studies, such as ``'N200...ens1-bug'``.

        Returns
        -------
        list of pathlib.Path
            The files written, one per observable, in order.

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
        >>> universe.add_TopologyAttr('masses', np.ones(2))
        >>> universe.atoms.positions = np.array([[0.0, 0.0, 0.0],
        ...                                      [0.0, 0.0, 3.0]])
        >>> analysis = ChainSize(universe.atoms, ('rfloryT', 'fsdT')).run()
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     paths = analysis.save(tmp, 'N2ens1-bug')
        ...     [path.name for path in paths], np.load(paths[1]).tolist()
        (['N2ens1-bug-rfloryTMon.npy', 'N2ens1-bug-fsdTMon.npy'], [3.0])
        """
        if any(obs not in self.results for obs in self.observables):
            raise RuntimeError("Call run() before save().")
        paths = []
        for obs in self.observables:
            path = Path(save_to) / f"{name}-{obs}{self.chain}.npy"
            np.save(path, self.results[obs])
            paths.append(path)
        return paths
