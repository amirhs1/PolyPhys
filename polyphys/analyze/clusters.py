"""
==========================================================
:mod:`polyphys.analyze.clusters`
==========================================================


The :mod:`polyphys.analyze.clusters` module turns the positions of a group of
particles in one frame into contact and cluster statistics, and the contacts
between a polymer and two-patch binding proteins into binding statistics.
Every routine is a pure numpy or scipy kernel on one frame, except
:func:`threshold_sweep`, which recounts saved frames over a grid of cutoffs;
reading frames, choosing the particles, and accumulating over frames are
otherwise the caller's work.

Notation
========

``n_particles``
    The number of particles in the group, :math:`n`.
``positions``
    Their coordinates, shape ``(n_particles, n_dims)``, in the caller's length
    unit; every length below is in that unit.
``dist``
    The matrix of centre-to-centre distances, shape
    ``(n_particles, n_particles)``, symmetric with a zero diagonal.
``cutoff``
    A length :math:`r_c`. Two particles are in *direct contact* when their
    centre-to-centre distance is at most :math:`r_c` (``inclusive=True``) or
    below it (``inclusive=False``); this is the contact criterion of Sevick,
    Monson, and Ottino (1988).
``direct_contacts``
    The boolean *direct-contact matrix* :math:`A`, ``A[i, j]`` true when
    particles ``i`` and ``j`` are in direct contact. For a self-distance
    matrix it is symmetric and its diagonal is true, since every particle is
    at distance zero from itself. It is the adjacency matrix of the *contact
    graph*, whose vertices are the particles and whose edges are the direct
    contacts; ``n_contacts`` is its number of edges, :math:`E`.
*cluster*
    A connected component of the contact graph: a maximal set of particles
    joined by chains of direct contacts. Its *size* is its number of
    particles, from 1 to ``n_particles``.
``contact_matrix``
    The boolean *contact matrix* :math:`C`, ``C[i, j]`` true when ``i`` and
    ``j`` are in the same cluster: the transitive closure of :math:`A`, block
    all-ones after sorting the particles by cluster.
*bonds*
    The number of other particles a particle is in direct contact with: its
    degree in the contact graph, from 0 to ``n_particles - 1``.
*genomic distance*
    For two monomers of one polymer, the number of backbone bonds between them
    along the chain: ``|i - j|`` on a linear chain, and the shorter way round,
    ``min(|i - j|, n_backbone - |i - j|)``, on a ring of ``n_backbone``
    monomers.
*patch*, *binder*
    A binder is a protein with two binding sites, its *patches*; binder ``k``
    owns patches ``2k`` and ``2k + 1``. A binder is *free* when neither patch
    touches a monomer, *dangling* when one does, and *bridging* when both do.
    A bridge joins two monomers; it is *cis* when their genomic distance is at
    most ``cis_threshold`` bonds and *trans* when it is larger.

Contents
========

Distances and contacts:
:func:`pair_distances` builds ``dist`` from ``positions``, with the minimum
image in the periodic directions of an orthogonal box;
:func:`find_direct_contacts` applies the cutoff.

Clusters:
:func:`cluster_labels` labels the connected components of the contact graph
with :func:`scipy.sparse.csgraph.connected_components`;
:func:`count_clusters` histograms their sizes; :func:`contact_matrix` expands
the labels to :math:`C`; :func:`count_bonds` histograms the degrees.

Genomic distances and pair separations:
:func:`genomic_distance` and :func:`foci_pair_info` measure distances along
the backbone; :func:`pair_distance_histograms` bins the time series of pair
separations on one shared grid, and :func:`pair_separation_density` turns the
counts into the probability density of the separation vector.

Binding:
:func:`enforce_single_patch_contact`, :func:`monomer_binder_contacts`,
:func:`split_binders`, and :func:`binder_contacts` reshape a monomer-patch
contact matrix; :func:`binding_statistics` counts free, dangling, bridging,
cis, and trans binders in one frame as a :class:`FrameBinding`,
:class:`BindingStats` accumulates frames and the loop-length histogram, and
:func:`threshold_sweep` recounts the saved distance matrices of a whole at
every cutoff and cis threshold of a grid, the research code's
choice-of-threshold sweep.

Conventions
===========

- A contact matrix is boolean. Integer 0/1 input is accepted everywhere and
  converted; any nonzero entry counts as a contact.
- Cluster labels are integers ``0, ..., n_clusters - 1``; which cluster gets
  which label is an implementation detail of SciPy, so compare partitions
  (``labels[i] == labels[j]``), not label values.
- Size and degree histograms are indexed by the size or degree itself:
  ``count_clusters(A)[s]`` is the number of clusters of size ``s`` and has
  length ``n_particles + 1`` with ``[0]`` always 0; ``count_bonds(A)[d]`` is
  the number of particles with ``d`` bonds and has length ``n_particles``.
- Genomic distances are in bonds, never in monomers between.
- Contact cutoffs, cis thresholds, and bin ranges are project constants and
  have no defaults here.
- Complexities are stated per call, for one frame, in ``n_particles``
  (:math:`n`), ``n_contacts`` (:math:`E`), and the other sizes named in each
  docstring.

References
==========

- E. M. Sevick, P. A. Monson, and J. M. Ottino, "Monte Carlo calculations of
  cluster statistics in continuum models of composite morphology",
  J. Chem. Phys. 88, 1198-1206 (1988). https://doi.org/10.1063/1.454720
- J. Hopcroft and R. Tarjan, "Algorithm 447: efficient algorithms for graph
  manipulation", Commun. ACM 16, 372-378 (1973).
  https://doi.org/10.1145/362248.362272
- SciPy, :func:`scipy.sparse.csgraph.connected_components`.
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.connected_components.html
- M. Rubinstein and R. H. Colby, Polymer Physics (Oxford University Press,
  2003), Chapter 2.

.. versionadded:: 0.6.0
"""

from dataclasses import dataclass, field
from typing import Literal, Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

__all__ = [
    'BindingStats',
    'FrameBinding',
    'Topology',
    'binder_contacts',
    'binding_statistics',
    'cluster_labels',
    'contact_matrix',
    'count_bonds',
    'count_clusters',
    'enforce_single_patch_contact',
    'find_direct_contacts',
    'foci_pair_info',
    'genomic_distance',
    'monomer_binder_contacts',
    'pair_distance_histograms',
    'pair_distances',
    'pair_separation_density',
    'split_binders',
    'threshold_sweep',
]

Topology = Literal['linear', 'ring']
"""The polymer topologies the genomic distance is defined for."""

# The eight counts of a FrameBinding, in the order BindingStats.as_dict and
# threshold_sweep report them.
_BINDING_COUNTS = (
    'n_patch_monomer_contacts', 'n_patches_engaged',
    'n_patches_free', 'n_binders_free', 'n_binders_dangling',
    'n_binders_bridging', 'n_pairs_cis', 'n_pairs_trans',
)


def _as_contacts(matrix: npt.ArrayLike, name: str) -> npt.NDArray[np.bool_]:
    """
    Copy `matrix` to a two-dimensional boolean array.

    Raises
    ------
    ValueError
        If `matrix` is not two-dimensional.
    """
    contacts = np.asarray(matrix)
    if contacts.ndim != 2:
        raise ValueError(
            f'{name} must be a two-dimensional matrix; got {contacts.ndim}'
            ' dimension(s).'
        )
    return contacts.astype(np.bool_, copy=False)


def _check_topology(topology: str) -> None:
    """Raise ``ValueError`` unless `topology` is ``'linear'`` or ``'ring'``."""
    if topology not in ('linear', 'ring'):
        raise ValueError(
            f"topology must be 'linear' or 'ring'; got '{topology}'."
        )


def pair_distances(
    positions: npt.ArrayLike,
    box: npt.ArrayLike | None = None,
) -> npt.NDArray[np.float64]:
    """
    Compute the matrix of centre-to-centre distances between all pairs.

    Parameters
    ----------
    positions : array_like
        Coordinates of the particles, shape ``(n_particles, n_dims)``, in the
        caller's length unit.
    box : array_like, optional
        Lengths of an orthogonal periodic box, shape ``(n_dims,)``, in the
        same unit. A direction with length ``0`` or ``inf`` is not periodic.
        With ``None``, no direction is periodic.

    Returns
    -------
    numpy.ndarray
        ``dist``, shape ``(n_particles, n_particles)``, ``float64``: the
        distance between particles ``i`` and ``j`` at ``[i, j]``; symmetric,
        with a zero diagonal.

    Raises
    ------
    ValueError
        If `positions` is not two-dimensional, if `box` has not one length
        per dimension, or if a box length is negative or NaN.

    Notes
    -----
    In a periodic direction of length :math:`L` the separation :math:`\\Delta`
    is replaced by its minimum image, :math:`\\Delta - L\\,
    \\mathrm{round}(\\Delta / L)`, which lies in :math:`[-L/2, L/2]`: the
    distance to the nearest periodic image of the other particle. It is the
    distance a contact criterion should use, provided the cutoff is below
    :math:`L / 2`; for wrapped coordinates it is the only way to recover a
    pair that straddles the boundary. Positions need not be wrapped.

    Time is :math:`O(n^2 d)` for :math:`n` particles in :math:`d` dimensions.
    Space is two ``(n, n)`` ``float64`` arrays; the dimensions are accumulated
    one at a time, so no ``(d, n, n)`` temporary is built.

    Examples
    --------
    A 3-4-5 right triangle:

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import pair_distances
    >>> pair_distances([[0.0, 0.0], [3.0, 0.0], [3.0, 4.0]])
    array([[0., 3., 5.],
           [3., 0., 4.],
           [5., 4., 0.]])

    Two particles at ``x = 0.1`` and ``x = 9.9`` are ``0.2`` apart through
    the boundary of a box of length ``10`` and ``9.8`` apart without it:

    >>> positions = [[0.1, 0.0], [9.9, 0.0]]
    >>> round(float(pair_distances(positions, box=[10.0, 0.0])[0, 1]), 10)
    0.2
    >>> round(float(pair_distances(positions)[0, 1]), 10)
    9.8

    .. versionadded:: 0.6.0
    """
    pos = np.asarray(positions, dtype=np.float64)
    if pos.ndim != 2:
        raise ValueError(
            'positions must have shape (n_particles, n_dims); got'
            f' {pos.ndim} dimension(s).'
        )
    n_particles, n_dims = pos.shape
    if box is None:
        lengths = np.zeros(n_dims)
    else:
        lengths = np.asarray(box, dtype=np.float64)
        if lengths.shape != (n_dims,):
            raise ValueError(
                f'box must have one length per dimension, shape ({n_dims},);'
                f' got shape {lengths.shape}.'
            )
        if np.any(np.isnan(lengths)) or np.any(lengths < 0):
            raise ValueError('box lengths must be non-negative numbers.')
    dist_sq = np.zeros((n_particles, n_particles), dtype=np.float64)
    for axis in range(n_dims):
        delta = pos[:, axis, None] - pos[None, :, axis]
        length = lengths[axis]
        if 0 < length < np.inf:
            delta -= length * np.round(delta / length)
        dist_sq += delta * delta
    return np.sqrt(dist_sq)


def find_direct_contacts(
    dist: npt.ArrayLike,
    cutoff: float,
    inclusive: bool = True,
) -> npt.NDArray[np.bool_]:
    """
    Mark the pairs whose distance is within a cutoff as direct contacts.

    Parameters
    ----------
    dist : array_like
        Centre-to-centre distances, shape ``(n_rows, n_cols)``, in the
        caller's length unit: a self-distance matrix from
        :func:`pair_distances`, or the distances between two groups, such as
        monomers along the rows and binder patches along the columns.
    cutoff : float
        The contact distance :math:`r_c`, in the same unit; positive.
    inclusive : bool, default True
        Count a pair at exactly the cutoff as a contact (``dist <= cutoff``)
        or not (``dist < cutoff``).

    Returns
    -------
    numpy.ndarray
        ``direct_contacts``, the same shape as `dist`, ``bool``. For a
        self-distance matrix it is symmetric with a true diagonal.

    Raises
    ------
    ValueError
        If `cutoff` is not positive or `dist` is not two-dimensional.

    Notes
    -----
    This is the contact criterion of Sevick, Monson, and Ottino (1988),
    J. Chem. Phys. 88, 1198, https://doi.org/10.1063/1.454720, for which the
    clusters are the connected components of the resulting graph. The cutoff
    is a project constant: in the foci projects it is the sum of a large
    monomer's and a crowder's diameter, the largest separation at which no
    crowder fits between the two.

    Time and space are :math:`O(n_\\mathrm{rows} n_\\mathrm{cols})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.clusters import find_direct_contacts
    >>> dist = np.array([[0.0, 1.0, 2.5],
    ...                  [1.0, 0.0, 1.5],
    ...                  [2.5, 1.5, 0.0]])
    >>> find_direct_contacts(dist, 1.5).astype(int)
    array([[1, 1, 0],
           [1, 1, 1],
           [0, 1, 1]])
    >>> find_direct_contacts(dist, 1.5, inclusive=False).astype(int)
    array([[1, 1, 0],
           [1, 1, 0],
           [0, 0, 1]])

    .. versionadded:: 0.6.0
    """
    if not cutoff > 0:
        raise ValueError(f'cutoff must be positive; got {cutoff}.')
    distances = np.asarray(dist, dtype=np.float64)
    if distances.ndim != 2:
        raise ValueError(
            f'dist must be a two-dimensional matrix; got {distances.ndim}'
            ' dimension(s).'
        )
    if inclusive:
        return np.asarray(distances <= cutoff, dtype=np.bool_)
    return np.asarray(distances < cutoff, dtype=np.bool_)


def cluster_labels(
    direct_contacts: npt.ArrayLike,
) -> tuple[int, npt.NDArray[np.int64]]:
    """
    Label the clusters, the connected components of the contact graph.

    Parameters
    ----------
    direct_contacts : array_like
        The direct-contact matrix :math:`A`, shape
        ``(n_particles, n_particles)``, boolean or 0/1, symmetric. The
        diagonal is ignored.

    Returns
    -------
    n_clusters : int
        The number of clusters, from 0 (no particles) to ``n_particles``.
    labels : numpy.ndarray
        Shape ``(n_particles,)``, ``int64``: ``labels[i]`` is the cluster of
        particle ``i``, in ``0, ..., n_clusters - 1``. Which cluster gets
        which label is not specified; compare ``labels[i] == labels[j]``.

    Raises
    ------
    ValueError
        If `direct_contacts` is not a square symmetric matrix.

    Notes
    -----
    Two particles are in the same cluster when a chain of direct contacts
    joins them, so the clusters are the connected components of the graph
    with adjacency matrix :math:`A`. They are found with
    :func:`scipy.sparse.csgraph.connected_components` (``directed=False``),
    a depth-first search that runs in time and space proportional to
    :math:`n + E` once the graph is sparse (Hopcroft and Tarjan, 1973,
    Commun. ACM 16, 372, https://doi.org/10.1145/362248.362272). Building the
    sparse graph from the dense matrix costs :math:`O(n^2)` time and
    :math:`O(E)` space, so the call is :math:`O(n^2)` in all; the Python-level
    merge of matrix columns it replaces was :math:`O(n^3)` and could stop
    short of the closure.

    Examples
    --------
    A path 0-1-2 and an isolated particle 3:

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import cluster_labels
    >>> contacts = np.array([[1, 1, 0, 0],
    ...                      [1, 1, 1, 0],
    ...                      [0, 1, 1, 0],
    ...                      [0, 0, 0, 1]])
    >>> n_clusters, labels = cluster_labels(contacts)
    >>> n_clusters
    2
    >>> bool(labels[0] == labels[2]), bool(labels[0] == labels[3])
    (True, False)

    .. versionadded:: 0.6.0
    """
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_particles = contacts.shape[0]
    if contacts.shape != (n_particles, n_particles):
        raise ValueError(
            'direct_contacts must be square; got shape'
            f' {contacts.shape}.'
        )
    if not np.array_equal(contacts, contacts.T):
        raise ValueError('direct_contacts must be symmetric.')
    if n_particles == 0:
        return 0, np.zeros(0, dtype=np.int64)
    n_clusters, labels = connected_components(
        csr_array(contacts.astype(np.int8)),
        directed=False,
        return_labels=True,
    )
    return int(n_clusters), np.asarray(labels, dtype=np.int64)


def count_clusters(direct_contacts: npt.ArrayLike) -> npt.NDArray[np.int64]:
    """
    Histogram the cluster sizes of one frame.

    Parameters
    ----------
    direct_contacts : array_like
        The direct-contact matrix :math:`A`, shape
        ``(n_particles, n_particles)``, boolean or 0/1, symmetric.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_particles + 1,)``, ``int64``: ``[s]`` is the number of
        clusters of size ``s``. ``[0]`` is always 0, since no cluster is
        empty, and ``sum(s * hist[s])`` is ``n_particles``.

    Raises
    ------
    ValueError
        If `direct_contacts` is not a square symmetric matrix.

    Notes
    -----
    The sizes are the multiplicities of the labels from
    :func:`cluster_labels`, counted with :func:`numpy.bincount` twice: once to
    get each cluster's size, once to histogram the sizes. This replaces
    reading the sizes off the eigenvalues of the contact matrix, which is
    :math:`O(n^3)` and meaningful only when that matrix is exactly block
    all-ones. Time is :math:`O(n^2)`, dominated by :func:`cluster_labels`;
    space is :math:`O(n + E)`.

    Examples
    --------
    A path 0-1-2 and an isolated particle 3 give one cluster of size 3 and
    one of size 1:

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import count_clusters
    >>> contacts = np.array([[1, 1, 0, 0],
    ...                      [1, 1, 1, 0],
    ...                      [0, 1, 1, 0],
    ...                      [0, 0, 0, 1]])
    >>> count_clusters(contacts)
    array([0, 1, 0, 1, 0])

    .. versionadded:: 0.6.0
    """
    _, labels = cluster_labels(direct_contacts)
    n_particles = labels.shape[0]
    sizes = np.bincount(labels)
    return np.bincount(sizes, minlength=n_particles + 1).astype(np.int64)


def contact_matrix(direct_contacts: npt.ArrayLike) -> npt.NDArray[np.bool_]:
    """
    Expand the direct contacts to the full contact matrix of the clusters.

    Parameters
    ----------
    direct_contacts : array_like
        The direct-contact matrix :math:`A`, shape
        ``(n_particles, n_particles)``, boolean or 0/1, symmetric.

    Returns
    -------
    numpy.ndarray
        The contact matrix :math:`C`, the same shape, ``bool``: ``[i, j]``
        is true when ``i`` and ``j`` are in the same cluster. Symmetric, with
        a true diagonal; it is its own transitive closure, so
        ``(C @ C > 0) == C``.

    Raises
    ------
    ValueError
        If `direct_contacts` is not a square symmetric matrix.

    Notes
    -----
    :math:`C` is the transitive closure of :math:`A` in the sense of Sevick,
    Monson, and Ottino (1988): ``C[i, j]`` is true when some chain of direct
    contacts joins ``i`` to ``j``. It is built as ``labels[i] == labels[j]``
    from :func:`cluster_labels`. Time is :math:`O(n^2)`; space is one
    ``(n, n)`` boolean array.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.clusters import contact_matrix
    >>> contacts = np.array([[1, 1, 0, 0],
    ...                      [1, 1, 1, 0],
    ...                      [0, 1, 1, 0],
    ...                      [0, 0, 0, 1]])
    >>> contact_matrix(contacts).astype(int)
    array([[1, 1, 1, 0],
           [1, 1, 1, 0],
           [1, 1, 1, 0],
           [0, 0, 0, 1]])

    .. versionadded:: 0.6.0
    """
    _, labels = cluster_labels(direct_contacts)
    return np.asarray(labels[:, None] == labels[None, :], dtype=np.bool_)


def count_bonds(direct_contacts: npt.ArrayLike) -> npt.NDArray[np.int64]:
    """
    Histogram the number of bonds, the direct contacts per particle.

    Parameters
    ----------
    direct_contacts : array_like
        The direct-contact matrix :math:`A`, shape
        ``(n_particles, n_particles)``, boolean or 0/1, symmetric. A true
        diagonal is allowed and does not count as a bond.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_particles,)``, ``int64``: ``[d]`` is the number of
        particles with ``d`` bonds, ``d`` from 0 to ``n_particles - 1``.
        ``sum(hist)`` is ``n_particles`` and ``sum(d * hist[d])`` is twice
        the number of direct contacts.

    Raises
    ------
    ValueError
        If `direct_contacts` is not a square symmetric matrix.

    Notes
    -----
    A particle's bonds are its degree in the contact graph: its row sum
    minus its diagonal entry. Time and space are :math:`O(n^2)`.

    Examples
    --------
    A star: particle 0 touches 1, 2, and 3, which touch nothing else, so one
    particle has three bonds and three have one.

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import count_bonds
    >>> star = np.array([[1, 1, 1, 1],
    ...                  [1, 1, 0, 0],
    ...                  [1, 0, 1, 0],
    ...                  [1, 0, 0, 1]])
    >>> count_bonds(star)
    array([0, 3, 0, 1])

    .. versionadded:: 0.6.0
    """
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_particles = contacts.shape[0]
    if contacts.shape != (n_particles, n_particles):
        raise ValueError(
            'direct_contacts must be square; got shape'
            f' {contacts.shape}.'
        )
    if not np.array_equal(contacts, contacts.T):
        raise ValueError('direct_contacts must be symmetric.')
    degrees = contacts.sum(axis=1) - np.diagonal(contacts).astype(np.int64)
    return np.bincount(degrees, minlength=n_particles).astype(np.int64)


def genomic_distance(
    atoms_a: npt.ArrayLike,
    atoms_b: npt.ArrayLike,
    topology: Topology,
    n_backbone: int,
) -> npt.NDArray[np.int64]:
    """
    Count the backbone bonds between paired monomers along the chain.

    Parameters
    ----------
    atoms_a, atoms_b : array_like
        Backbone indices of the paired monomers, integers in
        ``0, ..., n_backbone - 1``, in any two shapes that broadcast
        together.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed.
    n_backbone : int
        The number of monomers along the backbone; positive.

    Returns
    -------
    numpy.ndarray
        The genomic distance of each pair, the broadcast shape of the inputs,
        ``int64``, in bonds: ``|a - b|`` on a linear chain and
        ``min(|a - b|, n_backbone - |a - b|)`` on a ring, the shorter way
        round. The number of monomers between the pair is one less.

    Raises
    ------
    ValueError
        If `topology` is not ``'linear'`` or ``'ring'``, if `n_backbone` is
        not positive, or if an index is outside ``0, ..., n_backbone - 1``.

    Notes
    -----
    Time and space are :math:`O(m)` for :math:`m` pairs.

    Examples
    --------
    On a ring of 10 monomers, 0 and 9 are neighbours:

    >>> from polyphys.analyze.clusters import genomic_distance
    >>> genomic_distance([0, 1], [3, 9], 'linear', 10)
    array([3, 8])
    >>> genomic_distance([0, 1], [3, 9], 'ring', 10)
    array([3, 2])

    .. versionadded:: 0.6.0
    """
    _check_topology(topology)
    if n_backbone < 1:
        raise ValueError(f'n_backbone must be positive; got {n_backbone}.')
    a = np.asarray(atoms_a, dtype=np.int64)
    b = np.asarray(atoms_b, dtype=np.int64)
    for name, idx in (('atoms_a', a), ('atoms_b', b)):
        if idx.size and (idx.min() < 0 or idx.max() >= n_backbone):
            raise ValueError(
                f'{name} must hold indices in 0, ..., {n_backbone - 1}.'
            )
    separation = np.abs(a - b)
    if topology == 'ring':
        separation = np.minimum(separation, n_backbone - separation)
    return np.asarray(separation, dtype=np.int64)


def foci_pair_info(
    genomic_pos: npt.ArrayLike,
    n_monomers: int,
    topology: Topology,
) -> npt.NDArray[np.int64]:
    """
    List every pair of foci with its genomic distance.

    Parameters
    ----------
    genomic_pos : array_like
        Backbone indices of the foci, the tagged monomers, shape
        ``(n_foci,)``, integers in ``0, ..., n_monomers - 1``.
    n_monomers : int
        The number of monomers along the backbone, foci included.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_pairs, 3)``, ``int64``, one row per pair in the order of
        ``numpy.triu_indices(n_foci, 1)``: the two backbone indices and their
        genomic distance in bonds. ``n_pairs`` is ``n_foci (n_foci - 1) / 2``.

    Raises
    ------
    ValueError
        As :func:`genomic_distance`, or if `genomic_pos` is not
        one-dimensional.

    Notes
    -----
    The row order matches the upper triangle of a foci distance matrix, so
    ``dist[np.triu_indices(n_foci, 1)]`` gives the separations in the same
    order. Time and space are :math:`O(n_\\mathrm{foci}^2)`.

    Examples
    --------
    Three foci on a ring of 10 monomers:

    >>> from polyphys.analyze.clusters import foci_pair_info
    >>> foci_pair_info([0, 4, 8], 10, 'ring')
    array([[0, 4, 4],
           [0, 8, 2],
           [4, 8, 4]])

    .. versionadded:: 0.6.0
    """
    positions = np.asarray(genomic_pos, dtype=np.int64)
    if positions.ndim != 1:
        raise ValueError(
            'genomic_pos must be one-dimensional; got'
            f' {positions.ndim} dimension(s).'
        )
    first, second = np.triu_indices(positions.shape[0], 1)
    pos_a = positions[first]
    pos_b = positions[second]
    distance = genomic_distance(pos_a, pos_b, topology, n_monomers)
    return np.column_stack((pos_a, pos_b, distance)).astype(np.int64)


def pair_distance_histograms(
    distances: npt.ArrayLike,
    bin_size: float,
    r_max: float,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """
    Histogram the separation time series of each pair on one shared grid.

    Parameters
    ----------
    distances : array_like
        Pair separations over the frames, shape ``(n_pairs, n_frames)`` or
        ``(n_frames,)`` for one pair, in the caller's length unit;
        non-negative and finite.
    bin_size : float
        The bin width :math:`\\Delta r`, positive.
    r_max : float
        The upper end of the grid, positive. The grid runs from 0 to the
        first multiple of `bin_size` at or above `r_max`, so every pair,
        whole, and ensemble shares the same bins.

    Returns
    -------
    counts : numpy.ndarray
        Shape ``(n_pairs, n_bins)``, or ``(n_bins,)`` for a one-dimensional
        input, ``int64``: the number of frames whose separation falls in each
        bin. Each row sums to ``n_frames``.
    edges : numpy.ndarray
        Shape ``(n_bins + 1,)``, ``float64``: the bin edges,
        ``k * bin_size``. Bin ``k`` is ``[edges[k], edges[k + 1])``, except
        the last, which also includes its upper edge, as
        :func:`numpy.histogram` does.

    Raises
    ------
    ValueError
        If `bin_size` or `r_max` is not positive, if a separation is
        negative, NaN, or above the last edge, or if `distances` has more
        than two dimensions.

    Notes
    -----
    A separation above the grid raises rather than being dropped, because a
    dropped frame would silently break the normalization of
    :func:`pair_separation_density`. Choose `r_max` from the system, such as
    the largest separation the box allows, not from the data: a grid set by
    each whole's own maximum gives every whole a different number of bins,
    and an average over wholes is then taken over different sets of wholes
    in different bins.

    Time is :math:`O(n_\\mathrm{pairs} n_\\mathrm{frames} \\log
    n_\\mathrm{bins})`, from a binary search of the edges per frame; space is
    :math:`O(n_\\mathrm{pairs} (n_\\mathrm{frames} + n_\\mathrm{bins}))`.

    Examples
    --------
    >>> from polyphys.analyze.clusters import pair_distance_histograms
    >>> counts, edges = pair_distance_histograms(
    ...     [[0.5, 1.5, 1.5, 3.0], [2.0, 2.0, 2.0, 2.0]], bin_size=1.0,
    ...     r_max=3.0)
    >>> counts
    array([[1, 2, 1],
           [0, 0, 4]])
    >>> edges
    array([0., 1., 2., 3.])

    .. versionadded:: 0.6.0
    """
    if not bin_size > 0:
        raise ValueError(f'bin_size must be positive; got {bin_size}.')
    if not r_max > 0:
        raise ValueError(f'r_max must be positive; got {r_max}.')
    series = np.asarray(distances, dtype=np.float64)
    if series.ndim > 2:
        raise ValueError(
            'distances must have shape (n_pairs, n_frames) or (n_frames,);'
            f' got {series.ndim} dimensions.'
        )
    one_pair = series.ndim == 1
    rows = series.reshape(1, -1) if one_pair else series
    n_pairs, n_frames = rows.shape
    n_bins = int(np.ceil(r_max / bin_size))
    edges = bin_size * np.arange(n_bins + 1, dtype=np.float64)
    if np.any(np.isnan(rows)) or np.any(rows < 0) or np.any(rows > edges[-1]):
        raise ValueError(
            f'distances must lie in [0, {edges[-1]}], the histogram grid.'
        )
    # Bin k holds edges[k] <= r < edges[k + 1]; a separation on the last edge
    # goes into the last bin, as numpy.histogram does.
    bins = np.searchsorted(edges, rows, side='right') - 1
    bins = np.minimum(bins, n_bins - 1)
    flat = bins + n_bins * np.arange(n_pairs)[:, None]
    flat_counts = np.bincount(flat.ravel(), minlength=n_pairs * n_bins)
    counts = flat_counts.reshape(n_pairs, n_bins).astype(np.int64)
    if one_pair:
        return counts[0], edges
    return counts, edges


def pair_separation_density(
    counts: npt.ArrayLike,
    edges: npt.ArrayLike,
) -> npt.NDArray[np.float64]:
    """
    Estimate the probability density of a pair's separation vector.

    Parameters
    ----------
    counts : array_like
        Histogram counts from :func:`pair_distance_histograms`, shape
        ``(..., n_bins)``, one row per pair; every row must have a positive
        total.
    edges : array_like
        The bin edges, shape ``(n_bins + 1,)``, non-negative and strictly
        increasing, in the caller's length unit.

    Returns
    -------
    numpy.ndarray
        :math:`\\bar P(r)`, the same shape as `counts`, ``float64``, in units
        of one over the length unit cubed. For each row,
        ``sum(density * shell_volume)`` is 1.

    Raises
    ------
    ValueError
        If `edges` is not a strictly increasing non-negative vector of length
        ``n_bins + 1``, or if a row of `counts` is all zero.

    Notes
    -----
    Let :math:`\\vec R` be the separation vector of one pair of monomers and
    :math:`P(\\vec R)` its probability density over the chain's
    configurations, :math:`\\int P\\, d^3R = 1`. The scalar separation
    :math:`R = |\\vec R|` has density :math:`P(R) = 4\\pi R^2 \\bar P(R)`,
    where :math:`\\bar P(R)` is the average of :math:`P(\\vec R)` over the
    directions of :math:`\\vec R`; the factor :math:`4\\pi R^2` is the
    Jacobian of spherical coordinates, nothing else. With :math:`H_k` frames
    out of :math:`N_f` in the shell between :math:`R_k` and :math:`R_{k+1}`,
    whose volume is :math:`\\Delta V_k = \\tfrac{4\\pi}{3}(R_{k+1}^3 -
    R_k^3)`,

    .. math::

        \\bar P_k = \\frac{H_k}{N_f\\, \\Delta V_k},
        \\qquad \\sum_k \\bar P_k\\, \\Delta V_k = 1 .

    :math:`-k_B T \\ln \\bar P(R)` is the free energy of the pair at
    separation :math:`R`, up to a constant; for an ideal chain of :math:`n`
    Kuhn segments of length :math:`b`, :math:`P(\\vec R)` is the Gaussian
    :math:`(3 / 2\\pi n b^2)^{3/2} \\exp(-3R^2 / 2 n b^2)` (Rubinstein and
    Colby, 2003, Chapter 2), a reference the estimate can be compared with
    directly, since it carries no free constant. The estimate itself assumes
    no chain model and no geometry: the accessible volume of a confinement
    shapes :math:`\\bar P`, not its definition. The exact shell volume is
    used rather than :math:`4\\pi R_c^2 \\Delta r` at the bin centre, which
    is 25% too small in the first bin. Time and space are
    :math:`O(n_\\mathrm{pairs} n_\\mathrm{bins})`.

    Examples
    --------
    Four frames all in the shell between 2 and 3: the density is one over
    that shell's volume, :math:`\\tfrac{4\\pi}{3}(27 - 8)`.

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import pair_separation_density
    >>> density = pair_separation_density([0, 0, 4], [0.0, 1.0, 2.0, 3.0])
    >>> np.round(density, 5)
    array([0.     , 0.     , 0.01256])
    >>> shells = 4 * np.pi / 3 * np.diff(np.array([0.0, 1.0, 2.0, 3.0])**3)
    >>> float(np.sum(density * shells))
    1.0

    .. versionadded:: 0.6.0
    """
    hist = np.asarray(counts, dtype=np.float64)
    bounds = np.asarray(edges, dtype=np.float64)
    if bounds.ndim != 1 or bounds.shape[0] != hist.shape[-1] + 1:
        raise ValueError(
            'edges must have length n_bins + 1 ='
            f' {hist.shape[-1] + 1}; got shape {bounds.shape}.'
        )
    if bounds[0] < 0 or np.any(np.diff(bounds) <= 0):
        raise ValueError(
            'edges must be non-negative and strictly increasing.'
        )
    n_frames = hist.sum(axis=-1, keepdims=True)
    if np.any(n_frames <= 0):
        raise ValueError('every row of counts must have a positive total.')
    shell_volume = 4.0 * np.pi / 3.0 * np.diff(bounds**3)
    return np.asarray(hist / (n_frames * shell_volume), dtype=np.float64)


def enforce_single_patch_contact(
    direct_contacts: npt.ArrayLike,
) -> npt.NDArray[np.bool_]:
    """
    Keep at most one monomer per patch: the lowest-index one.

    Parameters
    ----------
    direct_contacts : array_like
        The monomer-patch contact matrix, shape ``(n_monomers, n_patches)``,
        boolean or 0/1: ``[i, j]`` true when monomer ``i`` touches patch
        ``j``.

    Returns
    -------
    numpy.ndarray
        The same shape, ``bool``, with every column reduced to its first true
        entry; columns without a contact stay empty.

    Raises
    ------
    ValueError
        If `direct_contacts` is not two-dimensional.

    Notes
    -----
    A patch is a single binding site and can hold one monomer, but a
    distance cutoff occasionally marks two or three adjacent monomers as its
    contacts. The project rule keeps the one with the lowest backbone index
    and treats the rest as artefacts. Applying the rule twice changes
    nothing. Time and space are :math:`O(n_\\mathrm{monomers}
    n_\\mathrm{patches})`.

    Examples
    --------
    Patch 0 touches monomers 1 and 2; only monomer 1 is kept.

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import enforce_single_patch_contact
    >>> contacts = np.array([[0, 1],
    ...                      [1, 0],
    ...                      [1, 0]])
    >>> enforce_single_patch_contact(contacts).astype(int)
    array([[0, 1],
           [1, 0],
           [0, 0]])

    .. versionadded:: 0.6.0
    """
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_monomers = contacts.shape[0]
    first = np.argmax(contacts, axis=0)
    engaged = contacts.any(axis=0)
    rows = np.arange(n_monomers)[:, None]
    return np.asarray((rows == first[None, :]) & engaged[None, :],
                      dtype=np.bool_)


def monomer_binder_contacts(
    direct_contacts: npt.ArrayLike,
    patches_per_binder: int,
) -> npt.NDArray[np.bool_]:
    """
    Collapse the patches of each binder into one column.

    Parameters
    ----------
    direct_contacts : array_like
        The monomer-patch contact matrix, shape ``(n_monomers, n_patches)``,
        boolean or 0/1, with the patches of binder ``k`` in columns
        ``k * patches_per_binder`` to ``(k + 1) * patches_per_binder - 1``.
    patches_per_binder : int
        The number of patches per binder; positive and dividing
        ``n_patches``.

    Returns
    -------
    numpy.ndarray
        The monomer-binder contact matrix, shape ``(n_monomers, n_binders)``,
        ``bool``: ``[i, k]`` true when monomer ``i`` touches any patch of
        binder ``k``.

    Raises
    ------
    ValueError
        If `direct_contacts` is not two-dimensional, or if
        `patches_per_binder` is not positive or does not divide the number of
        columns.

    Notes
    -----
    Time and space are :math:`O(n_\\mathrm{monomers} n_\\mathrm{patches})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.clusters import monomer_binder_contacts
    >>> contacts = np.array([[1, 0, 0, 0],
    ...                      [0, 1, 0, 0],
    ...                      [0, 0, 0, 1]])
    >>> monomer_binder_contacts(contacts, 2).astype(int)
    array([[1, 0],
           [1, 0],
           [0, 1]])

    .. versionadded:: 0.6.0
    """
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_monomers, n_patches = contacts.shape
    if patches_per_binder < 1:
        raise ValueError(
            f'patches_per_binder must be positive; got {patches_per_binder}.'
        )
    if n_patches % patches_per_binder != 0:
        raise ValueError(
            f'the number of columns, {n_patches}, must be a multiple of'
            f' patches_per_binder, {patches_per_binder}.'
        )
    n_binders = n_patches // patches_per_binder
    grouped = contacts.reshape(n_monomers, n_binders, patches_per_binder)
    return np.asarray(grouped.any(axis=2), dtype=np.bool_)


def split_binders(
    direct_contacts: npt.ArrayLike,
) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.bool_]]:
    """
    Separate the dangling binders from the bridging ones.

    Parameters
    ----------
    direct_contacts : array_like
        The monomer-binder contact matrix, shape ``(n_monomers, n_binders)``,
        boolean or 0/1, as :func:`monomer_binder_contacts` returns it.

    Returns
    -------
    dangling : numpy.ndarray
        The same shape, ``bool``: the columns that touch exactly one monomer,
        every other column cleared.
    bridging : numpy.ndarray
        The same shape, ``bool``: the columns that touch exactly two
        monomers, every other column cleared.

    Raises
    ------
    ValueError
        If `direct_contacts` is not two-dimensional.

    Notes
    -----
    A binder touching no monomer is free and appears in neither output; one
    touching three or more monomers, which a two-patch binder cannot do once
    :func:`enforce_single_patch_contact` has been applied, appears in
    neither either. Time and space are :math:`O(n_\\mathrm{monomers}
    n_\\mathrm{binders})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.clusters import split_binders
    >>> contacts = np.array([[1, 1, 0],
    ...                      [0, 1, 0],
    ...                      [0, 0, 0]])
    >>> dangling, bridging = split_binders(contacts)
    >>> dangling.astype(int)
    array([[1, 0, 0],
           [0, 0, 0],
           [0, 0, 0]])
    >>> bridging.astype(int)
    array([[0, 1, 0],
           [0, 1, 0],
           [0, 0, 0]])

    .. versionadded:: 0.6.0
    """
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_monomers_bound = contacts.sum(axis=0)
    dangling = contacts & (n_monomers_bound == 1)[None, :]
    bridging = contacts & (n_monomers_bound == 2)[None, :]
    return (np.asarray(dangling, dtype=np.bool_),
            np.asarray(bridging, dtype=np.bool_))


def binder_contacts(
    direct_contacts: npt.ArrayLike,
    topology: Topology,
) -> npt.NDArray[np.bool_]:
    """
    Link binders that touch the same or neighbouring monomers.

    Parameters
    ----------
    direct_contacts : array_like
        The monomer-binder contact matrix, shape ``(n_monomers, n_binders)``,
        boolean or 0/1, with the monomers in backbone order along the rows.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed; on a ring the first and last
        monomers are neighbours.

    Returns
    -------
    numpy.ndarray
        The binder-binder direct-contact matrix, shape
        ``(n_binders, n_binders)``, ``bool``: ``[k, l]`` true when binders
        ``k`` and ``l`` touch the same monomer or two monomers adjacent along
        the backbone. Symmetric, with a true diagonal for every binder, bound
        or not, so :func:`count_clusters` sees every binder as a cluster of
        at least size 1.

    Raises
    ------
    ValueError
        If `direct_contacts` is not two-dimensional or `topology` is invalid.

    Notes
    -----
    With :math:`M` the monomer-binder matrix and :math:`N_i = M_i \\lor
    M_{i+1}` the matrix of the monomer pairs ``(i, i + 1)`` that are
    neighbours, the result is :math:`M^\\top M > 0 \\lor N^\\top N > 0`
    plus the identity: the first term links binders sharing a monomer, the
    second links binders on neighbouring monomers. Time is
    :math:`O(n_\\mathrm{monomers} n_\\mathrm{binders}^2)` from the two
    products; space is :math:`O(n_\\mathrm{monomers} n_\\mathrm{binders} +
    n_\\mathrm{binders}^2)`.

    Examples
    --------
    Binder 0 on monomer 0 and binder 1 on monomer 3 are linked on a ring of
    four monomers, where 0 and 3 are neighbours, but not on a linear chain:

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import binder_contacts
    >>> contacts = np.array([[1, 0],
    ...                      [0, 0],
    ...                      [0, 0],
    ...                      [0, 1]])
    >>> binder_contacts(contacts, 'ring').astype(int)
    array([[1, 1],
           [1, 1]])
    >>> binder_contacts(contacts, 'linear').astype(int)
    array([[1, 0],
           [0, 1]])

    .. versionadded:: 0.6.0
    """
    _check_topology(topology)
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_binders = contacts.shape[1]
    if topology == 'ring':
        neighbours = contacts | np.roll(contacts, -1, axis=0)
    else:
        neighbours = contacts[:-1] | contacts[1:]
    weights = contacts.astype(np.int64)
    shared = weights.T @ weights > 0
    adjacent_weights = neighbours.astype(np.int64)
    adjacent = adjacent_weights.T @ adjacent_weights > 0
    linked = shared | adjacent | np.eye(n_binders, dtype=np.bool_)
    return np.asarray(linked, dtype=np.bool_)


@dataclass(frozen=True, eq=False)
class FrameBinding:
    """
    The binding statistics of two-patch binders in one frame.

    Attributes
    ----------
    n_patch_monomer_contacts : int
        Direct contacts between patches and monomers, before the single-patch
        rule; a patch touching two monomers contributes two.
    n_patches_engaged : int
        Patches touching at least one monomer.
    n_patches_free : int
        Patches touching none; ``n_patches_engaged + n_patches_free`` is the
        number of patches.
    n_binders_free : int
        Binders with both patches free.
    n_binders_dangling : int
        Binders with exactly one patch engaged.
    n_binders_bridging : int
        Binders with both patches engaged. The three binder counts sum to the
        number of binders.
    n_pairs_cis : int
        Distinct monomer pairs bridged at a genomic distance of at least 1 and
        at most ``cis_threshold`` bonds.
    n_pairs_trans : int
        Distinct monomer pairs bridged at more than ``cis_threshold`` bonds.
    bridged_pairs : numpy.ndarray
        Shape ``(n_bridged_pairs, 3)``, ``int64``: for each distinct bridged
        pair, the two monomer indices, lower first, and their genomic
        distance in bonds. ``n_pairs_cis + n_pairs_trans`` is its number of
        rows.

    Notes
    -----
    The pair counts are per distinct monomer pair, so two binders bridging
    the same pair count once, and a binder whose two patches hold the same
    monomer is bridging but bridges no pair.

    Examples
    --------
    >>> from polyphys.analyze.clusters import binding_statistics
    >>> frame = binding_statistics([[1, 0], [0, 1]], 'linear',
    ...                            cis_threshold=1)
    >>> frame.n_binders_bridging, frame.n_pairs_cis
    (1, 1)

    .. versionadded:: 0.6.0
    """

    n_patch_monomer_contacts: int
    n_patches_engaged: int
    n_patches_free: int
    n_binders_free: int
    n_binders_dangling: int
    n_binders_bridging: int
    n_pairs_cis: int
    n_pairs_trans: int
    bridged_pairs: npt.NDArray[np.int64]


def binding_statistics(
    direct_contacts: npt.ArrayLike,
    topology: Topology,
    cis_threshold: int,
) -> FrameBinding:
    """
    Count free, dangling, bridging, cis, and trans binders in one frame.

    Parameters
    ----------
    direct_contacts : array_like
        The monomer-patch contact matrix, shape ``(n_monomers, n_patches)``,
        boolean or 0/1, with the monomers in backbone order along the rows
        and the two patches of binder ``k`` in columns ``2k`` and ``2k + 1``.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed.
    cis_threshold : int
        The largest genomic distance, in bonds, of a cis bridge;
        non-negative. A project constant.

    Returns
    -------
    FrameBinding
        The counts and the bridged pairs; see its attributes.

    Raises
    ------
    ValueError
        If `direct_contacts` is not two-dimensional or has an odd number of
        columns, if `topology` is invalid, or if `cis_threshold` is negative.

    Notes
    -----
    The patch and binder counts use the contacts as given, so a patch
    marked on two monomers still counts as one engaged patch and its binder
    as one dangling or bridging binder. The bridged pairs use the contacts
    after :func:`enforce_single_patch_contact`, so each bridging binder
    joins exactly two monomers; the pairs are then made distinct, their
    genomic distance is taken from :func:`genomic_distance`, and a pair of
    distance 0 (both patches on one monomer) is dropped.

    Time is :math:`O(n_\\mathrm{monomers} n_\\mathrm{patches} +
    n_\\mathrm{binders} \\log n_\\mathrm{binders})`; space is
    :math:`O(n_\\mathrm{monomers} n_\\mathrm{patches})`.

    Examples
    --------
    A ring of six monomers and four binders. Binder 0 bridges monomers 0 and
    2; binder 1 has patch 2 on monomers 3 and 4 (an artefact) and patch 3
    free, so it dangles; binder 2 is free; binder 3 bridges monomers 5 and
    0, neighbours on the ring. With ``cis_threshold=1`` the second bridge is
    cis and the first, at 2 bonds, is trans.

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import binding_statistics
    >>> contacts = np.zeros((6, 8), dtype=int)
    >>> contacts[0, 0] = contacts[2, 1] = 1
    >>> contacts[3, 2] = contacts[4, 2] = 1
    >>> contacts[5, 6] = contacts[0, 7] = 1
    >>> frame = binding_statistics(contacts, 'ring', cis_threshold=1)
    >>> (frame.n_patch_monomer_contacts, frame.n_patches_engaged,
    ...  frame.n_patches_free)
    (6, 5, 3)
    >>> (frame.n_binders_free, frame.n_binders_dangling,
    ...  frame.n_binders_bridging)
    (1, 1, 2)
    >>> frame.n_pairs_cis, frame.n_pairs_trans
    (1, 1)
    >>> frame.bridged_pairs
    array([[0, 2, 2],
           [0, 5, 1]])

    .. versionadded:: 0.6.0
    """
    _check_topology(topology)
    if cis_threshold < 0:
        raise ValueError(
            f'cis_threshold must be non-negative; got {cis_threshold}.'
        )
    contacts = _as_contacts(direct_contacts, 'direct_contacts')
    n_monomers, n_patches = contacts.shape
    if n_patches % 2 != 0:
        raise ValueError(
            'direct_contacts must have two columns per binder; got'
            f' {n_patches} columns.'
        )
    n_binders = n_patches // 2
    per_patch = contacts.sum(axis=0)
    engaged = per_patch > 0
    n_patches_engaged = int(np.count_nonzero(engaged))
    first_engaged = engaged[0::2]
    second_engaged = engaged[1::2]
    n_binders_free = int(np.count_nonzero(~first_engaged & ~second_engaged))
    n_binders_bridging = int(
        np.count_nonzero(first_engaged & second_engaged)
    )
    # Bridged pairs: one monomer per patch, then the two patches of each
    # bridging binder.
    single = enforce_single_patch_contact(contacts)
    monomer_of_patch = np.argmax(single, axis=0)
    bridging = first_engaged & second_engaged
    pairs: npt.NDArray[np.int64] = np.column_stack(
        (monomer_of_patch[0::2][bridging], monomer_of_patch[1::2][bridging])
    )
    pairs = np.sort(pairs, axis=1)
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    pairs = np.unique(pairs, axis=0).reshape(-1, 2)
    distance = genomic_distance(
        pairs[:, 0], pairs[:, 1], topology, n_monomers
    )
    bridged_pairs = np.column_stack((pairs, distance)).astype(np.int64)
    n_pairs_cis = int(np.count_nonzero(distance <= cis_threshold))
    return FrameBinding(
        n_patch_monomer_contacts=int(per_patch.sum()),
        n_patches_engaged=n_patches_engaged,
        n_patches_free=n_patches - n_patches_engaged,
        n_binders_free=n_binders_free,
        n_binders_dangling=n_binders - n_binders_free - n_binders_bridging,
        n_binders_bridging=n_binders_bridging,
        n_pairs_cis=n_pairs_cis,
        n_pairs_trans=int(distance.shape[0]) - n_pairs_cis,
        bridged_pairs=bridged_pairs,
    )


@dataclass
class BindingStats:
    """
    Binding statistics accumulated over frames, with the loop-length
    histogram.

    Parameters
    ----------
    n_monomers : int
        The number of monomers along the backbone; positive.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed.
    cis_threshold : int
        The largest genomic distance, in bonds, of a cis bridge;
        non-negative.

    Attributes
    ----------
    frames : list of FrameBinding
        The frames appended so far, in order.
    loop_length_histogram : numpy.ndarray
        ``int64``: ``[d]`` is the number of distinct bridged pairs over all
        frames at a genomic distance of ``d`` bonds. Its length is
        ``n_monomers`` on a linear chain, where the largest distance is
        ``n_monomers - 1``, and ``n_monomers // 2 + 1`` on a ring, where it
        is ``n_monomers // 2``. ``[0]`` is always 0.

    Notes
    -----
    Building the object costs :math:`O(n_\\mathrm{monomers})` space for the
    histogram; each :meth:`append` costs :math:`O(n_\\mathrm{bridged\\
    pairs})`.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.clusters import BindingStats
    >>> stats = BindingStats(n_monomers=4, topology='linear',
    ...                      cis_threshold=1)
    >>> stats.update([[1, 0], [0, 1], [0, 0], [0, 0]])
    >>> stats.update([[1, 0], [0, 0], [0, 0], [0, 1]])
    >>> stats.n_frames
    2
    >>> stats.as_dict()['n_pairs_cis']
    array([1, 0])
    >>> stats.loop_length_histogram
    array([0, 1, 0, 1])

    .. versionadded:: 0.6.0
    """

    n_monomers: int
    topology: Topology
    cis_threshold: int
    frames: list[FrameBinding] = field(default_factory=list)
    loop_length_histogram: npt.NDArray[np.int64] = field(init=False)

    def __post_init__(self) -> None:
        _check_topology(self.topology)
        if self.n_monomers < 1:
            raise ValueError(
                f'n_monomers must be positive; got {self.n_monomers}.'
            )
        if self.cis_threshold < 0:
            raise ValueError(
                f'cis_threshold must be non-negative; got'
                f' {self.cis_threshold}.'
            )
        if self.topology == 'ring':
            n_lengths = self.n_monomers // 2 + 1
        else:
            n_lengths = self.n_monomers
        self.loop_length_histogram = np.zeros(n_lengths, dtype=np.int64)

    @property
    def n_frames(self) -> int:
        """The number of frames appended so far."""
        return len(self.frames)

    def append(self, frame: FrameBinding) -> None:
        """
        Add one frame's statistics and its bridged pairs to the histogram.

        Parameters
        ----------
        frame : FrameBinding
            As :func:`binding_statistics` returns it, computed with this
            object's `topology` and `cis_threshold` for a chain of
            `n_monomers`.

        Raises
        ------
        ValueError
            If a bridged pair's genomic distance does not fit the histogram,
            which means the frame came from a longer chain or another
            topology.
        """
        distance = frame.bridged_pairs[:, 2]
        if distance.size and distance.max() >= self.loop_length_histogram.size:
            raise ValueError(
                'a bridged pair has a genomic distance beyond the chain;'
                ' the frame does not match n_monomers and topology.'
            )
        np.add.at(self.loop_length_histogram, distance, 1)
        self.frames.append(frame)

    def update(self, direct_contacts: npt.ArrayLike) -> None:
        """
        Compute one frame's statistics and append them.

        Parameters
        ----------
        direct_contacts : array_like
            The monomer-patch contact matrix of the frame, shape
            ``(n_monomers, n_patches)``, as :func:`binding_statistics` takes
            it.

        Raises
        ------
        ValueError
            As :func:`binding_statistics`, or if the matrix has not
            `n_monomers` rows.
        """
        contacts = _as_contacts(direct_contacts, 'direct_contacts')
        if contacts.shape[0] != self.n_monomers:
            raise ValueError(
                f'direct_contacts must have n_monomers = {self.n_monomers}'
                f' rows; got {contacts.shape[0]}.'
            )
        self.append(
            binding_statistics(contacts, self.topology, self.cis_threshold)
        )

    def as_dict(self) -> dict[str, npt.NDArray[np.int64]]:
        """
        Return the eight counts as time series.

        Returns
        -------
        dict of numpy.ndarray
            One ``int64`` array of shape ``(n_frames,)`` per count, keyed by
            the :class:`FrameBinding` attribute name; `bridged_pairs` is not
            included.
        """
        return {
            name: np.array(
                [getattr(frame, name) for frame in self.frames],
                dtype=np.int64,
            )
            for name in _BINDING_COUNTS
        }


def threshold_sweep(
    distances: npt.ArrayLike,
    cutoffs: Sequence[float] | npt.ArrayLike,
    cis_thresholds: Sequence[int],
    topology: Topology,
    *,
    inclusive: bool,
) -> tuple[pd.DataFrame, npt.NDArray[np.float64]]:
    """
    Recount the binding statistics of saved distance matrices at every
    cutoff and cis threshold of a grid.

    Parameters
    ----------
    distances : array_like
        The monomer-to-patch distance matrices of one whole, shape
        ``(n_frames, n_monomers, n_patches)``, in the caller's length unit,
        with the monomers in backbone order along the second axis and the
        two patches of binder ``k`` at ``2k`` and ``2k + 1`` along the
        third; the ``distances`` that
        :class:`~polyphys.probe.binders.BinderStatistics` stores.
    cutoffs : sequence of float
        The contact distances :math:`r_c` to sweep, in the same unit, each
        positive; non-empty. Swept in the order given.
    cis_thresholds : sequence of int
        The largest genomic distances, in bonds, of a cis bridge to sweep;
        each non-negative, non-empty. Swept in the order given, inside each
        cutoff.
    topology : {'linear', 'ring'}
        Whether the chain is open or closed.
    inclusive : bool
        Count a pair at exactly the cutoff as a contact (``dist <= cutoff``)
        or not (``dist < cutoff``), as :func:`find_direct_contacts`. Keyword
        only, with no default, because the comparison decides whether the
        sweep reproduces stored counts: a sweep agrees with the series that
        :class:`~polyphys.probe.binders.BinderStatistics` saved only at the
        spec's cutoff, cis threshold, and comparison.

    Returns
    -------
    counts : pandas.DataFrame
        One row per ``(cutoff, cis_threshold)`` pair, cutoffs outer and cis
        thresholds inner, in the order given, with no row repeated. Columns:
        ``cutoff`` (float64), ``cis_threshold`` (int64), and the eight
        counts of :class:`FrameBinding`, each the mean over the frames
        (float64), under the attribute names ``n_patch_monomer_contacts``,
        ``n_patches_engaged``, ``n_patches_free``, ``n_binders_free``,
        ``n_binders_dangling``, ``n_binders_bridging``, ``n_pairs_cis``, and
        ``n_pairs_trans``.
    loop_length_histogram : numpy.ndarray
        Shape ``(n_rows, n_lengths)``, ``float64``, one row per row of
        `counts`: ``[k, d]`` is the mean over the frames of the number of
        distinct bridged pairs at a genomic distance of ``d`` bonds at the
        cutoff of row ``k``. ``n_lengths`` is ``n_monomers`` on a linear
        chain and ``n_monomers // 2 + 1`` on a ring, as
        :class:`BindingStats`; ``[:, 0]`` is always 0. The cis threshold
        does not enter the histogram, so the rows of one cutoff are equal.

    Raises
    ------
    ValueError
        If `distances` is not three-dimensional or has no frames, if
        `cutoffs` or `cis_thresholds` is empty, if a cutoff is not positive
        or a cis threshold is negative, if the number of patches is odd, or
        if `topology` is invalid.

    Notes
    -----
    This is the research code's choice-of-threshold sweep on the saved
    ``distMatTMonHnsPatch`` matrices, with its two tables, the
    frame-averaged counts and the frame-averaged loop-length histogram,
    returned side by side rather than written with repeated rows. The
    research sweep compared strictly, while the H-NS study specs count with
    ``inclusive=True`` at 0.685; the strict sweep's row at 0.685 therefore
    omits the pairs at exactly the cutoff and does not match the stored
    series. The research grid, ``numpy.arange(0.535, 0.685 + 0.01, 0.01)``,
    also ran one step past its end to 0.695; build a grid with a count,
    ``0.535 + 0.01 * numpy.arange(16)``, to stop at 0.685.

    Every count and the histogram are means over the frames, which treats
    the frames as exchangeable; they are not uncertainties. An error bar on
    a swept count needs the per-frame series of :class:`BindingStats` and
    an autocorrelation-aware estimate such as
    :func:`polyphys.analyze.stats.block_average`.

    The contacts and the bridged pairs of a frame depend on the cutoff only;
    the cis threshold splits the bridged pairs into cis and trans. Each
    frame is therefore counted once per cutoff, with
    :func:`find_direct_contacts` and :func:`binding_statistics`, and its
    bridged pairs are split once per cis threshold. With :math:`n_c`
    cutoffs, :math:`n_t` cis thresholds, :math:`n_f` frames,
    :math:`N_m` monomers, :math:`N_p` patches, and :math:`N_b = N_p / 2`
    binders, time is :math:`O(n_c n_f (N_m N_p + N_b \\log N_b) + n_c n_t
    n_f N_b)` and space, beyond the input, is :math:`O(N_m N_p + n_c n_t
    n_\\mathrm{lengths})`.

    Examples
    --------
    A linear chain of three monomers and one binder, in two frames. In the
    first frame, patch 0 sits 0.5 from monomer 0 and patch 1 sits 0.6 from
    monomer 2; in the second, patch 1 is far from every monomer. Below 0.6
    the binder dangles in both frames; at 0.65 it bridges monomers 0 and 2,
    two bonds apart, in the first frame, a trans bridge for a cis threshold
    of 1.

    >>> import numpy as np
    >>> from polyphys.analyze.clusters import threshold_sweep
    >>> distances = np.full((2, 3, 2), 2.0)
    >>> distances[0, 0, 0] = distances[1, 0, 0] = 0.5
    >>> distances[0, 2, 1] = 0.6
    >>> counts, loops = threshold_sweep(
    ...     distances, [0.55, 0.65], [1], 'linear', inclusive=True)
    >>> counts[['cutoff', 'n_binders_dangling', 'n_binders_bridging',
    ...         'n_pairs_trans']]
       cutoff  n_binders_dangling  n_binders_bridging  n_pairs_trans
    0    0.55                 1.0                 0.0            0.0
    1    0.65                 0.5                 0.5            0.5
    >>> loops
    array([[0. , 0. , 0. ],
           [0. , 0. , 0.5]])

    .. versionadded:: 0.7.0
    """
    _check_topology(topology)
    matrices = np.asarray(distances, dtype=np.float64)
    if matrices.ndim != 3:
        raise ValueError(
            'distances must have shape (n_frames, n_monomers, n_patches);'
            f' got {matrices.ndim} dimension(s).'
        )
    n_frames, n_monomers, n_patches = matrices.shape
    if n_frames == 0:
        raise ValueError('distances must hold at least one frame.')
    if n_patches % 2 != 0:
        raise ValueError(
            'distances must have two patches per binder; got'
            f' {n_patches} patches.'
        )
    cutoff_grid = np.asarray(cutoffs, dtype=np.float64).reshape(-1)
    if cutoff_grid.size == 0:
        raise ValueError('cutoffs must not be empty.')
    if not np.all(cutoff_grid > 0):
        raise ValueError(
            f'every cutoff must be positive; got {cutoff_grid.tolist()}.'
        )
    cis_grid = np.asarray(cis_thresholds, dtype=np.int64).reshape(-1)
    if cis_grid.size == 0:
        raise ValueError('cis_thresholds must not be empty.')
    if not np.all(cis_grid >= 0):
        raise ValueError(
            'every cis threshold must be non-negative; got'
            f' {cis_grid.tolist()}.'
        )
    if topology == 'ring':
        n_lengths = n_monomers // 2 + 1
    else:
        n_lengths = n_monomers
    n_rows = cutoff_grid.size * cis_grid.size
    means = np.zeros((n_rows, len(_BINDING_COUNTS)), dtype=np.float64)
    histogram = np.zeros((n_rows, n_lengths), dtype=np.float64)
    # The six counts that do not depend on the cis threshold, then the
    # cis and trans pair counts, in the order of _BINDING_COUNTS.
    n_fixed = len(_BINDING_COUNTS) - 2
    for i, cutoff in enumerate(cutoff_grid):
        fixed = np.zeros(n_fixed, dtype=np.int64)
        n_cis = np.zeros(cis_grid.size, dtype=np.int64)
        n_pairs = 0
        loops = np.zeros(n_lengths, dtype=np.int64)
        for matrix in matrices:
            contacts = find_direct_contacts(matrix, cutoff, inclusive)
            # Any cis threshold gives the same contacts and bridged pairs;
            # 0 makes every pair trans, and the split is redone below.
            frame = binding_statistics(contacts, topology, 0)
            fixed += [getattr(frame, name)
                      for name in _BINDING_COUNTS[:n_fixed]]
            distance = frame.bridged_pairs[:, 2]
            n_pairs += distance.size
            n_cis += np.count_nonzero(
                distance[None, :] <= cis_grid[:, None], axis=1)
            loops += np.bincount(distance, minlength=n_lengths)
        rows = slice(i * cis_grid.size, (i + 1) * cis_grid.size)
        means[rows, :n_fixed] = fixed / n_frames
        means[rows, n_fixed] = n_cis / n_frames
        means[rows, n_fixed + 1] = (n_pairs - n_cis) / n_frames
        histogram[rows] = loops / n_frames
    counts = pd.DataFrame(means, columns=list(_BINDING_COUNTS))
    counts.insert(0, 'cis_threshold', np.tile(cis_grid, cutoff_grid.size))
    counts.insert(0, 'cutoff', np.repeat(cutoff_grid, cis_grid.size))
    return counts, histogram
