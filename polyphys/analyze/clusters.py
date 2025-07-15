"""
==========================================================
Cluster Analysis --- :mod:`polyphys.analyze.clusters`
==========================================================

The :mod:`polyphys.analyze.clusters` module provides a collection of functions
for analyzing and detecting particle clusters in polymer systems. These
functions are tailored to specific analysis needs, particularly for H-NS
proteins and other binding interactions in polymer simulations.

Functions
=========

.. autofunction:: count_foci_bonds
.. autofunction:: count_foci_clusters
.. autofunction:: foci_info
.. autofunction:: foci_rdf
.. autofunction:: foci_histogram
.. autofunction:: whole_dist_mat_foci
.. autofunction:: genomic_distance
.. autofunction:: hns_binding
.. autofunction:: enforce_single_patch_dir_contact
.. autofunction:: generate_mon_bind_direct
.. autofunction:: split_binder_matrix
.. autofunction:: find_binder_clusters_need_attention
.. autofunction:: find_binder_clusters_need_attention_old

Notes
=====
Many of the functions in this submodule are highly specific to certain analysis
tasks within the broader polymer physics simulations.

References
==========
- Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
  calculations of cluster statistics in continuum models of composite
  morphology. The Journal of Chemical Physics, 89(1), 668-676.
  https://doi.org/10.1063/1.454720
"""
from typing import Dict, List, Tuple, Callable
from itertools import combinations
import numpy as np
import pandas as pd
from ..manage.types import TopologyT, ParserType
from ..manage.utils import invalid_keyword
from .contacts import generate_contact_matrix


def count_foci_bonds(dir_contacts: np.ndarray) -> np.ndarray:
    """
    Compute the number and frequency of bonds in a direct contact matrix.

    This function analyzes the direct contact matrix to infer the distribution
    of bonds per atom. Each atom's number of bonds is determined by summing
    the contact matrix along its rows or columns, subtracting the diagonal
    self-bonds. Thus, number of bonds ranges from 0 to :math:`n_{atoms}-1`

    Parameters
    ----------
    dir_contacts : np.ndarray, shape (n_atoms, n_atoms)
        A symmetric binary matrix where element `(i, j)` is 1 if a contact
        exists between atoms `i` and `j`, and 0 otherwise. Diagonal elements
        must be 1, indicating self-bonds.

    Returns
    -------
    bind_list: np.ndarray, shape (n_atoms-1,)
        A 1D array where the index represents bond size (number of bonds per
        atom), and the value is the frequency of that bond size in the system.

    Raises
    ------
    ValueError
         If the input `dir_contacts` matrix is not square.

    Examples
    --------
    >>> dir_contacts = np.array([[1, 1, 0],
    ...                          [1, 1, 1],
    ...                          [0, 1, 1]])
    >>> count_foci_bonds(dir_contacts)
    array([0, 1, 2])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720
    """
    # Ensure the input is a square matrix
    n_atoms, n_cols = dir_contacts.shape

    if n_atoms != n_cols:
        raise ValueError("The direct contact matrix must be square.")

    # Sum along rows to count bonds for each atom and exclude self-bonds
    bonds = dir_contacts.sum(axis=1) - 1  # Removing diagonal self-bonds

    # Create a histogram of bond sizes
    bond_sizes, bond_counts = np.unique(bonds, return_counts=True)
    bond_list = np.zeros(n_atoms, dtype=int)  # Array for bond frequencies
    np.put(bond_list, bond_sizes, bond_counts)
    return bond_list


def count_foci_clusters(contacts: np.ndarray) -> np.ndarray:
    """
    Infer the number and size of clusters from a contact matrix.

    This function calculates the eigenvalues of the provided contact matrix.
    The number of non-zero eigenvalues corresponds to the number of clusters,
    and their values indicate the sizes of these clusters.

    Before computing eigenvalues, the function checks if the contact matrix is
    symmetric and positive semi-definite. Eigenvalues with negligible imaginary
    parts are treated as real.

    Parameters
    ----------
    contacts : np.ndarray
        A binary square matrix representing contacts between atoms. Each
        element is either 0 or 1, indicating the absence or presence of a
        contact, respectively. Diagonal elements are assumed to be 1,

    Returns
    -------
    cluster_list: np.ndarray, shape (cluster_size-1,)
        A 1D array where each index represents :math:`cluster_size - 1`, and
        the values are the frequency of clusters of that size.

    Raises
    ------
    ValueError
        If the input matrix is not symmetric, not positive semi-definite,
        or contains invalid cluster sizes.

    Notes
    -----
    An eigenvalue of 0 does not have a direct physical meaning; however, the
    frequency of this eigenvalue corresponds to the total number of repeated
    clusters of size 2 or larger. In terms of linear algebra, this is
    equivalent to the number of dependent columns in the contact matrix.

    Examples
    --------
    >>> contacts = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])
    >>> count_foci_clusters(contacts)
    array([1, 1, 1])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720
    """
    n_atoms, _ = contacts.shape
    # Compute eigenvalues and round for stability
    clusters_row = np.linalg.eigvalsh(contacts)
    clusters = np.asarray(np.round(clusters_row), dtype=int)

    # Validate eigenvalues for physical consistency
    if np.any(clusters < 0):
        raise ValueError("A cluster with size smaller than 0 found!")

    # Count occurrences of each cluster size
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    cluster_list = np.zeros(n_atoms + 1, dtype=int)
    np.put(cluster_list, cluster_sizes, cluster_counts)

    return cluster_list


def foci_info(
    genomic_pos: np.ndarray,
    nmon_large: int,
    nmon_small: int
) -> np.ndarray:
    """
    Generate information about all possible pairs of large monomers in a
    circular polymer.

    This function calculates the genomic indices and distances for all pairs of
    large monomers in a circular heterogeneous polymer composed of two types of
    monomers: small and large. The genomic distance is defined as the minimum
    distance along the ring topology.

    Parameters
    ----------
    genomic_pos : np.ndarray
        A 1D array of size `nmon_large` containing the genomic positions of
        large monomers.
    nmon_large : int
        Number of large monomers.
    nmon_small : int
        Number of small monomers.

    Returns
    -------
    pair_info: np.ndarray
        An array of shape `(n_pairs, 3)` where:
        - The first two columns represent the indices of monomer pairs along
          the backbone.
        - The third column contains the genomic distance between the pair.

    Notes
    -----
    The genomic distance accounts for the circular topology of the polymer. It
    is the minimum of the direct distance and the distance around the circle.

    Examples
    --------
    >>> genomic_pos = np.array([0, 5, 10, 15])
    >>> foci_info(genomic_pos, nmon_large=4, nmon_small=6)
    array([[ 0,  5,  4],
           [ 0, 10,  3],
           [ 0, 15,  2],
           [ 5, 10,  4],
           [ 5, 15,  3],
           [10, 15,  4]])
    """
    nmon = nmon_small + nmon_large

    # Generate all combinations of large monomer pairs
    pairs = np.array(list(combinations(range(nmon_large), 2)), dtype=np.int_)

    # Retrieve genomic positions for each pair
    # genomic_idx = np.choose(pairs, genomic_pos)
    # to-do: check wither choose is better or take.
    genomic_idx = np.take(genomic_pos, pairs)

    # Calculate direct genomic distances
    genomic_dist = np.abs(np.diff(genomic_idx, axis=1))
    # Applying the topological contain on the genomic distance, i.e. finding
    # minimum distance between a pair on a circle. The pair genomic distance is
    # inclusive, so their difference is one unit more than the actual number of
    # genomic units between the pair, but it is equal to the number of bonds
    # (not the total length of the bonds).

    # Adjust distances for circular topology
    # Number of small monomers between to large monomers:
    genomic_dist = np.minimum.reduce((genomic_dist, nmon - genomic_dist)) - 1

    # Combine pair indices and genomic distances
    # pair_info = np.concatenate((genomic_idx, genomic_dist), axis=1)
    pair_info = np.hstack((genomic_idx, genomic_dist))

    return pair_info


def foci_rdf(
    hist: np.ndarray,
    bin_centers: np.ndarray,
    binsize: float
) -> np.ndarray:
    """
    Compute the radial distribution function (RDF).

    The RDF provides the probability of finding a particle at a given radial
    distance relative to another particle. It is calculated by normalizing
    particle frequencies in each radial bin with respect to the shell volume.

    Parameters
    ----------
    hist : np.ndarray
        An array containing the frequencies of particles in each radial bin.
    bin_centers : np.ndarray
        An array containing the positions of bin centers along the radial
        direction.
    binsize : float
        The width of each radial bin.

    Returns
    -------
    np.ndarray
        The computed RDF values for each radial bin.

    Notes
    -----
    The RDF is normalized such that the integral over all bins equals 1.

    Examples
    --------
    >>> hist = np.array([10, 20, 30])
    >>> bin_centers = np.array([1, 2, 3])
    >>> binsize = 0.5
    >>> foci_rdf(hist, bin_centers, binsize)
    array([0.31830989, 0.15915494, 0.1061033 ])
    """
    shell_volumes = 4 * np.pi * bin_centers**2 * binsize
    rdf = hist / shell_volumes
    rdf /= np.sum(rdf)  # Normalize RDF
    return rdf


def foci_histogram(
    tags: List[str],
    pairs_dist: np.ndarray,
    binsize: float
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Generate histograms and RDFs for pairwise distances.

    This function computes histograms and radial distribution functions (RDFs)
    for all pairwise distances of large monomers along a heterogeneous polymer.
    It uses the same number of bins and bin range for all pairs. The lower bin
    range is 0 while the upper one is the maximum of physical distances between
    any pairs.

    Parameters
    ----------
    tags : List[str]
        A list of strings containing pair names.
    pairs_dist : np.ndarray
        An array of shape `(n_pairs, n_frames)` containing the time series of
        pairwise distances for each pair.
    binsize : float
        The width of each histogram bin.

    Returns
    -------
    Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]
        - `pair_hists`: A dictionary where keys are pair tags, and values are
          histograms of pair distances. Includes a special key 'bin_center'
          for bin centers.
        - `pair_rdfs`: A dictionary where keys are pair tags, and values are
          RDFs of pair distances. Includes a special key 'bin_center' for bin
          centers.

    Examples
    --------
    >>> tags = ['pair1', 'pair2']
    >>> pairs_info = np.array([[0, 1, 2], [1, 2, 3]])
    >>> pairs_dist = np.array([[1.0, 1.5], [2.0, 2.5]])
    >>> binsize = 0.5
    >>> foci_histogram(tags, pairs_info, pairs_dist, binsize)
    ({'pairDistHistFoci-pair1': array([...]), ...},
     {'pairDistRdfFoci-pair1': array([...]), ...})
    """
    # Determine histogram range and bin centers
    # the maximum pair distance used as upper limit of histogram range in all
    # pair histogram. This value is rounded to be a multiple of binsize.
    upper_range = np.ceil(np.max(pairs_dist) / binsize) * binsize
    max_range = (0, upper_range)
    max_bins = int(np.ceil((max_range[1] - max_range[0]) / binsize))
    bin_edges = np.arange(max_range[0], max_range[1] + binsize, binsize)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # Compute histograms for all pairs
    hists = [
        np.histogram(data, bins=max_bins, range=max_range)[0]
        for data in pairs_dist
    ]
    hist_tags = [f"pairDistHistFoci-{tag}" for tag in tags]
    pair_hists = dict(zip(hist_tags, hists))
    pair_hists['bin_center'] = bin_centers

    # Compute RDFs for all pairs
    rdfs = [foci_rdf(hist, bin_centers, binsize) for hist in hists]
    rdf_tags = [f"pairDistRdfFoci-{tag}" for tag in tags]
    pair_rdfs = dict(zip(rdf_tags, rdfs))
    pair_rdfs['bin_center'] = bin_centers

    return pair_hists, pair_rdfs


def whole_dist_mat_foci(
    whole_path: str,
    whole_info: ParserType
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate time series, histograms, and RDFs of foci distances from a
    distance matrix.

    This function processes the distance matrix of large monomers (foci) in a
    polymer system. It computes time series of pairwise distances, histograms,
    and radial distribution functions (RDFs) for all foci pairs.

    Parameters
    ----------
    whole_path : str
        Path to the distance matrix file. The matrix should be a 3D array where
        the first axis represents time frames.
    whole_info : ClusterParserInstance
        A instance from a parser class in one of the projects that needs
        cluster information about particle entities. The parser instance
        contains metadata about the molecular dynamics simulation project,
        including the number of large and small monomers and monomer
        properties.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        - `pair_hists`: DataFrame where columns are pair tags, and values are
          pair histograms.
        - `pair_rdfs`: DataFrame where columns are pair tags, and values are
          pair RDFs.
        - `pair_tseries`: DataFrame where columns are pair tags, and values are
          time series of pair distances.

    Examples
    --------
    >>> whole_path = "distance_matrix.npy"
    >>> whole_info = TransFociCub(dmon_small=1.0, nmon_large=4, nmon_small=6)
    >>> pair_hists, pair_rdfs, pair_tseries = \
            whole_dist_mat_foci(whole_path, whole_info)

    Notes
    -----
    The distance matrix is assumed to be upper-triangular for efficiency. Tags
    are generated in the format "locXlocYgenDistZ" where X and Y are monomer
    indices, and Z is their genomic distance.
    """
    # Load the distance matrix
    whole = np.load(whole_path)

    # Extract parameters from whole_info
    dmon_small = getattr(whole_info, 'dmon_small')
    nmon_large = getattr(whole_info, 'nmon_large')
    nmon_small = getattr(whole_info, 'nmon_small')

    # Extract genomic positions and compute pair information
    diagonal_idxs = np.diag_indices(nmon_large)
    genomic_pos = whole[0, diagonal_idxs[0], diagonal_idxs[1]]
    pairs_info = foci_info(genomic_pos, nmon_large, nmon_small)

    # Generate tags for each pair
    tags = [
        f"loc{x[0]}loc{x[1]}genDist{x[2]}"
        for x in pairs_info[:, :3]
    ]

    # Extract upper-triangular indices for distances
    triu_indices = np.triu_indices(nmon_large, k=1)
    pairs_dist = whole[:, triu_indices[0], triu_indices[1]].T

    # Compute histograms and RDFs
    binsize = 0.5 * dmon_small
    pair_hists_dict, pair_rdfs_dict = foci_histogram(tags, pairs_dist, binsize)

    # Convert results to DataFrames
    pair_hists = pd.DataFrame.from_dict(pair_hists_dict)
    pair_rdfs = pd.DataFrame.from_dict(pair_rdfs_dict)

    # Time series DataFrame
    tseries_tags = [f"pairDistTFoci-{tag}" for tag in tags]
    pair_tseries = pd.DataFrame(dict(zip(tseries_tags, pairs_dist)))

    return pair_hists, pair_rdfs, pair_tseries


def enforce_single_patch_dir_contact(
    contact_matrix: np.ndarray
) -> np.ndarray:
    """
    Enforces single-patch contact for each column in a contact matrix.

    This function modifies a binary contact matrix to ensure that each patch
    (column) is only in contact with one monomer (row). If a patch is in
    contact with multiple monomers, only the first (topmost) contact is
    preserved, and the rest are set to zero.

    Parameters
    ----------
    contact_matrix : np.ndarray
        Binary contact matrix of shape `(n_monomers, n_patches)`, where 1
        indicates a direct contact between a monomer and a patch.

    Returns
    -------
    contact_matrix: np.ndarray
        Modified contact matrix where all values after the first contact in
        each column are set to 0.

    Raises
    ------
    ValueError
        If the input matrix is not 2D.

    Notes
    -----
    Each H-NS patch can be in contact with one monomer. However, visual
    inspection has revealed that a single patch can be in contact with more
    than one monomers, which is an artifacts and does not have physical
    grounds. The inspection has revealed these monomers are adjacent monomers
    along the chain backbone. Therefore, it is decided to only consider the
    contact of the H-NS patch with the first of these monomers as a real
    contact and consider the rest of contacts as artifacts. This means setting
    the non-zero elements below the first non-zero elements along a column
    equal to zero.

    Examples
    --------
    >>> contact_matrix = np.array([[1, 0, 1],
    ...                            [1, 1, 0],
    ...                            [0, 1, 0]])
    >>> enforce_single_patch_dir_contact(contact_matrix)
    array([[1, 0, 1],
           [0, 1, 0],
           [0, 0, 0]])
    """
    if contact_matrix.ndim != 2:
        raise ValueError("Input contact matrix must be a 2D array.")

    # Create a copy to avoid modifying the input matrix in place
    contact_matrix = contact_matrix.copy()

    # Iterate over columns to enforce single-patch contact
    for col_idx in range(contact_matrix.shape[1]):
        # Find the first contact index in the column
        first_contact_idx = np.argmax(contact_matrix[:, col_idx] != 0)
        # Zero out all contacts below the first contact
        contact_matrix[first_contact_idx + 1:, col_idx] = 0

    return contact_matrix


def genomic_distance(
    atoms_a: np.ndarray,
    atoms_b: np.ndarray,
    topology: TopologyT,
    num_backbone_atoms: int
) -> np.ndarray:
    """
    Calculate the genomic distance between pairs of atoms along a polymer
    backbone.

    Genomic distance is the intra-monomer distance between atoms along the
    backbone of a polymer. It is defined based on the polymer topology
    ('linear' or 'ring').

    Parameters
    ----------
    atoms_a : np.ndarray
        Array of atom indices in the first group along the backbone, with shape
        `(num_backbone_atoms,)` or `(num_backbone_atoms, 1)`.
    atoms_b : np.ndarray
        Array of atom indices in the second group along the backbone, with
        shape `(num_backbone_atoms,)` or `(num_backbone_atoms, 1)`.
    topology : TopologyT
        Topology of the polymer ('linear' or 'ring').
    num_backbone_atoms : int
        Total number of atoms along the backbone.

    Returns
    -------
    distances: np.ndarray
        Array of element-wise genomic distances between corresponding pairs of
        atoms along the backbone.

    Raises
    ------
    ValueError
        If `atoms_a` and `atoms_b` have different shapes.
        If `topology` is not 'linear' or 'ring'.

    Notes
    -----
    - For 'linear' topology, the genomic distance is the absolute difference
      between atom indices.
    - For 'ring' topology, the genomic distance is the minimum of the direct
      distance and the distance around the circle.
    - Genomic distance is measured in terms of the number of bonds, not atoms.
      It is 1 unit more than the number of atoms between a pair of atoms.

    Examples
    --------
    >>> atoms_a = np.array([0, 2, 4])
    >>> atoms_b = np.array([1, 3, 5])
    >>> genomic_distance(atoms_a, atoms_b, topology='linear',
                         num_backbone_atoms=6)
    array([1, 1, 1])
    >>> genomic_distance(atoms_a, atoms_b, topology='ring',
                         num_backbone_atoms=6)
    array([1, 1, 1])
    """
    if atoms_a.shape != atoms_b.shape:
        raise ValueError("'atoms_a' and 'atoms_b' must have the same shape.")
    if topology not in ('linear', 'ring'):
        raise ValueError(
            f"Invalid topology: '{topology}'. Must be 'linear' or 'ring'.")

    # Define topology-specific distance calculation
    topology_functions: Dict[TopologyT, Callable] = {
        'linear': lambda a, b: np.abs(a - b),
        'ring': lambda a, b, n: np.minimum(np.abs(a - b), n - np.abs(a - b)),
    }
    distances: np.ndarray = \
        topology_functions[topology](atoms_a, atoms_b, num_backbone_atoms)
    return distances


def hns_binding(
    direct_contact: np.ndarray,
    topology: TopologyT,
    cis_threshold: int = 4,
) -> Tuple[Dict[str, int], np.ndarray]:
    """
    Calculate the binding statistics of H-NS proteins in a polymer system.

    This function evaluates the binding behavior of H-NS proteins, which
    consist of a core and two patches capable of binding monomers in the
    polymer. It computes statistics for free, bridged, and dangling H-NS
    proteins, and determines genomic distances between monomers bridged by H-NS
    cores.

    The genomic distance is classified as 'cis' or 'trans' based on the
    `cis_threshold`, and a histogram of genomic distances (loop sizes) is
    updated accordingly.

    Parameters
    ----------
    direct_contact : np.ndarray
        A binary matrix of shape `(n_mon, n_hpatch)` where:
        - `n_mon` is the number of monomers.
        - `n_hpatch` is the number of H-NS patches (twice the number of H-NS
          cores).
        - Each element `(i, j)` is `1` if monomer `i` is in contact with patch
          `j`, and `0` otherwise.

    topology : TopologyT
        The topology of the polymer, either 'linear' or 'ring'.

    cis_threshold : int, optional
        The genomic distance threshold (in bonds) below which an H-NS binding
        is classified as 'cis'. Default is `4`.

    Returns
    -------
    Tuple[Dict[str, int], np.ndarray]
        - `stats`: A dictionary containing various binding statistics:
          - 'nBoundHnsPatch': Total number of monomer-patch contacts.
          - 'nEngagedHnsPatch': Total number of engaged patches.
          - 'nFreeHnsPatch': Total number of free patches.
          - 'nFreeHnsCore': Total number of free H-NS cores.
          - 'nBridgeHnsCore': Total number of bridging H-NS cores.
          - 'nDangleHnsCore': Total number of dangling H-NS cores.
          - 'nCisHnsCore': Total number of cis-binding H-NS cores.
          - 'nTransHnsCore': Total number of trans-binding H-NS cores.

        - `hist`: A 1D array representing a histogram of genomic distances
          (loop sizes).

    Raises
    ------
    ValueError
        If the provided topology is invalid ('linear' or 'ring' are the only
        valid options).

    Notes
    -----
    - An H-NS core is considered free if neither of its two patches is bound.
    - A core is considered bridging if both patches are bound to different
      monomers.
    - A core is considered dangling if only one of its patches is bound.
    - Genomic distances are computed between monomer pairs bridged by H-NS
      cores.
    - For a 'linear' topology, the genomic distance is the absolute difference
      between monomer indices.
    - For a 'ring' topology, the genomic distance is the minimum of the direct
      distance and the distance around the circle.
    """
    invalid_keyword(topology, ['linear', 'ring'])
    n_mon, n_hpatch = direct_contact.shape
    n_hcore = n_hpatch // 2

    # Initialize stats and histogram
    stats: Dict[str, int] = {
        'nBoundTHnsPatch': 0,
        'nEngagedHnsPatch': 0,
        'nFreeHnsPatch': 0,
        'nFreeHnsCore': 0,
        'nBridgeHnsCore': 0,
        'nDangleHnsCore': 0,
        'nCisHnsCore': 0,
        'nTransHnsCore': 0,
    }
    max_gen_distance = {'linear': n_mon, 'ring': n_mon // 2 + 1}
    hist = np.zeros(max_gen_distance[topology], dtype=int)

    # Compute patch-level and core-level stats
    d_cont_per_hpatch = np.sum(direct_contact, axis=0)
    stats['nBoundHnsPatch'] = np.sum(direct_contact)
    stats['nEngagedHnsPatch'] = np.count_nonzero(d_cont_per_hpatch)
    stats['nFreeHnsPatch'] = n_hpatch - stats['nEngagedHnsPatch']

    free_cores = \
        np.sum((d_cont_per_hpatch[0::2] == 0) & (d_cont_per_hpatch[1::2] == 0))
    bridged_cores = \
        np.sum((d_cont_per_hpatch[0::2] > 0) & (d_cont_per_hpatch[1::2] > 0))
    dangling_cores = n_hcore - free_cores - bridged_cores

    stats['nFreeHnsCore'] = free_cores
    stats['nBridgeHnsCore'] = bridged_cores
    stats['nDangleHnsCore'] = dangling_cores

    # Enforce single-patch binding and calculate genomic distances
    single_patch_dir_contact = enforce_single_patch_dir_contact(direct_contact)
    cont_m_hpatch = generate_contact_matrix(single_patch_dir_contact)
    cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
    cont_m_m = np.matmul(cont_m_hcore, cont_m_hcore.T)
    cont_m_m_nonzeros = np.array(np.where(np.triu(cont_m_m, 1) > 0)).T
    # Calculate genomic distances for monomer pairs bridged by an H-NS protein:
    genomic_distances = genomic_distance(
        cont_m_m_nonzeros[:, 0],
        cont_m_m_nonzeros[:, 1],
        topology,
        n_mon
    )
    # A 2D array with shape `(n_pairs, 3)`: Column 0 is the index of the first
    # monomer in the pair. Column 1 is the index of the second monomer in the
    # pair. And, column 3 is the genomic distance between the two monomers:
    m_m_gen_dist = np.column_stack((cont_m_m_nonzeros, genomic_distances))

    # Classify bindings and update histogram
    cis_bindings = \
        np.count_nonzero(
            (m_m_gen_dist[:, 2] > 0) & (m_m_gen_dist[:, 2] <= cis_threshold)
            )
    trans_bindings = np.count_nonzero(m_m_gen_dist[:, 2] > cis_threshold)
    stats['nCisHnsCore'] = cis_bindings
    stats['nTransHnsCore'] = trans_bindings
    np.add.at(hist, m_m_gen_dist[:, 2], 1)

    return stats, hist


def generate_mon_bind_direct(
    mon_patch_dir: np.ndarray,
    patch_per_binder: int
) -> np.ndarray:
    """
    Create a monomer-binder direct contact matrix from a monomer-patch direct
    contact matrix.

    Parameters
    ----------
    mon_patch_dir : np.ndarray
        Binary matrix of shape `(n_monomer, patch_per_binder * n_binder)`,
        where each element `(i, j)` is `1` if monomer `i` is in contact
        with patch `j`, and `0` otherwise.
    patch_per_binder : int
        Number of patches per binder.

    Returns
    -------
    np.ndarray
        Binary matrix of shape `(n_monomer, n_binder)` where each element
        `(i, j)` is `1` if monomer `i` is in contact with binder `j`, and `0`
        otherwise.

    Raises
    ------
    ValueError
        If the input matrix dimensions are invalid or not a multiple of
        `patch_per_binder`.
    """
    if mon_patch_dir.ndim != 2:
        raise ValueError("Input monomer-patch matrix must be 2D.")
    if mon_patch_dir.shape[1] % patch_per_binder != 0:
        raise ValueError(
            f"Number of columns ({mon_patch_dir.shape[1]}) is not a multiple"
            f" of patch_per_binder ({patch_per_binder})."
        )

    # Reshape and apply any operation to determine binder contact
    n_bind = mon_patch_dir.shape[1] // patch_per_binder
    mon_bind_dir = np.zeros((mon_patch_dir.shape[0], n_bind), dtype=np.integer)

    for i in range(n_bind):
        start_col = patch_per_binder * i
        end_col = start_col + patch_per_binder
        binder_patches = mon_patch_dir[:, start_col:end_col]
        binder = np.any(binder_patches, axis=1).astype(int)
        mon_bind_dir[:, i] = binder

    return mon_bind_dir


def split_binder_matrix(m_ij: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split the monomer-binder direct contact matrix into two matrices based
    on whether a binder is a bridger or a dangler.

    Parameters
    ----------
    m_ij : np.ndarray
        Binary matrix of shape `(n_monomer, n_binder)` where each element
        `(i, j)` is `1` if monomer `i` is in contact with binder `j`, and `0`
        otherwise.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        - `d_ij`: Matrix representing dangling binders (sum = 1 per column).
        - `b_ij`: Matrix representing bridging binders (sum = 2 per column).

    Raises
    ------
    ValueError
        If the input matrix is not binary.
    """
    if not np.array_equal(m_ij, m_ij.astype(bool)):
        raise ValueError("Input matrix must be binary.")

    column_sums = np.sum(m_ij, axis=0)
    danglers = column_sums == 1
    bridgers = column_sums == 2

    d_ij = m_ij * danglers
    b_ij = m_ij * bridgers

    return d_ij, b_ij


def find_binder_clusters_need_attention(
    m_ij: np.ndarray,
    topology: TopologyT
) -> np.ndarray:
    """
    Identify clusters of binders in a polymer system based on their
    binding to the same or adjacent monomers.

    This function sets all diagonal elements of the output matrix to 1,
    indicating each binder forms a cluster with itself. It then iterates over
    all pairs of binders, checking for direct connections to the same monomer
    or connections via adjacent monomers.

    Parameters
    ----------
    m_ij : np.ndarray
        Binary matrix of shape `(n_monomer, n_binder)` where each element
        `(i, j)` is `1` if monomer `i` is in contact with binder `j`, and `0`
        otherwise.
    topology : TopologyT
        The topology of the polymer ('linear' or 'ring').

    Returns
    -------
    np.ndarray
        Binary square matrix of shape `(n_binder, n_binder)`. Each element
        `(i, j)` is `1` if binders `i` and `j` are part of the same cluster,
        and `0` otherwise.

    Notes
    -----
    - Diagonal elements indicate whether a binder is bound (`1`) or unbound
      (`0`).
    - Adjacent monomers in a 'ring' topology include the first and last
      monomers.
    """
    print('This function has not yet tested')
    _, n_binders = m_ij.shape
    b_ij = np.zeros((n_binders, n_binders), dtype=int)

    # Diagonal: 1 if the binder is bound to at least one monomer, 0 otherwise
    b_ij[np.arange(n_binders), np.arange(n_binders)] = \
        m_ij.any(axis=0).astype(int)

    # Same-monomer clusters
    binder_pairs = m_ij.T @ m_ij  # Shape: (n_binder, n_binder)
    b_ij |= binder_pairs > 0

    # Adjacent-monomer clusters (linear topology)
    adjacent_binders = m_ij[:-1] | m_ij[1:]
    adjacent_pairs = adjacent_binders.T @ adjacent_binders
    b_ij |= adjacent_pairs > 0

    # Wrap-around clusters (ring topology)
    if topology == 'ring':
        wrap_binders = m_ij[0] | m_ij[-1]
        wrap_pairs = np.outer(wrap_binders, wrap_binders)
        b_ij |= wrap_pairs

    return b_ij


def find_binder_clusters_need_attention_old(
    m_ij: np.ndarray,
    topology: TopologyT
) -> np.ndarray:
    """
    Identify clusters of binders in a polymer system based on their
    binding to the same or adjacent monomers.

    This function sets all diagonal elements of the output matrix to 1,
    indicating each binder forms a cluster with itself. It then iterates over
    all pairs of binders, checking for direct connections to the same monomer
    or connections via adjacent monomers.

    Parameters
    ----------
    m_ij : np.ndarray
        Binary matrix of shape `(n_monomer, n_binder)` where each element
        `(i, j)` is `1` if monomer `i` is in contact with binder `j`, and `0`
        otherwise.
    topology : TopologyT
        The topology of the polymer ('linear' or 'ring').

    Returns
    -------
    np.ndarray
        Binary square matrix of shape `(n_binder, n_binder)`. Each element
        `(i, j)` is `1` if binders `i` and `j` are part of the same cluster,
        and `0` otherwise.

    Notes
    -----
    - Diagonal elements indicate whether a binder is bound (`1`) or unbound
      (`0`).
    - Adjacent monomers in a 'ring' topology include the first and last
      monomers.
    """
    print('This function has not yet tested')
    n_binders = m_ij.shape[1]
    n_monomers = m_ij.shape[0]
    b_ij = np.zeros((n_binders, n_binders), dtype=int)
    # diagonal elements have meaning here; if it is 0 means unbound hns/binder
    # but if 1 means it bridged two monomers or dangle from on monomer.
    np.fill_diagonal(b_ij, 1)
    # Check for binders connected to the same monomer
    for i in range(n_binders):
        for j in range(i + 1, n_binders):
            if np.any(m_ij[:, i] & m_ij[:, j]):
                b_ij[i, j] = b_ij[j, i] = 1

    # Check for binders connected to adjacent monomers
    for monomer in range(n_monomers - 1):
        adjacent_binders = np.where(m_ij[monomer] | m_ij[monomer + 1])[0]
        for binder in adjacent_binders:
            b_ij[binder, adjacent_binders] = 1

    if topology == 'ring':
        adjacent_binders = np.where(m_ij[0] | m_ij[-1])[0]
        for binder in adjacent_binders:
            b_ij[binder, adjacent_binders] = 1

    return b_ij
