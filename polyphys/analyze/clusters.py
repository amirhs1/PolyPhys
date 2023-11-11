"""A submodule to detect and analyze particle clusters in a given system. It is
worth noting there is a difference analyzing particle clusters and
trajectory (or ensemble) clusters.
Below, different algorithms found in literature are implemented in scientific
Python.
"""
from typing import Dict, List, Tuple, Optional, Union
from itertools import combinations
from collections import defaultdict
import numpy as np
import pandas as pd
import numpy.linalg as npla

from ..manage.parser import TransFociCub, TransFociCyl


def apply_pbc_orthogonal(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: Dict[int, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """Updates the periodic boundary condition (PBC) information in the
    dimensions (keys) by the length (values) passed by `pbc`. The key
    are indexes of the `pbc_lengths` and `pbc_lengths_inverse` arrays.

    Parameters
    ----------
    pbc_lengths: numpy.ndarray
        An array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        An array of the inverses of lengths/sizes in different directions.
    pbc: dict
        A dictionary of dimensions (keys) and lengths (values) in which
        the pbc exists.

    Return
    ------
    pbc_lengths: numpy.ndarray
        The updated array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        The updated array of the inverses of lengths/sizes in different
        directions.
    """
    for dim, length in pbc.items():
        pbc_lengths[dim] = length
        pbc_lengths_inverse[dim] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def self_dist_array(
    positions: np.ndarray,
    pbc: Optional[Dict[int, float]] = None
) -> np.ndarray:
    """
    Calculate pairwise distances between all pairs of atoms for a given list
    of atoms' `positions`.

    Parameters
    ----------
    positions : np.ndarray, shape (n_atoms, n_dims)
        An array of atomic positions, where each row represents the position of
        an atom in n_dims
        Cartesian dimensions.
    pbc : dict, optional
        A dictionary of dimensions (keys) and lengths (values) in which the
        periodic boundary conditions exist. Keys are integers where 0 means 'x'
        dimension, 1 means 'y' dimension, and 2 means 'z'
        dimension in the Cartesian coordinate system. Default is None.

    Returns
    -------
    np.ndarray, shape (n_atoms, n_atoms)
        A symmetric matrix where each element (i, j) is the center-to-center
        distance between atoms i and j.

    References
    ----------
    Albanie S. (2019) Euclidean Distance Matrix Trick. In: Advances in Computer
    Vision. EECV 2018. Lecture Notes in Computer Science, vol 11131. Springer,
    Cham. https://doi.org/10.1007/978-3-030-01446-4_2
    """
    n_atoms, n_dims = positions.shape
    # Calculate pairwise distance matrix
    pos_t = positions.T
    pos_i = np.reshape(pos_t, (n_dims, 1, n_atoms))
    pos_j = np.reshape(pos_t, (n_dims, n_atoms, 1))
    dist_ji = pos_j - pos_i

    # Calculate PBC lengths and inverses, and apply PBCs if PBCs are present
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
            )
        pbc_lengths = pbc_lengths.reshape(n_dims, 1, 1)
        pbc_lengths_inverse = pbc_lengths_inverse.reshape(n_dims, 1, 1)
        dist_ji -= pbc_lengths * np.around(pbc_lengths_inverse * dist_ji)

    # Calculate squared distances and return square root of matrix
    dist_ji = np.sqrt(np.sum(np.square(dist_ji), axis=0))
    np.fill_diagonal(dist_ji, 0)
    return dist_ji


def self_dist_array_opt(
    positions: np.ndarray,
    pbc: Optional[Dict[int, float]] = None
) -> np.ndarray:
    """Calculate all possible distances between all the pair of atoms for a
    given list of atoms' `positions`. In 3 dimensions, this method is not
    much faster the `self_dist_array`.

    To-do
    -----
    Check the pbc applied below; is it correct?

    Parameters
    ----------
    positions : np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different dimensions
        or components.

    pbc : dict, defualt None
        A dictionary of dimensions (keys) and lengths (values) in which
        the pbc exists. keys are integers where 0 means 'x' dimension,
        1 means 'y' dimension, and 2 means 'z' dimension in the cartesian
        coordinate system.

    Return
    ------
    np.ndarray
        A square matrix of size n_atoms * n_atoms in which each element are
        the center-to-center distances between different pairs of atoms. All
        the diagonal elements of this matrix is 0.

    References
    ----------
    Albanie S. (2019) Euclidean Distance Matrix Trick. In: Advances in Computer
    Vision. EECV 2018. Lecture Notes in Computer Science, vol 11131. Springer,
    Cham. https://doi.org/10.1007/978-3-030-01446-4_2
    """
    n_atoms, n_dims = positions.shape
    # Redefining the shapes of different arrays to combine linear algebra
    # with numpy broadcasting.
    gram = np.matmul(positions, positions.T)
    dist_sq = np.diag(gram).reshape(n_atoms, 1) + \
        np.diag(gram).reshape(1, n_atoms) - 2 * gram

    # Calculate PBC lengths and inverses, and apply PBCs if PBCs are present
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
        )
        pbc_lengths = pbc_lengths**2
        pbc_lengths_inverse = pbc_lengths_inverse**2
        dist_sq = dist_sq - pbc_lengths * np.around(
            pbc_lengths_inverse * dist_sq
            )
    return np.sqrt(dist_sq)


def find_direct_contacts(
    dist: np.ndarray,
    cutoff: float,
    inclusive: bool = True
) -> np.ndarray:
    """
    Determines whether pairs of atoms are in direct contact based on a cutoff
    distance.

    Direct contact is defined as the center-to-center distance between two
    atoms being less than or equal to the cutoff distance, or less than the
    cutoff distance depending on the value of the `inclusive` parameter.

    Parameters
    ----------
    dist : np.ndarray, shape (n_atoms, m_atoms)
        A matrix of center-to-center distances between all pairs of atoms.
    cutoff : float
        The cutoff distance for defining direct contact.
    inclusive : bool, optional
        If True, includes pairs of atoms where the center-to-center distance
        is equal to the cutoff distance. If False, only includes pairs where
        the distance is strictly less than the cutoff distance. Default is
        True.

    Returns
    -------
    np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where each element (i, j) is 1 if the center-to-center
        distance between atoms i and j is within the cutoff distance, and 0
        otherwise. If the input dist matrix is the symmteric self-distance
        sqaure matrix, then the output is also a symmetric square matrix,
        with all diagonal elements set to 1.

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://aip.scitation.org/doi/10.1063/1.454720 .
    """
    if cutoff <= 0:
        raise ValueError("'cutoff' must be a positive number.")
    if inclusive:
        direct_contacts = dist <= cutoff
    else:
        direct_contacts = dist < cutoff

    return np.asarray(direct_contacts, dtype=int)


def generate_contact_matrix(direct_contacts: np.ndarray) -> np.ndarray:
    """
    Generates a symmetric matrix of all contacts between pairs of atoms
    including both direct and indirect contacts, from a matrix of direct
    contacts between pairs of atoms.

    The `direct_contacts` matrix is a binary matrix of the shape
    (n_atoms, m_atoms), where any of 'm_atoms' along the columns any be in
    direct contact with any of 'n_atoms' along the rows.
    If the `direct_contact` matrix is the symmetric self-direct square matrix
    of a system with 'n_atoms', then the output is also symmetric.
    See `find_direct_contacts` function.

    Parameters
    ----------
    direct_contacts : np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where each element (i, j) is 1 if atoms i and j have a
        direct contact, and 0 otherwise.

    Returns
    -------
    contact_matrix : np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where each element (i, j) is 1 if atoms i and j have a
        contact, and 0 otherwise.

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics,
    89(1), 668-676. https://doi.org/10.1063/1.454720
    """
    _, n_cols = direct_contacts.shape  # shape: n_atoms * m_atoms
    contact_matrix = direct_contacts.copy()
    # Loop over all pairs of column indices
    for i, j in combinations(range(n_cols), 2):
        # Check if columns i and j have any common nonzero elements
        if np.any(contact_matrix[:, i] & contact_matrix[:, j]):
            # If so, set all nonzero elements in columns i and j to 1
            transient = contact_matrix[:, i] | contact_matrix[:, j]
            contact_matrix[:, i] = transient
            contact_matrix[:, j] = transient

    return contact_matrix


def random_direct_contacts(
    n_atoms: int,
    seed: Optional[int] = 0,
    save_to: Optional[Tuple[str, str]] = None
) -> np.ndarray:
    """Generates a direct-contact matrix with size `n_atoms` which can be
    used to test a clustering algorithm. See `direct_contact`function
    description for more info.

    All the diagonal elements of the direct contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many
    repeated columns.

    Note
    ----
    This function can be used for testing the `count_bonds` and
    `count_clusters` methods.

    Parameters
    ----------
    n_atoms: int
        Number of atoms or nodes.
    seed: int, default 0
        The seed for the random number generator.
    save_to: tuple of str, default None
        A tuple of two strings in which the first element is an absolute or
        a relative path of a directory to which the random direct
        contact matrix is saved and the second one is the matrix name.

    Return
    ------
    dir_contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.
    """
    # Creating direct-contacts matrix
    np.random.seed(seed)
    dir_contacts = np.random.binomial(n=1, p=0.35, size=[n_atoms * n_atoms])
    dir_contacts = np.reshape(dir_contacts, (n_atoms, n_atoms))
    dir_contacts_up = np.triu(dir_contacts)  # upper triangle
    dir_contacts_diag = np.identity(n_atoms, dtype=int)
    # dropping diagonal elements:
    dir_contacts_up = dir_contacts_up - dir_contacts_diag
    # ensuring all diagonal elements are either 0 or 1:
    dir_contacts_up[dir_contacts_up < 0] = 0
    dir_contacts_low = dir_contacts_up.T  # lower triangle: symmetry!
    dir_contacts = dir_contacts_up + dir_contacts_low + dir_contacts_diag
    # Creating the full direct contacts matrix
    dir_contacts = generate_contact_matrix(dir_contacts)
    if save_to is not None:
        output = save_to[0] + \
            f"random_direct_contact-natoms{n_atoms}" + save_to[1]
        np.save(output, dir_contacts)
    return dir_contacts


def count_foci_bonds(dir_contacts: np.ndarray) -> np.ndarray:
    """Infers the number and size of bonds from the direct contact matrix
    `direct_contacts`.

    The direct contact matrix is symmetric and its diagonal elements are 1.
    The column- or row-wise sum of the matrix shows the number of bonds
    an atom has with itself and other atoms. A self-bond happens due to the
    definition of the direct contact matrix (in which the diagonal elements
    are 1). Such self-bonds are subtracted below to give the correct number of
    bonds, thus number of bonds per atom ranges from 0 to $n_{atoms}-1$.

    Parameters
    ----------
    dir_contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    Return
    ------
    cluster_list: array-like
        A 1D array in which the indexes are "bond" sizes and the values
        are frequency of bonds with that size.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = dir_contacts.shape
    bond_list = np.zeros(n_atoms, dtype=int)
    # removing self-bonds:
    bonds = dir_contacts.sum(axis=1) - 1
    bond_sizes, bond_counts = np.unique(bonds, return_counts=True)
    np.put(bond_list, bond_sizes, bond_counts)
    return bond_list


def count_foci_clusters(contacts: np.ndarray) -> np.ndarray:
    """Infers the number and size of clusters from a contact matrix `contacts`.

    All the diagonal elements of the direct contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many repeated
    columns. While the number of non-zero eigenvalues of the contact
    matrix is equal to the number of clusters, the values of eigenvalues
    show the sizes of clusters. The frequencies of eigenvalues are equal to
    the frequencies of clusters with different sizes. An eigenvalue with value
    0 does not have any physical meaning, but the frequency of this eigenvalue
    shows the total number of repeated clusters with sizes equal to or larger
    than 2 (that is equal to the number of dependant columns in the
    contact matrix in the linear algebra lingo).

    Parameters
    ----------
    contacts: np.ndarray
        A matrix of size n_atoms * n_atoms in which each element is either
        0 or 1 and indicates whether a contact is between a pair of
        atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    Return
    ------
    cluster_list: array-like
        A 1D array in which the indexes are "cluster_size - 1" sizes and the
        values are frequency of clusters with that size.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = contacts.shape
    # While eigenvalue 0 is not physically important, eigenvalue n_atoms
    # is. Here, a cluster list with values ranging from 0 to n_atoms is created
    # to properly collect all the eigenvalues.
    cluster_list = np.zeros(n_atoms + 1, dtype=int)
    clusters, _ = npla.eig(contacts)
    # numpy returns infinitesimal numbers instead of zero when an eigenvalue
    # is 0. Below, these eigenvalues are rounded to 0. The non-zero
    # values are already integer, so rounding does not affect them.
    clusters = np.asarray(np.round(clusters), dtype=int)
    # Sanity check for having meaning cluster sizes:
    if np.any(clusters < 0):
        raise ValueError("A cluster with size smaller than 0 found!")
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    # Above, the eigenvalue with value 0 can be interpreted in this way:
    # The system has clusters with size equal to or larger than 1
    # The frequency of eigenvalue 0 is the total number of clusters (or
    # dependant columns in the contact matrix) with size equal to
    # or larger than 2.
    np.put(cluster_list, cluster_sizes, cluster_counts)
    return cluster_list


def foci_info(
    genomic_pos: np.ndarray,
    nmon_large: int,
    nmon_small: int
) -> np.ndarray:
    """Generates a dictionary of the information about all the possible pairs
    of large monomers in a circular/ring heterogeneous polymer composed of two
    different monomer types (small and large monomers).

    This information is an array of shape(n_pairs,3) where the first two
    columns are the pair index along the backbone (a.k.a genomic position), the
    2nd one is the genomic distance (i.e. the minimum number of genomic indices
    between a pair in a circular polymer).

    Parameters
    ----------
    genomic_pos: np.ndarray
        A vector of size `nmon_large` that contains the genomic position of
        large monomers
    nmon_large: int
        Number of large monomers
    nmon_small: int
        Number of small monomers

    Return
    ------
    pair_info: np.ndarray
        An array of shape(n_pairs,3) that contains the genomic indices, and
        genomic distance of each pair.
    """
    nmon = nmon_small + nmon_large
    pairs = np.array(list(combinations(range(nmon_large), 2)), dtype=np.int_)
    genomic_idx = np.choose(pairs, genomic_pos)
    genomic_dist = np.diff(genomic_idx, axis=1)
    # Applying the topological contain on the genomic distance, i.e. finding
    # minimum distance between a pair on a circle.
    # the pair is inclusive, so their difference is one unit more than the
    # actual number of genomic units between them, but it is equal to the
    # number of bonds (not the length of bonds).
    # Number of small monomers between to large monomers:
    genomic_dist = np.minimum.reduce((nmon-genomic_dist, genomic_dist)) - 1
    pair_info = np.concatenate(
        (genomic_idx, genomic_dist),
        axis=1
        )
    return pair_info


def v_histogram(
    data: np.ndarray,
    bins: int,
    bin_range: Tuple[float, float],
    **kwargs
) -> np.ndarray:
    """Vectorize `n.histogram` so it can be mapped to a collection of data
    with different `nbins` and `range` kwargs. The `range` is divided into
    `nbins` bins, so there are `nbins` equal bins.

    Parameters
    ----------
    data: np.ndarray
        1D array-like input data
    bins: int
        Number of bins
    bin_range: float
        The lower and upper range of bins.

    Return
    ------
    hist: np.ndarray
        The frequencies of `data` in bins.
    """
    hist, _ = np.histogram(data, bins=bins, range=bin_range, **kwargs)
    return hist


def foci_rdf(
    hist: np.ndarray,
    bin_centers: np.ndarray,
    binsize: float
) -> np.ndarray:
    """Compute the radial distribution function to find the center of a
    particle in a pair given position at a radial distance r from the center of
    the other particle in that pair.

    Parameters
    ----------
    hist: np.ndarray
       The frequencies of particles in bins with equal size.
    bin_centers
        The positions of bin centers along the radial direction.
    binsize: float
        Size of each bin

    Return
    ------
    rdf: np.ndarray
        The rdf in bins with equal size.
    """
    rdf = hist / (bin_centers**2 * 4 * np.pi * binsize)
    rdf = rdf / np.sum(rdf)  # normalization
    return rdf


def foci_histogram(
    tags: List[str],
    pairs_info: np.ndarray,
    pairs_dist: np.ndarray,
    binsize: float
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Generates histograms of pair distances, use the same number of bins and
    the same bin range for all pairs.

    The lower bin range is 0 while the upper one is the maximum of physical
    distances between any pairs.

    Parameters
    ----------
    tags: list of str
        A list that contains pair names.
    pairs_info: np.ndarray
        An array of shape(n_pairs,3) that contains the genomic indices, and
        genomic distance.
    pairs_dist: np.ndarray
        An array of shape (n_pairs, n_frames) that contains the time series
        of each pair distances.
    binsize: float
        Size of each bin

    Return
    ------
    pair_hists: dict
        A dictionary in which the keys are pair tags and the values of pair
        histograms. There is also extra key named 'bin_center' which contains
        the center of bins.
    pair_hists: dict
        A dictionary in which the keys are pair tags and the values of pair
        rdfs. There is also extra key named 'bin_center' which contains
        the center of bins.
    """
    n_pairs = pairs_info.shape[0]
    # the maximum pair distance used as upper limit of histogram range in all
    # pair histogram. This value is rounded to be a multiple of binsize.
    upper_range = np.ceil(np.max(pairs_dist) / binsize) * binsize
    max_range = (0, upper_range)
    max_bins = np.ceil((max_range[1] - max_range[0]) / binsize).astype(np.int_)
    pairs_range = [max_range] * n_pairs
    pairs_bins = [max_bins] * n_pairs
    bin_edges = np.arange(max_range[0], max_range[1] + binsize, binsize)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    hists = list(map(v_histogram, pairs_dist, pairs_bins, pairs_range))
    hist_tags = ['pairDistHistFoci-' + tag for tag in tags]
    pair_hists = dict(zip(hist_tags, hists))
    pair_hists['bin_center'] = bin_centers
    pairs_bin_centers = [bin_centers] * n_pairs
    pairs_binsizes = [binsize] * n_pairs
    rdfs = list(map(foci_rdf, hists, pairs_bin_centers, pairs_binsizes))
    rdf_tags = ['pairDistRdfFoci-' + tag for tag in tags]
    pair_rdfs = dict(zip(rdf_tags, rdfs))
    pair_rdfs['bin_center'] = bin_centers
    return pair_hists, pair_rdfs


def whole_dist_mat_foci(
    whole_path: str,
    whole_info: Union[TransFociCub, TransFociCyl]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generates time series and histograms of foci distances from distance
    matrix of (large) monomers.

    A foci are a pair of large monomer. The distance matrix is upper-triangular
    matrix containing all the spatial distances between any possible foci of
    large monomers.

    Parameters
    ----------
    whole_path: str
        Path to a "whole" distance matrix.
    whole_info: TransFociT
        A TransFociT object contains information about a "whole" simulation

    Return
    ------
    pairs_hists: dict
        A dataframe in which the columns are pair tags and the values of pair
        histograms. There is also extra column named 'bin_center' that contains
        the center of bins.
    pairs_hists: dict
        A dataframe in which the columns are pair tags and the values of pair
        rdfs. There is also extra column named 'bin_center' that contains
        the center of bins.
    pairs_tseries: dict
        A dataframe in which the columns are pair tags and the values of pair
        timeseries.
    """
    whole = np.load(whole_path)
    dmon_small = whole_info.dmon_small
    nmon_large = whole_info.nmon_large
    nmon_small = whole_info.nmon_small
    diagonal_idxs = np.diag_indices(nmon_large)
    genomic_pos = whole[0, diagonal_idxs[0], diagonal_idxs[1]]
    pairs_info = foci_info(genomic_pos, nmon_large, nmon_small)
    tags = list(
        map(
            lambda x: "loc{0}loc{1}genDist{2}".format(x[0], x[1], x[2]),
            pairs_info[:, :3]
        )
    )
    triu_indices = np.triu_indices(nmon_large, k=1)
    pairs_dist = whole[:, triu_indices[0], triu_indices[1]].T
    binsize = 0.5 * dmon_small
    pair_hists, pair_rdfs = foci_histogram(
        tags, pairs_info, pairs_dist, binsize
    )
    pair_hists = pd.DataFrame.from_dict(pair_hists)
    pair_rdfs = pd.DataFrame.from_dict(pair_rdfs)
    tseries_tags = ['pairDistTFoci-' + tag for tag in tags]
    pair_tseries = dict(zip(tseries_tags, pairs_dist))
    pair_tseries = pd.DataFrame.from_dict(pair_tseries)
    return pair_hists, pair_rdfs, pair_tseries


def enforce_single_patch_dir_contact(
    contact_matrix: np.ndarray
) -> np.ndarray:
    """Given a binary contact matrix of shape (n_monomers, n_patches), where
    each element represents whether a given monomer and patch are in direct
    contact or not, sets the values of all the elements after the first
    contact (the first non-zero element) in each column to 0.

    Each H-NS patch can be in contact with one monomer. However, visual
    inspection has revealed that a single patch can be in contact with more
    than one monomers, which is an artifacts and does not have physical
    grounds. The inspection has revealed these monomers are adjacent monomers
    along the chain backbone. Therefore, it is decided to only consider the
    contact of the H-NS patch with the first of these monomers as a real
    contact and conisder the rest of contacts as artifacts. This means setting
    the non-zero elements below the first
    non-zero elements along a column equal to zero.


    Parameters
    ----------
    contact_matrix : numpy.ndarray
        The binary contact matrix of shape (n_monomers, n_patches).

    Returns
    -------
    numpy.ndarray
        The modified contact matrix, where the values of all elements after
        the first contact in each column have been set to 0.
    """
    # make a copy of the input matrix to avoid modifying it in place
    contact_matrix = contact_matrix.copy()

    # loop over columns
    for j in range(contact_matrix.shape[1]):
        # find index of first contact in column j
        i = np.argmax(contact_matrix[:, j] != 0)
        # set values of elements after first contact to 0
        contact_matrix[i+1:, j] = 0

    return contact_matrix


def genomic_distance(
    atoms_a: np.ndarray,
    atoms_b: np.ndarray,
    topology: str,
    num_backbone_atoms: int
) -> np.ndarray:
    """Calculate the intra-monomer distance, a.k.a the genomic distance,
    between atoms in `atoms_a` and `atoms_b` along the backbone of a polymer
    with the given `topology`.

    Parameters
    ----------
    atoms_a : np.ndarray
        Array of atom indices in the first group along the backbone, with shape
        (num_backbone_atoms,) or (num_backbone_atoms, 1).
    atoms_b : np.ndarray
        Array of atom indices in the second group along the backbone, with
        shape (num_backbone_atoms,) or (num_backbone_atoms, 1).
    topology : str
        Topology of the polymer ('linear' or 'ring').
    num_backbone_atoms : int
        Total number of atoms along the backbone.

    Returns
    -------
    np.ndarray
        Array of element-wise genomic distances between corresponding pairs of
        atoms along the backbone.

    Raises
    ------
    ValueError
        If `atoms_a` and `atoms_b` have different shapes or if `topology` is
        not 'linear' or 'ring'.

    Notes
    -----
    The genomic distance between two atoms is the absolute value of the
    difference in their indices along the backbone of the polymer. If the
    polymer has a 'ring' topology, the genomic distance is the minimum of the
    absolute difference and the difference in index when traversing the ring
    in the opposite direction.
    Genomic distance is measured in number of bonds, not number of atoms, so
    it is 1 unit more than the number of atoms between a pair of atoms.
    """
    if atoms_a.shape != atoms_b.shape:
        raise ValueError("'atoms_a' and 'atoms_b' must have the same shape.")
    if topology not in ('linear', 'ring'):
        raise ValueError(
            f"The genomic distance is not defined for '{topology}' topology.")
    topology_functions = {
        'linear': lambda a, b: np.abs(a - b),
        'ring': lambda a, b, n: np.minimum(np.abs(a - b), n - np.abs(a - b)),
    }
    distances: np.ndarray = \
        topology_functions[topology](atoms_a, atoms_b, num_backbone_atoms)
    return distances


def hns_genomic_distance(
    contact_nonzeros: np.ndarray,
    topology: str,
    n_mons: int,
) -> np.ndarray:
    """Calculate the intra-monomer distance between the pair of monomers
    bridged by a H-NS protein in a polymer.

    Parameters
    ----------
    contact_m_m_nonzeros: np.ndarray
        Pair of indices of non-zero elements in the monomer-monomer contact
        matrix.
    topology: str
        Polymer topology.
    n_mons: int
        Number of monomers in the polymer.

    Return
    ------
    contact_m_m_gen_dist: np.ndarray
        A matrix with 3 columns anv varing number of rows, where the first
        column contain the indices of the first monomer, second column
        contaions the indices of the second monomer, and the third one contains
        the genomic distance betweem the two monnomers in pairs of bridged
        monomers.
    """
    genomic_distances: np.ndarray = genomic_distance(
        contact_nonzeros[:, 0],
        contact_nonzeros[:, 1],
        topology,
        n_mons
    )
    contact_m_m_gen_dist: np.ndarray = np.column_stack((
        contact_nonzeros,
        genomic_distances
    ))  # axis 0: monomer i, monomer j, genomic_dist_ij
    return contact_m_m_gen_dist


def hns_binding(
    direct_contact: np.ndarray,
    topology: str,
    cis_threshold: int = 4,
    binding_stats: Dict[str, List[int]] = defaultdict(List[int]),
    loop_length_hist: np.ndarray = np.array([])
) -> Tuple[Dict[str, List[int]], np.ndarray]:
    """Calculate the binding statistics of H-NS proteins to monomers.

    A H-NS protein consist of a core and two patches at the poles of the core.
    In the `direct_contact` matrix, every two columns represent two patches of
    the same H-NS proteins. We use this data to determine whether an H-NS
    protein is in the unbound, dangled, or bridged mode.

    In principal, a patch cannot be in direct contact with more than one
    monomer. However, due to many-body nature of the system and non-directional
    nature of the forces, a patch can be in contact with 2 or 3 monomers in
    very rare situations. Consequently, the some of rows along each column in
    the direct_contact matrix can be more than one.

    Parameters
    ----------
    direct_contact: np.ndarray, shape (n_mon, n_hpatch)
        A binary matrix where each element (i, j) is 1 if atoms i and j have a
        contact, and 0 otherwise.

    topology: str
        Topology of the polymer containing monomers.

    cis_threshold: int, default 4
        The genomic distance in number of bonds (not number of monomers) less
        than or equal (inclusive) to which a H-NS is of cis-binding type.

    binding_stats: dict
        A dictionary of binding statistics of H-NS proteins to which new
        statisitics is appended.

    loop_length_hist: list of array
        A list of arrays where each array is of shape (n_bridged_pairs, 3).
        Each row of the an array contaion the index of first monomer, index of
        second one and the genomic distance between them.

    Return
    ------
    binding_stats: dict
        A new or updated dictionary contains the statistical data about the
        binding of H-NS proteins to monomers.

    loop_size_hist: np.ndarray
        A new or updated arrat contains the histogram of loop sizes.
    """
    # Initialize results if not provided
    n_mon, n_hpatch = direct_contact.shape
    n_hcore = n_hpatch // 2
    if topology not in ('ring', 'linear'):
        raise ValueError(
            f"The genomic distance is not defined for '{topology}' topology"
            )
    # For both, ring with odd or even number of monomers, the maximum
    # intra-chain distance is the same.
    if loop_length_hist.size == 0:
        max_gen_distance = {
            'linear': n_mon,
            'ring': n_mon//2+1
        }
        loop_length_hist = np.zeros(max_gen_distance[topology])
    # A H-NS patch can bind more than one monomers, so the total number of
    # of bindings between H-NS patches and monomers can be more than the total
    # number of H-Ns patches. However, the total number of engaged and
    # free patches is always less than or equal to the total number of H-NS
    # patches:
    # These two 'n_m_hpatch_bound' and 'n_hpatch_engaged' have contributions
    # from in valid double and triple binding of a patch.
    binding_stats['n_m_hpatch_bound'].append(np.sum(direct_contact))
    d_cont_per_hpatch = np.sum(direct_contact, axis=0)
    binding_stats['n_hpatch_engaged'].append(
        np.count_nonzero(d_cont_per_hpatch))
    binding_stats['n_hpatch_free'].append(
        n_hpatch - binding_stats['n_hpatch_engaged'][-1])
    # The correct statistics for counting briding, dangling, and free H-NS
    # proteins even if a patch is incorrectely in contact with 2 or more
    # monomers:
    binding_stats['n_hcore_free'].append(np.sum(
        (d_cont_per_hpatch[0::2] == 0) & (d_cont_per_hpatch[1::2] == 0)
        ))
    binding_stats['n_hcore_bridge'].append(
        np.sum((d_cont_per_hpatch[0::2] > 0) & (d_cont_per_hpatch[1::2] > 0)
               ))
    binding_stats['n_hcore_dangle'].append(
        n_hcore-binding_stats['n_hcore_free'][-1] -
        binding_stats['n_hcore_bridge'][-1]
        )
    single_patch_dir_contact = enforce_single_patch_dir_contact(direct_contact)
    cont_m_hpatch = generate_contact_matrix(single_patch_dir_contact)
    cont_m_hcore = np.zeros((n_mon, n_hcore))
    # Every two columns belongs to the same H-NS protein:
    cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
    # Symmetric squared monomer-monomer matrix:
    cont_m_m = np.matmul(cont_m_hcore, cont_m_hcore.T)  # symmetric
    cont_m_m_triu = np.triu(cont_m_m, 1)
    cont_m_m_nonzeros = np.array(np.where(cont_m_m_triu > 0)).T
    m_m_gen_dist = hns_genomic_distance(cont_m_m_nonzeros, topology, n_mon)
    binding_stats['n_hcore_cis'].append(np.count_nonzero(
        (m_m_gen_dist[:, 2] > 0) & (m_m_gen_dist[:, 2] <= cis_threshold)
        ))
    binding_stats['n_hcore_trans'].append(np.count_nonzero(
        m_m_gen_dist[:, 2] > cis_threshold))
    np.add.at(loop_length_hist, m_m_gen_dist[:, 2], 1)
    return binding_stats, loop_length_hist


def count_hns_clusters(
    direct_contact: np.ndarray,
) -> Tuple[Dict[str, List[int]], np.ndarray]:
    """Calculate the clustering statistics of H-NS proteins along the polymer
    chain.

    The H-NS-core-H-nS-core contact symmetric squared matrix is created from the asymmetric unsquared distance marix monomer-H-NS-patch matrix.

    Parameters
    ----------
    direct_contact: np.ndarray, shape (n_mon, n_hpatch)
        A binary matrix where each element (i, j) is 1 if atoms i and j have a
        contact, and 0 otherwise.

    topology: str
        Topology of the polymer containing monomers.

    cis_threshold: int, default 4
        The genomic distance in number of bonds (not number of monomers) less
        than or equal (inclusive) to which a H-NS is of cis-binding type.

    binding_stats: dict
        A dictionary of binding statistics of H-NS proteins to which new
        statisitics is appended.

    loop_length_hist: list of array
        A list of arrays where each array is of shape (n_bridged_pairs, 3).
        Each row of the an array contaion the index of first monomer, index of
        second one and the genomic distance between them.

    Return
    ------
    binding_stats: dict
        A new or updated dictionary contains the statistical data about the
        binding of H-NS proteins to monomers.

    loop_size_hist: np.ndarray
        A new or updated arrat contains the histogram of loop sizes.
    """
    #
    n_mon, n_hpatch = direct_contact.shape
    n_hcore = n_hpatch // 2
    single_patch_dir_contact = enforce_single_patch_dir_contact(direct_contact)
    cont_m_hpatch = generate_contact_matrix(single_patch_dir_contact)
    # Asymmetric unsquared monomer-H-NS-core matrix:
    # Every two columns belongs to the same H-NS protein
    cont_m_hcore = np.zeros((n_mon, n_hcore))
    cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
    # Symmetric squared H-NS-core-H-NS-core matrix:
    cont_hcore_hcore = np.matmul(cont_m_hcore.T, cont_m_hcore)
    return cont_hcore_hcore
