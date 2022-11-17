"""A sub-module to detect and analyze particle clusters in a given system. It
worth noting there is a difference analyzing particle clusters and
trajectory (or ensemble) clusters.
Below, different algorithms found in literature are implemented in scientific
Python.
"""
from typing import Dict, List, Tuple, Optional
from itertools import combinations
import numpy as np
import pandas as pd
import numpy.linalg as npla
from ..manage.typer import TransFociT


def apply_pbc(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: dict
) -> Tuple[np.ndarray, np.ndarray]:
    """Updates the periodic boundary condition (PBC) information in the
    dimentions (keys) by the length (values) passed by `pbc`. The key
    are indexes of the `pbc_lengths` and `pbc_lengths_inverse` arrays.

    Parameters
    ----------
    pbc_lengths: numpy.ndarray
        An array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        An array of the invesres of lengths/sizes in different directions.
    pbc: dict
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Return
    ------
    pbc_lengths: numpy.ndarray
        The updated array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        The updated array of the invesres of lengths/sizes in different
        directions.
    """
    for dim, length in pbc.items():
        pbc_lengths[dim] = length
        pbc_lengths_inverse[dim] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def dist_sq_matrix(
    positions: np.ndarray,
    pbc: Optional[dict] = None
) -> np.ndarray:
    """Finds the squares of distances between all the pair of atoms
    (even self-distances where i=j) for a given list of atoms' `positions`.

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different dimentions
        or components.

    pbc: dict, defualt None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Rutern
    ------
    dist_sq: np.ndarray
        A square matrix of size n_atoms * n_atoms in which each elements are
        the squares of the center-to-center distances between different pairs
        of atoms. All the diagonal elements of this materix is 0.

    References
    ----------
    2019 - Albanie S - Euclidean Distance Matrix Trick
    """
    n_atoms, n_dims = positions.shape
    # Defining pbc lengths with value zero in all directions:
    # The approach allows us to combine linear agebra and numpy broadcasting
    # and efficently apply pbc.
    pbc_lengths = np.zeros(n_dims)
    pbc_lengths_inverse = np.zeros(n_dims)
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc(
            pbc_lengths, pbc_lengths_inverse, pbc
        )
    # Redefining the shapes of different arrays to combine linear algebra
    # with numpy broadcasting.
    pbc_lengths = np.reshape(pbc_lengths, (n_dims, 1, 1))
    pbc_lengths_inverse = np.reshape(pbc_lengths_inverse, (n_dims, 1, 1))
    pos_T = positions.T
    pos_j_element = np.reshape(pos_T, (n_dims, n_atoms, 1))
    pos_i_element = np.reshape(pos_T, (n_dims, 1, n_atoms))
    # differences:
    dist_sq = pos_j_element - pos_i_element
    # applying pbc
    dist_sq = dist_sq - pbc_lengths * np.around(pbc_lengths_inverse * dist_sq)
    # squaring each elements
    dist_sq = dist_sq ** 2
    # sum over qaxis=0 means sum over d in x_ijd^2 where d is the number of
    # dimentions and x_ijd are d elements of a given r_ij^2.
    dist_sq = np.sum(dist_sq, axis=0)
    return dist_sq


def dist_sq_matrix_opt(
    positions: np.ndarray,
    pbc: Optional[dict] = None
) -> np.ndarray:
    """Finds the squares of distances between all the pair of atoms
    (even self-distances where i=j) for a given list of atoms' `positions`.
    In 3 dimesions, this method is not much faster the `dist_sq_matrix`.

    To-do
    -----
    Check the pbc applied below; is it correct?

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different dimentions
        or components.

    pbc: dict, defualt None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Rutern
    ------
    dist_sq: np.ndarray
        A square matrix of size n_atoms * n_atoms in which each elements are
        the squares of the center-to-center distances between different pairs
        of atoms. All the diagonal elements of this materix is 0.

    References
    ----------
    2019 - Albanie S - Euclidean Distance Matrix Trick
    """
    n_atoms, n_dims = positions.shape
    # Defining pbc lengths with value zero in all directions:
    # The approach allows us to combine linear agebra and numpy broadcasting
    # and efficently apply pbc.
    pbc_lengths = np.zeros(n_dims)
    pbc_lengths_inverse = np.zeros(n_dims)
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc(
            pbc_lengths, pbc_lengths_inverse, pbc
        )
    pbc_lengths = pbc_lengths**2
    pbc_lengths_inverse = pbc_lengths_inverse**2
    # Redefining the shapes of different arrays to combine linear algebra
    # with numpy broadcasting.
    gram = np.matmul(positions, positions.T)
    dist_sq = np.diag(gram).reshape(n_atoms, 1) + \
        np.diag(gram).reshape(1, n_atoms) - 2 * gram
    dist_sq = dist_sq - pbc_lengths * np.around(pbc_lengths_inverse * dist_sq)
    return dist_sq


def direct_contact(dist_sq: np.ndarray, cutoff: float) -> np.ndarray:
    """Finds direct contancts between pairs of atoms.

    Direct contact refers to a conditon in which the center-to-center
    distance between two atoms is equal to or less than the `cutoff` distance.
    This definition means that an atom is in direct contact with itself.
    Such a definition is correct since clusters with size 1 exist (i. e. a
    signle atom that is not clustered with other atoms); moreover, this
    definition esnures there is a correct number of atoms in a cluster.

    Parameters
    ----------
    dist_sq: np.ndarray
        A square matrix of size n_atoms * n_atoms in which each elements are
        the squares of the center-to-center distances between different pairs
        of atoms. All the diagonal elements of this materix is 0.

    cutoff: float
        The cutoff center-to-center distance for defining a direct
        contact.

    Return
    ------
    np.ndarray:
        A symmeteric matrix of size n_atoms * n_atoms in which each element is
        either 0 or 1 and indicates whether an direct contact is between a
        pair of atoms or not. 1 means direct contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition given above.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    return np.asarray(dist_sq <= (cutoff * cutoff), dtype=int)


def contact_matrix(dir_contacts: np.ndarray) -> np.ndarray:
    """Finds the direct and indirect contants between atoms in a system
    from the direct contacts between these atoms `dir_contacts`.

    Parameters
    ----------
    dir_contacts: np.ndarray
        A symmeteric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether an direct contact is
        between a pair of atoms or not. 1 means direct contact exists
        while 0 means it does not. All the diagonal elements of this matrix
        are 1 by the definition given above.

    Return
    ------
    contacts: np.ndarray
        A block diagonal matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a pair
        of atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the definition
        in `direct_contacts` function.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = dir_contacts.shape
    col_pairs = list(combinations(range(n_atoms), 2))
    contacts = dir_contacts.copy()
    # TODO Is this approach good or converting dir_contacts to
    # contact_matrix?
    for (i, j) in col_pairs:
        if np.any(contacts[:, i] & contacts[:, j]):
            transient = contacts[:, i] | contacts[:, j]
            contacts[:, i] = transient
            contacts[:, j] = transient
    return contacts


def contact_from_sq_dists(
    positions: np.ndarray,
    cutoff: float,
    pbc: Optional[dict] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Creates the direct contact matrix for a group of atoms using a
    cutoff for the center-to-center distance as the ceritrion for clustering.

    All the diagonal elemnts of the contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many
    repeated columns.

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different
        dimentions or components.

    cutoff: float
        The cutoff center-to-center distance for defining a direct
        contact.

    pbc: dict, default None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Return
    ------
    dir_contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    contacts: np.ndarray
        A matrix of size n_atoms * n_atoms in which each element is either
        0 or 1 and indicates whether a contact is between a pair of
        atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    foci_pair_dist: np.ndarray
        An upper triangluar matrix with 0 zero diagonal elements that
        contains the center-to-center distance between all the pairs of
        atoms.
    """
    dist_sq_mat = dist_sq_matrix(positions, pbc)
    dir_contacts = direct_contact(dist_sq_mat, cutoff)
    contacts = contact_matrix(dir_contacts)
    return dir_contacts, contacts


def random_direct_contacts(
    n_atoms: int,
    seed: Optional[int] = 0,
    save_to: Optional[Tuple[str, str]] = None
) -> np.ndarray:
    """Generates a direct-contact matrix with size `n_atoms` which can be
    used to test a clustering algorithm. See `direct_contact`function
    description for more info.

    All the diagonal elemnts of the direct contact matrix are 1; moreover,
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
        A tumple of two strings in which the first element is an/a
        absolute/relative path of a directory to which the random direct
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
    n_atoms = 5
    # Creating direct-contacts matrix
    np.random(seed)
    dir_contacts = np.random.binomial(n=1, p=0.35, size=[n_atoms * n_atoms])
    dir_contacts = np.reshape(dir_contacts, (n_atoms, n_atoms))
    dir_contacts_up = np.triu(dir_contacts)  # upper triangle
    dir_contacts_diag = np.identity(n_atoms, dtype=int)
    # dropping diagnal elements:
    dir_contacts_up = dir_contacts_up - dir_contacts_diag
    # ensuring all diagonal elements are either 0 or 1:
    dir_contacts_up[dir_contacts_up < 0] = 0
    dir_contacts_low = dir_contacts_up.T  # lower traingle: symmetry!
    dir_contacts = dir_contacts_up + dir_contacts_low + dir_contacts_diag
    # Creating the full direct contacts matrix
    dir_contacts = contact_matrix(dir_contacts)
    if save_to is not None:
        output = save_to[0] + \
            f"random_direct_contact-natoms{n_atoms}" + save_to[1]
        np.save(output, dir_contacts)
    return dir_contacts


def count_bonds(dir_contacts: np.ndarray):
    """Infers the number and size of bonds from the direct contact matrix
    `dir_contacts`.

    The direct contact matrix is symmetric and its diagonal elements are 1.
    The column- or row-wise sum of the matrix shows the number of bonds
    an atom has with itself and other atoms. A self-bond happens due to the
    definition of the direct contact matirx (in which the diagonal elements
    are 1). Such self-bonds are substracted below to give the correct number of
    bonds, thus number of bonds per atom ranges from 0 to $n_{atoms}-1$.

    Parameters
    ----------
    contacts: np.ndarray
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


def count_clusters(contacts: np.ndarray) -> np.ndarray:
    """Infers the number and size of clusters from a contact matrix `contacts`.

    All the diagonal elemnts of the direct contact matrix are 1; moreover,
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
    # While Tte eigenvalue 0 is not physically important, the eigenvaue n_atoms
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
    nmon_small: int,
) -> np.ndarray:
    """Generates a dictionany of the information about all the possible pairs
    of large monomers in a circular/ring heterogenous polymner composed of two
    different monomer types (small and large monomers).

    This information is an array of shape(n_pairs,3) where the first two
    columns are the pair index along the backbone (a.k.a genomic position), the
    2nd one is the genomic distance (i.e. the minimum number of genomic indices
    between a pair in a circular polymer).

    Parameters
    ----------
    genomic_pos: np.ndaaray
        A vector of size `nmon_large` that contains the genomic position of
        large monomers
    nmon_large: int
        Number of large monomers
    nmon_large: int
        Number of small monomers

    Return
    ------
    pair_info: np.ndarray
        An array of shape(n_pairs,3) that contains the genomic indices, and
        genomic distance of each pair.
    """
    nmon = nmon_small + nmon_large
    pairs = np.array(list(combinations(range(nmon_large), 2)), dtype=np.int)
    genomic_idx = np.choose(pairs, genomic_pos)
    genomic_dist = np.diff(genomic_idx, axis=1)
    # Applying the topological contrain on the genmoic distance, i. e. finding
    # minimum distance between a pair on a circle.
    # the pair is inclusive, so their difference is one unit more than the
    # actual number of genomic units between them, thus sabtracting by 1.0.
    genomic_dist = np.minimum.reduce((nmon-genomic_dist, genomic_dist)) - 1.0
    pair_info = np.concatenate(
        (genomic_idx, genomic_dist),
        axis=1
        )
    return pair_info


def v_histogram(
    data: np.ndarray,
    bins: int,
    range: Tuple[float, float],
    **kwargs
) -> np.ndarray:
    """Vectorizes `n.histogram` so it can be mapped to a collection of data
    with different `nbins` and `range` kwargs. The `range` is divided into
    `nbins` bins, so there are `nbins` equal bins.

    Parameters
    ----------
    data: np.ndarray
        1D array-like input data
    bins: int
        Number of bins
    range: float
        The lower and upper range of bins.

    Return
    ------
    hist: np.ndarray
        The frequencies of `data` in bins.
    """
    hist, _ = np.histogram(data, bins=bins, range=range, **kwargs)
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
    hist: np.ndaaray
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
) -> Dict[str, np.ndarray]:
    """Generates histograms of pair distances, use the same number of bins and
    the same bin range for all pairs.

    The lower bin range is 0 while the upper one is the maximum pf the maximum
    physical distances between any pairs.

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
    upper_range = np.ceil(pairs_dist.max() / binsize) * binsize
    max_range = (0, upper_range)
    max_bins = np.ceil((max_range[1] - max_range[0]) / binsize).astype(np.int)
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


def whole_distMat_foci(
    whole_path: str,
    whole_info: TransFociT
) -> Tuple[
    Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]
]:
    """Generates time series and histograms of foci distances from distance
    matrix of (large) monomers.

    A foci is a pair of large monomer. The ditance matrix is upper-triangular
    matrix containing all the spatial distances between any possible foci of
    large monomers.

    Parameters
    ----------
    whole_path: str
        Path to a "whole" diatance matrix.
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
    piars_tseries: dict
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
