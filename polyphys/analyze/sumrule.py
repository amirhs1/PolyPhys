"""
==========================================================
:mod:`polyphys.analyze.sumrule`
==========================================================


The :mod:`polyphys.analyze.sumrule` module tests the density sum rule of
Jeon et al. (2016) on the ensemble-averaged density profiles of a space: it
scales each species' profile by its particle size, normalizes it by a
reference value, folds a longitudinal profile about the polymer's centre,
and sums the species' scaled profiles state point by state point. The
profiles come from the ``ensAvg`` stage of the artifact lineage that
:mod:`polyphys.manage.organizer` builds; the per-frame densities behind them
are the caller's work.

The profiles must be the raw ensemble averages, in their physical units: a
volume fraction or a number density, not a profile already divided by its
sum, its integral, or its maximum. The sum rule relates the species to each
other through their physical scales, so a species normalized on its own
enters the sum :math:`\\sum_i \\phi_i / a_i` with an arbitrary weight and
``<prop>-scale`` loses its unit and its comparability with the nominal bulk
value. Only the ``<prop>-norm`` columns of such a profile survive, because
dividing by the profile's own ``'first'`` or ``'max'`` value is
scale-invariant. :func:`space_profiles` warns when a profile sums to one,
the signature of a profile normalized by its sum, but cannot detect the
other normalizations from one profile.

Notation
========

:math:`a_i`
    The particle diameter of species :math:`i`, in the monomer diameter: the
    parser attribute that the study's ``species_sizes`` names for it.
:math:`\\phi_i(x)`, :math:`\\rho_i(x)`
    The local volume fraction (dimensionless) and local number density (per
    cubic monomer diameter) of species :math:`i`, related by
    :math:`\\phi_i = \\pi a_i^3 \\rho_i / 6`.
:math:`r`
    The radial coordinate: the distance from the cylinder axis in a
    cylinder and from the polymer's centre of mass in a periodic cube. The
    bins of an ``r`` profile run outward from :math:`r = 0`.
:math:`z`
    The longitudinal coordinate: the axial distance from the polymer's
    centre of mass in a cylinder, with bins symmetric about :math:`z = 0`.
*scaled profile*
    One term of the sum rule, the column ``<prop>-scale``:
    :math:`\\phi_i / a_i` for a volume fraction and :math:`a_i^2 \\rho_i` for
    a number density.
*normalizer*
    The reference value a scaled profile, or the species sum, is divided by:
    see :py:data:`~polyphys.manage.types.NormalizerT`.
:math:`D`, :math:`a_c`
    The cylinder diameter and a crowder diameter; the ratio
    :math:`\\kappa = (D - a_c) / a_c` sets how far the crowder layering on
    the wall reaches inward.

Contents
========

- :func:`normalize_r` and :func:`normalize_z`: normalize one scaled profile,
  folding a longitudinal one about :math:`z = 0`.
- :func:`space_profiles`: scale, normalize, and annotate every
  ensemble-averaged profile of one species across a space.
- :func:`sum_rule_profiles`: sum the species' scaled profiles and normalize
  the sum.
- :data:`PROFILE_NORMALIZERS` and :data:`SUM_NORMALIZERS`: the default
  normalizers by geometry and direction.

Dependencies
============
- `numpy` and `pandas`.
- :mod:`polyphys.manage.organizer`:
  :func:`~polyphys.manage.organizer.make_database` and
  :func:`~polyphys.manage.organizer.create_fullname`.
- :mod:`polyphys.manage.utils` and :mod:`polyphys.manage.types`.

References
==========
- C. Jeon, C. Hyeon, Y. Jung, and B.-Y. Ha, "How are molecular crowding and
  the spatial organization of a biopolymer interrelated", Soft Matter 12,
  9786-9796 (2016). https://doi.org/10.1039/c6sm01924b
- J. Kim, C. Jeon, H. Jeong, Y. Jung, and B.-Y. Ha, "A polymer in a crowded
  and confined space: effects of crowder size and poly-dispersity", Soft
  Matter 11, 1877-1888 (2015). https://doi.org/10.1039/c4sm02198c
- A. H. Sadeghi Isfahani, PhD thesis, University of Waterloo (2024),
  Chapter 3.

.. versionadded:: 0.6.0
"""

import logging
import warnings
from glob import glob
from typing import Any, Callable, Literal, Mapping, Sequence, TYPE_CHECKING

import numpy as np
import pandas as pd

from ..manage.organizer import create_fullname, make_database
from ..manage.types import (
    DirectionT,
    EntityT,
    GeometryT,
    GroupT,
    NormalizerT,
    ParserClassT,
)
from ..manage.utils import invalid_keyword, round_to_multiple, sort_filenames

if TYPE_CHECKING:
    # Unused here, but sphinx_autodoc_typehints runs this block to resolve
    # the "ParserBase" forward reference inside ParserClassT.
    from ..manage.parser import ParserBase  # noqa: F401

logger = logging.getLogger(__name__)


# The normalizer each density profile takes by default, keyed by (geometry,
# direction). In a cylinder, r runs from the axis to the wall, whose crowder
# layering makes the maximum the wrong reference, so the axis value is used;
# along z and in a periodic cube the maximum is the far-field value of a
# crowder profile and the centre value of a monomer profile.
PROFILE_NORMALIZERS: dict[tuple[GeometryT, DirectionT], NormalizerT] = {
    ('cylindrical', 'r'): 'first',
    ('cylindrical', 'z'): 'max',
    ('cubic', 'r'): 'max',
}
"""The default normalizer of a density profile, by (geometry, direction):
``'first'`` for ``r`` in a cylinder, ``'max'`` otherwise.

.. versionadded:: 0.6.0
"""

# The normalizer the species sum takes by default: its own value at the
# axis in a cylinder, and its far-field value along z and in a cube, where
# only crowders remain.
SUM_NORMALIZERS: dict[tuple[GeometryT, DirectionT], NormalizerT] = {
    ('cylindrical', 'r'): 'first',
    ('cylindrical', 'z'): 'last',
    ('cubic', 'r'): 'last',
}
"""The default normalizer of the species sum, by (geometry, direction):
``'first'`` for ``r`` in a cylinder, ``'last'`` otherwise.

.. versionadded:: 0.6.0
"""


def _normalizer_value(
    prop: str, profile: pd.DataFrame, normalizer: NormalizerT
) -> float:
    """
    Return the value that `normalizer` picks from the ``<prop>-scale``
    column of `profile`: its first or last value, its maximum, the given
    number, or the result of the given function of `profile`.

    Raises
    ------
    ValueError
        If `normalizer` is a string other than ``'first'``, ``'last'``, or
        ``'max'``, or neither a string, a number, nor callable.
    """
    scale = profile[f"{prop}-scale"]
    if isinstance(normalizer, str):
        invalid_keyword(normalizer, ["first", "last", "max"])
        if normalizer == "first":
            return float(scale.iloc[0])
        if normalizer == "last":
            return float(scale.iloc[-1])
        return float(scale.max())
    if isinstance(normalizer, (int, float)) and not isinstance(
        normalizer, bool
    ):
        return float(normalizer)
    if callable(normalizer):
        return float(normalizer(profile))
    raise ValueError(
        f"'normalizer' must be 'first', 'last', 'max', a number, or a "
        f"function of the profile, not {normalizer!r}."
    )


def _with_parser(
    func: Callable[[pd.DataFrame, Any], float], info: Any
) -> Callable[[pd.DataFrame], float]:
    """
    Bind the parser instance `info` as the second argument of `func`, a
    two-argument normalizer of :func:`space_profiles`.
    """
    def normalizer(profile: pd.DataFrame) -> float:
        return func(profile, info)
    return normalizer


def _normalize(
    prop: str, profile: pd.DataFrame, normalizer: NormalizerT
) -> pd.DataFrame:
    """
    Return a copy of `profile` with the ``<prop>-normalizer`` column set to
    the value `normalizer` picks and ``<prop>-norm`` to ``<prop>-scale``
    divided by it; a zero value warns and sets ``<prop>-norm`` to zero.
    """
    out = profile.copy()
    value = _normalizer_value(prop, out, normalizer)
    out[f"{prop}-normalizer"] = value
    if value != 0:
        out[f"{prop}-norm"] = out[f"{prop}-scale"] / value
    else:
        warnings.warn(
            f"The normalizer of '{prop}' is zero; its normalized values are"
            " set to zero.",
            UserWarning,
        )
        out[f"{prop}-norm"] = 0.0
    return out


def normalize_r(
    prop: str,
    ens_avg: pd.DataFrame,
    normalizer: NormalizerT = "first",
) -> pd.DataFrame:
    """
    Normalize the scaled radial profile ``<prop>-scale`` of an
    ensemble-averaged distribution by a reference value.

    The radial coordinate :math:`r` is the distance from the cylinder axis
    in a cylinder and from the polymer's centre of mass in a periodic cube,
    in the monomer diameter; the profile's bins run outward from
    :math:`r = 0`.

    Parameters
    ----------
    prop : str
        The full property name, such as ``'rPhiCrd'``; the input has the
        column ``<prop>-scale``.
    ens_avg : pd.DataFrame
        The ensemble-averaged distribution, one row per bin. It is not
        modified.
    normalizer : :py:data:`~polyphys.manage.types.NormalizerT`, \
default 'first'
        The reference value: ``'first'`` is the value at :math:`r = 0`,
        ``'last'`` the value of the outermost bin, ``'max'`` the maximum;
        a number is used as given; a function is called with `ens_avg` and
        returns the value.

    Returns
    -------
    pd.DataFrame
        A copy of `ens_avg` with two more columns: ``<prop>-normalizer``,
        the reference value, and ``<prop>-norm``, ``<prop>-scale`` divided
        by it.

    Raises
    ------
    ValueError
        If `normalizer` is none of the accepted forms.

    Warns
    -----
    UserWarning
        If the reference value is zero; ``<prop>-norm`` is then zero.

    Examples
    --------
    >>> profile = pd.DataFrame({'bin_center': [0.5, 1.5, 2.5],
    ...                         'rPhiCrd-scale': [2.0, 4.0, 1.0]})
    >>> normalize_r('rPhiCrd', profile)['rPhiCrd-norm'].tolist()
    [1.0, 2.0, 0.5]
    >>> normalize_r('rPhiCrd', profile, 'max')['rPhiCrd-norm'].tolist()
    [0.5, 1.0, 0.25]
    >>> normalize_r('rPhiCrd', profile, 8.0)['rPhiCrd-norm'].tolist()
    [0.25, 0.5, 0.125]

    A function sees the whole frame, so it can average the inner bins:

    >>> normalize_r('rPhiCrd', profile,
    ...             lambda df: df['rPhiCrd-scale'].iloc[:2].mean()
    ...             )['rPhiCrd-normalizer'].iloc[0]
    np.float64(3.0)

    Notes
    -----
    Time and space are :math:`O(n_\\text{bins} \\, n_\\text{cols})` for the
    copy. The references of the sum rule are in :func:`space_profiles`.

    .. versionadded:: 0.6.0
    """
    return _normalize(prop, ens_avg, normalizer)


def normalize_z(
    prop: str,
    ens_avg: pd.DataFrame,
    normalizer: NormalizerT = "max",
    fold: bool = True,
) -> pd.DataFrame:
    """
    Normalize the scaled longitudinal profile ``<prop>-scale`` of an
    ensemble-averaged distribution by a reference value and, by default,
    fold it about :math:`z = 0`.

    The longitudinal coordinate :math:`z` is the axial distance from the
    polymer's centre of mass, in the monomer diameter. The bins are assumed
    symmetric about :math:`z = 0`, so folding averages bin :math:`i` with
    its mirror bin :math:`n - 1 - i` and keeps the centre bin of an odd
    count as it is; the folded ``bin_center`` is :math:`|z|`.

    Parameters
    ----------
    prop : str
        The full property name, such as ``'zRhoMon'``; the input has the
        column ``<prop>-scale``.
    ens_avg : pd.DataFrame
        The ensemble-averaged distribution, one row per bin, with a
        ``bin_center`` column when `fold` is True. It is not modified.
    normalizer : :py:data:`~polyphys.manage.types.NormalizerT`, default 'max'
        The reference value, as in :func:`normalize_r`; ``'max'`` is the
        far-field value of a crowder profile and the centre value of a
        monomer profile when each is monotone in :math:`|z|`.
    fold : bool, default True
        Whether to fold the profile about :math:`z = 0`.

    Returns
    -------
    pd.DataFrame
        A copy of `ens_avg` with ``<prop>-normalizer`` and ``<prop>-norm``
        added, as in :func:`normalize_r`. With `fold`, every numeric column
        is averaged with its mirror bin, the other columns are those of the
        :math:`z \\geq 0` half, and the result has ``ceil(n / 2)`` rows
        ordered by :math:`|z|`.

    Raises
    ------
    ValueError
        If `normalizer` is none of the accepted forms, or `fold` is True
        and `ens_avg` has no ``bin_center`` column.

    Warns
    -----
    UserWarning
        If the reference value is zero; ``<prop>-norm`` is then zero.

    Examples
    --------
    An even number of bins folds to half as many:

    >>> profile = pd.DataFrame({'bin_center': [-1.5, -0.5, 0.5, 1.5],
    ...                         'zRhoMon-scale': [1.0, 3.0, 3.0, 1.0]})
    >>> normalize_z('zRhoMon', profile)
       bin_center  zRhoMon-scale  zRhoMon-normalizer  zRhoMon-norm
    0         0.5            3.0                 3.0      1.000000
    1         1.5            1.0                 3.0      0.333333

    An odd number keeps its centre bin:

    >>> profile = pd.DataFrame({'bin_center': [-1.0, 0.0, 1.0],
    ...                         'zRhoMon-scale': [2.0, 4.0, 2.0]})
    >>> normalize_z('zRhoMon', profile)['zRhoMon-norm'].tolist()
    [1.0, 0.5]

    Notes
    -----
    The ``-var`` and ``-sem`` columns of the input are averaged like every
    other numeric column, which is what the research code did; the result
    is not the uncertainty of the folded mean, which would need the
    per-replicate folded profiles. Time and space are
    :math:`O(n_\\text{bins} \\, n_\\text{cols})`.

    .. versionadded:: 0.6.0
    """
    out = _normalize(prop, ens_avg, normalizer)
    if not fold:
        return out
    if "bin_center" not in out.columns:
        raise ValueError("Folding needs a 'bin_center' column.")
    n_bins = len(out)
    half = n_bins // 2
    # The z > 0 half in increasing |z|, and the z < 0 half reversed to match
    # it bin for bin; the centre bin of an odd count stays on its own.
    pos = out.iloc[n_bins - half:].reset_index(drop=True)
    neg = out.iloc[:half].iloc[::-1].reset_index(drop=True).copy()
    neg["bin_center"] = -neg["bin_center"]
    numeric = out.select_dtypes(include="number").columns
    folded = pos.copy()
    folded[numeric] = 0.5 * (pos[numeric] + neg[numeric])
    if n_bins % 2:
        centre = out.iloc[[half]].copy()
        centre["bin_center"] = centre["bin_center"].abs()
        folded = pd.concat([centre, folded], axis=0)
    return folded.reset_index(drop=True)


def space_profiles(
    input_database: str,
    quantity: Literal["Rho", "Phi"],
    direction: DirectionT,
    species: EntityT,
    parser: ParserClassT,
    hierarchy: str,
    physical_attrs: list[str],
    group: GroupT,
    species_sizes: Mapping[EntityT, str],
    *,
    crowders: Sequence[EntityT] = ("Crd",),
    normalizer: (
        Literal["first", "last", "max"] | float
        | Callable[[pd.DataFrame, Any], float] | None
    ) = None,
    divisor: float = 0.025,
    round_to: int = 3,
    is_save: bool = False,
    crowder_attr: str = "phi_bulk_c",
) -> pd.DataFrame:
    """
    Scale and normalize the ensemble-averaged density profiles of one
    species across a space for the density sum rule, add the physical
    attributes of each state point, and concatenate them into one
    DataFrame.

    The density sum rule of Jeon et al. (2016), eqns (5) to (7), states that
    in a crowded polymer solution the volume fractions :math:`\\phi_i(x)`
    of the species :math:`i`, each rescaled by its particle diameter
    :math:`a_i`, sum to the far-field value of the crowder term:

    .. math::

        \\sum_i \\frac{\\phi_i(x)}{a_i} \\approx \\frac{\\phi_c(\\infty)}{a_c},
        \\qquad
        \\sum_i a_i^2 \\rho_i(x) \\approx a_c^2 \\rho_c(\\infty),

    where the number density :math:`\\rho_i = 6 \\phi_i / (\\pi a_i^3)`
    makes the two forms the same statement, and :math:`x` is the radial
    coordinate :math:`r` (the distance from the cylinder axis in a cylinder
    and from the polymer's centre of mass in a periodic cube) or the
    longitudinal coordinate :math:`z` (the axial distance from the
    polymer's centre of mass). This function builds one term of the sum:
    for ``quantity='Phi'`` the scaled profile ``<prop>-scale`` is the
    volume fraction divided by the species' diameter, and for
    ``quantity='Rho'`` the number density times the diameter squared.
    :func:`sum_rule_profiles` builds the sum.

    The input profiles must be the raw ensemble averages, in their physical
    units. A profile normalized beforehand, by its sum, its integral, or its
    maximum, enters the species sum with an arbitrary weight, and its
    ``<prop>-scale`` loses its unit; only its ``<prop>-norm`` survives,
    since dividing by the profile's own reference value is scale-invariant.
    A profile whose ``<prop>-mean`` sums to one, as one normalized by its
    sum does, raises a :class:`UserWarning` naming the file; the other
    normalizations cannot be told from one profile.

    Parameters
    ----------
    input_database : str
        Path to the directory of the ensemble-averaged profiles, named
        ``<ensemble_long>-<group>-<direction><quantity><species>-ensAvg.csv``,
        ending with a path separator. Each file has a ``bin_center`` column
        and a ``<direction><quantity><species>-mean`` column holding the
        raw ensemble-averaged volume fraction or number density.
    quantity : {'Rho', 'Phi'}
        The number density or the volume fraction.
    direction : {'r', 'z'}
        The coordinate of the profile; ``'z'`` exists in a cylinder only.
    species : EntityT
        The species label of the profile, such as ``'Mon'`` or ``'Crd'``; a
        key of `species_sizes`.
    parser : ParserClassT
        The study's parser class; its geometry chooses the default
        normalizer and the recentred columns.
    hierarchy : str
        Pattern prefix of the filenames, such as ``'N*'``.
    physical_attrs : list of str
        The parser attributes added as columns to the result.
    group : GroupT
        The particle group of the profiles, usually ``'all'``.
    species_sizes : mapping of (EntityT -> str)
        Each species label mapped to the parser attribute holding its
        particle diameter, in the monomer diameter: the ``species_sizes`` of
        the study's :class:`~polyphys.manage.projects.ProjectSpec`.
    crowders : sequence of EntityT, default ('Crd',)
        The crowder species, each a key of `species_sizes`; one set of
        ``bin_center`` columns in crowder units is added per crowder.
    normalizer : NormalizerT or callable or None, default None
        The reference value of the normalization. None takes the default of
        :data:`PROFILE_NORMALIZERS` for the parser's geometry and
        `direction`. Otherwise ``'first'``, ``'last'``, ``'max'``, or a
        number, as in :func:`normalize_r`; or a function of two arguments,
        the scaled profile and the parser instance of the file, returning
        the value, so that ``lambda df, info: info.phi_bulk_c / info.dcrowd``
        normalizes a ``'Phi'`` profile by the nominal bulk volume fraction
        of crowders in scale units. An explicit value is logged at INFO
        with the default it replaces.
    divisor : float, default 0.025
        Rounding step of the bulk volume fraction of crowders named by
        `crowder_attr`.
    round_to : int, default 3
        Number of decimal places of that rounding.
    is_save : bool, default False
        If True, save the result as
        ``<space>-<group>-<direction><quantity><species>-normalizedRescaled-space.csv``
        in the ``analysis`` phase's ``space`` directory, which
        :func:`make_database` creates next to `input_database`, and log the
        path at INFO.
    crowder_attr : str, default 'phi_bulk_c'
        The parser attribute holding the bulk volume fraction of crowders,
        added as a column of that name and rounded into
        ``phi_bulk_c_round``.

    Returns
    -------
    pd.DataFrame
        The profiles of every file, concatenated in filename order, with
        ``<prop> = <direction><quantity><species>`` and these columns added
        to each file's own:

        - ``<prop>-scaler``: the diameter :math:`a` for ``'Phi'`` and
          :math:`a^2` for ``'Rho'``; ``<prop>-sumrule_constant`` repeats
          it.
        - ``<prop>-scale``: ``<prop>-mean`` divided by :math:`a` for
          ``'Phi'`` and times :math:`a^2` for ``'Rho'``.
        - ``<prop>-scale-normalized_curve``: ``<prop>-scale`` divided by
          its sum over the bins.
        - ``<prop>-normalizer`` and ``<prop>-norm``: from
          :func:`normalize_r` or, for ``'z'``, :func:`normalize_z` with
          folding, so that a ``'z'`` profile has ``ceil(n / 2)`` rows.
        - ``<prop>-norm-sem``, when the file has it and `normalizer` is
          None: the error of ``<prop>-norm``, carried through from the
          file and, for ``'z'``, folded like every other numeric column.
          :mod:`polyphys.analyze.pipeline` writes it as the standard error
          over the wholes of each whole's profile divided by its own value
          in the bin the default normalizer names, since the file's
          ``-sem`` is the error of the mean, not of the ratio. With an
          explicit `normalizer` the column no longer describes
          ``<prop>-norm``, so it is dropped with a warning.
        - ``bin_center-norm``: ``bin_center`` divided by its maximum.
        - ``ensemble_long``, the file's state-point name; one column per
          entry of `physical_attrs`; `crowder_attr`; and
          ``phi_bulk_c_round``.
        - For each crowder of `crowders` with size attribute ``c``:
          ``bin_center-<c>``, which is :math:`2 b / a_c` for the bin centre
          :math:`b`, the bin centre in crowder radii. For ``'r'`` in a
          cylinder of diameter :math:`D` also
          ``bin_center-<c>-recentered``, which is
          :math:`2 b / a_c - (D \\bmod a_c) / a_c`, so that the wall sits at
          an integer number of crowder diameters and the layers of
          different state points overlay, and ``bin_center-recentered-norm``
          (``bin_center-recentered-norm-<c>`` with several crowders), which
          is :math:`b - (D \\bmod a_c)` divided by its maximum.

    Raises
    ------
    ValueError
        If `quantity` is not ``'Rho'`` or ``'Phi'``; if the parser's
        geometry and `direction` have no entry in
        :data:`PROFILE_NORMALIZERS`, such as ``'z'`` in a cube; if no file
        matches; or if `normalizer` is none of the accepted forms.
    KeyError
        If `species` or a crowder is not a key of `species_sizes`.

    Warns
    -----
    UserWarning
        If a file's ``<prop>-mean`` sums to one within ``1e-6``, which marks
        a profile normalized by its sum; the file is processed regardless,
        since a raw profile can sum to one by coincidence. If a file has a
        ``<prop>-norm-sem`` column and `normalizer` is not None; the column
        is dropped. Also if a file's reference value is zero; its
        ``<prop>-norm`` is then zero.

    Examples
    --------
    One ensemble-averaged radial volume fraction of crowders of diameter
    ``sig2.0`` in a cylinder, normalized by its axis value:

    >>> import tempfile
    >>> from pathlib import Path
    >>> from polyphys.manage.parser import SumRuleCyl
    >>> tmp = tempfile.TemporaryDirectory()
    >>> database = Path(tmp.name) / 'analysis' / 'N20-all-ensAvg'
    >>> database.mkdir(parents=True)
    >>> state = 'N20epsilon5.0r5.5lz25.0sig2.0nc100dt0.005bdump1000adump5000'
    >>> pd.DataFrame({'bin_center': [0.5, 1.5, 2.5],
    ...               'rPhiCrd-mean': [0.2, 0.3, 0.1]}).to_csv(
    ...     database / f'{state}-all-rPhiCrd-ensAvg.csv', index=False)
    >>> space = space_profiles(str(database) + '/', 'Phi', 'r', 'Crd',
    ...                        SumRuleCyl, 'N*', ['dcyl'], 'all',
    ...                        {'Mon': 'dmon', 'Crd': 'dcrowd'})
    >>> space[['bin_center', 'rPhiCrd-scale', 'rPhiCrd-norm',
    ...        'bin_center-dcrowd-recentered']]
       bin_center  rPhiCrd-scale  rPhiCrd-norm  bin_center-dcrowd-recentered
    0         0.5           0.10           1.0                           0.5
    1         1.5           0.15           1.5                           1.5
    2         2.5           0.05           0.5                           2.5
    >>> space['dcyl'].tolist()
    [10.0, 10.0, 10.0]
    >>> tmp.cleanup()

    Notes
    -----
    - Files are selected by globbing ``input_database + hierarchy +
      '-<group>-<prop>-ensAvg.csv'`` and sorted alphanumerically by name.
    - Let n_files be the number of matched files, n_bins the number of rows
      of each, and n_cols the number of columns of the result. Time is
      :math:`O(n_\\text{files} \\, n_\\text{bins} \\, n_\\text{cols} +
      n_\\text{files} \\log n_\\text{files})`; one file is in memory at a
      time while the finished frames accumulate, so peak memory is about
      twice the result, :math:`O(n_\\text{files} \\, n_\\text{bins} \\,
      n_\\text{cols})`.
    - The default normalizers follow the research code: the value at
      :math:`r = 0` in a cylinder, where the outermost bins are the crowder
      layers on the wall and no far field exists, and the maximum
      otherwise. The sum rule is a local relation and holds in the
      moderate-compaction regime to about 15 percent (Jeon et al., 2016,
      Section 3.1); under cylindrical confinement its radial form holds
      only where the wall layering is weak, which the ratio
      :math:`\\kappa = (D - a_c) / a_c` controls (Sadeghi Isfahani, 2024,
      Chapter 3).
    - ``<prop>-norm`` is the normalized ensemble average, the published
      estimator (Jeon et al., 2016, p. 9788: a time average within each
      run, then over runs, then normalized by the far-field value), not the
      average of per-whole normalized profiles. Its error is the
      ``<prop>-norm-sem`` of the file, from the spread of the independent
      wholes, each normalized at the same bin.

    References
    ----------
    - C. Jeon, C. Hyeon, Y. Jung, and B.-Y. Ha, "How are molecular crowding
      and the spatial organization of a biopolymer interrelated", Soft
      Matter 12, 9786-9796 (2016). https://doi.org/10.1039/c6sm01924b
    - J. Kim, C. Jeon, H. Jeong, Y. Jung, and B.-Y. Ha, "A polymer in a
      crowded and confined space: effects of crowder size and
      poly-dispersity", Soft Matter 11, 1877-1888 (2015).
      https://doi.org/10.1039/c4sm02198c
    - A. H. Sadeghi Isfahani, PhD thesis, University of Waterloo (2024),
      Chapter 3.

    .. versionadded:: 0.6.0
    .. versionchanged:: 0.7.0
       Carry ``<prop>-norm-sem`` through, and warn on a profile that sums
       to one.
    """
    invalid_keyword(quantity, ["Rho", "Phi"])
    if species not in species_sizes:
        raise KeyError(
            f"Species {species!r} is not in 'species_sizes' "
            f"{sorted(species_sizes)}."
        )
    missing = [label for label in crowders if label not in species_sizes]
    if missing:
        raise KeyError(
            f"Crowders {missing} are not in 'species_sizes' "
            f"{sorted(species_sizes)}."
        )
    geometry = parser._geometry
    try:
        default = PROFILE_NORMALIZERS[(geometry, direction)]
    except KeyError:
        raise ValueError(
            f"No density profile along {direction!r} in a {geometry} "
            f"geometry; the (geometry, direction) pairs are "
            f"{sorted(PROFILE_NORMALIZERS)}."
        ) from None
    if normalizer is not None:
        logger.info(
            "Normalizing %s%s%s profiles with %r instead of the %s %s "
            "default %r.", direction, quantity, species, normalizer, geometry,
            direction, default,
        )

    prop = f"{direction}{quantity}{species}"
    prop_ext = f"-{group}-{prop}-ensAvg.csv"
    ens_avg_csvs = sort_filenames(
        glob(input_database + hierarchy + prop_ext), ["csv"]
    )
    if not ens_avg_csvs:
        raise ValueError(
            f"No '{hierarchy}{prop_ext}' file in '{input_database}'."
        )
    size_attr = species_sizes[species]
    prop_csvs = []

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        prop_info = parser(ens_avg_csv[0], "ensemble_long", group, ispath=True)

        # A profile normalized by its sum adds up to one; the sum rule needs
        # the raw profile, so warn, but go on: a raw one can sum to one too.
        if abs(ens_avg[f"{prop}-mean"].sum() - 1.0) <= 1e-6:
            warnings.warn(
                f"The '{prop}-mean' column of '{ens_avg_csv[0]}' sums to one"
                " within 1e-6, as a profile normalized by its sum does; the"
                " sum rule needs the raw ensemble-averaged profile.",
                UserWarning,
            )

        # The pipeline's per-whole error belongs to the default normalizer;
        # with another, it would not describe <prop>-norm.
        norm_sem = f"{prop}-norm-sem"
        if norm_sem in ens_avg.columns and normalizer is not None:
            warnings.warn(
                f"Dropping '{norm_sem}' of '{ens_avg_csv[0]}': it is the"
                f" error under the default normalizer {default!r}, not under"
                f" {normalizer!r}.",
                UserWarning,
            )
            ens_avg = ens_avg.drop(columns=norm_sem)

        # One term of the sum rule: phi / a, or rho * a**2
        size = getattr(prop_info, size_attr)
        if quantity == "Phi":
            scaler = size
            ens_avg[f"{prop}-scale"] = ens_avg[f"{prop}-mean"] / scaler
        else:
            scaler = size ** 2
            ens_avg[f"{prop}-scale"] = ens_avg[f"{prop}-mean"] * scaler
        ens_avg[f"{prop}-scaler"] = scaler
        ens_avg[f"{prop}-scale-normalized_curve"] = (
            ens_avg[f"{prop}-scale"] / ens_avg[f"{prop}-scale"].sum()
        )

        # Normalize by the reference value; fold a z profile about z = 0
        file_normalizer: NormalizerT
        if normalizer is None:
            file_normalizer = default
        elif callable(normalizer):
            file_normalizer = _with_parser(normalizer, prop_info)
        else:
            file_normalizer = normalizer
        if direction == "z":
            ens_avg = normalize_z(prop, ens_avg, file_normalizer)
        else:
            ens_avg = normalize_r(prop, ens_avg, file_normalizer)
        ens_avg[f"{prop}-sumrule_constant"] = scaler
        ens_avg["bin_center-norm"] = (
            ens_avg["bin_center"] / ens_avg["bin_center"].max()
        )

        # Add the state point's name and physical attributes
        ens_avg["ensemble_long"] = getattr(prop_info, "ensemble_long")
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(prop_info, attr_name)
        ens_avg[crowder_attr] = getattr(prop_info, crowder_attr)
        ens_avg["phi_bulk_c_round"] = ens_avg[crowder_attr].apply(
            lambda x: round_to_multiple(x, divisor, round_to)
        )

        # Bin centres in crowder units, recentred on the wall in a cylinder
        for label in crowders:
            crowder_size_attr = species_sizes[label]
            crowder_size = getattr(prop_info, crowder_size_attr)
            ens_avg[f"bin_center-{crowder_size_attr}"] = (
                2 * ens_avg["bin_center"] / crowder_size
            )
            if geometry == "cylindrical" and direction == "r":
                leftover = getattr(prop_info, "dcyl") % crowder_size
                ens_avg[f"bin_center-{crowder_size_attr}-recentered"] = (
                    ens_avg[f"bin_center-{crowder_size_attr}"]
                    - leftover / crowder_size
                )
                recentered = ens_avg["bin_center"] - leftover
                suffix = f"-{crowder_size_attr}" if len(crowders) > 1 else ""
                ens_avg[f"bin_center-recentered-norm{suffix}"] = (
                    recentered / recentered.max()
                )

        prop_csvs.append(ens_avg)

    prop_db = pd.concat(prop_csvs, axis=0)
    prop_db = prop_db.reset_index(drop=True)

    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = (
            save_to_space
            + f"{create_fullname(space, group, prop)}"
            + "-normalizedRescaled-space.csv"
        )
        prop_db.to_csv(filepath, index=False)
        logger.info("Saved %s.", filepath)

    return prop_db


def sum_rule_profiles(
    profiles: Mapping[EntityT, pd.DataFrame],
    prop: str,
    geometry: GeometryT,
    direction: DirectionT,
    *,
    by: Sequence[str] = ("ensemble_long",),
    normalizer: NormalizerT | None = None,
) -> pd.DataFrame:
    """
    Sum the scaled density profiles of several species, state point by
    state point, and normalize the sum by its own reference value: the
    left side of the density sum rule.

    With the species' scaled profiles :math:`\\phi_i(x) / a_i` from
    :func:`space_profiles`, the sum :math:`\\sum_i \\phi_i(x) / a_i` is
    constant in :math:`x` where the sum rule holds (Jeon et al., 2016,
    eqn (5); eqn (7) for several crowders). Normalized by its own value at
    the reference, the sum is one there and its departure from one
    elsewhere measures the rule's failure.

    The species' profiles must descend from the raw ensemble averages, as
    :func:`space_profiles` requires: the sum weights each species by its
    physical scale, and a species normalized on its own, by its sum, its
    integral, or its maximum, enters with an arbitrary weight. Nothing here
    can detect that, so it is the caller's contract.

    Parameters
    ----------
    profiles : mapping of (EntityT -> pd.DataFrame)
        Each species label mapped to its :func:`space_profiles` result for
        the same `prop`; each has the columns of `by`, ``bin_center``, and
        ``<prop><species>-scale``, in the same row order, with ``-scale``
        built from the raw ensemble-averaged profile.
    prop : str
        The direction and quantity, such as ``'rPhi'`` or ``'zRho'``.
    geometry : {'cubic', 'cylindrical'}
        The geometry of the study, which with `direction` chooses the
        default normalizer.
    direction : {'r', 'z'}
        The coordinate of the profiles.
    by : sequence of str, default ('ensemble_long',)
        The columns that identify a state point; the sum is normalized
        within each. An empty sequence treats the whole input as one.
    normalizer : :py:data:`~polyphys.manage.types.NormalizerT` or None, \
default None
        The reference value of each state point's sum. None takes the
        default of :data:`SUM_NORMALIZERS`: the value at :math:`r = 0` for
        ``'r'`` in a cylinder, the far-field (last-bin) value otherwise,
        where the monomer term vanishes and the sum is the crowder term
        :math:`\\phi_c(\\infty) / a_c`. A function is called with the state
        point's frame of the columns below. An explicit value is logged at
        INFO with the default it replaces.

    Returns
    -------
    pd.DataFrame
        One row per bin and state point, in the input order: the `by`
        columns, ``bin_center``, ``<prop>Sum-scale`` (the sum over species),
        ``<prop>Sum-normalizer`` (the reference value of the state point),
        and ``<prop>Sum-norm``. When every species carries
        ``<prop><species>-norm-sem``, also ``<prop>Sum-norm-sem``, the
        species' errors added in quadrature,
        :math:`\\sqrt{\\sum_i s_i^2}`, under the assumption that the
        species' errors are independent. It is the error of the sum of the
        normalized species, :math:`\\sum_i \\phi_i / (a_i N_i)` for each
        species' own normalizer :math:`N_i`; it equals the error of
        ``<prop>Sum-norm``, :math:`\\sum_i \\phi_i / a_i` over the sum's
        normalizer, only when every species shares one normalizer value.

    Raises
    ------
    ValueError
        If `profiles` is empty; if the species disagree on the `by` columns
        or on ``bin_center`` in any row; if (`geometry`, `direction`) has no
        entry in :data:`SUM_NORMALIZERS`; or if `normalizer` is none of the
        accepted forms.

    Warns
    -----
    UserWarning
        If a state point's reference value is zero; its ``<prop>Sum-norm``
        is then zero.

    Examples
    --------
    Two species over three radial bins of one state point in a cube: the
    monomer term falls to zero where the crowder term reaches its far-field
    value, and the sum is one where the rule holds:

    >>> keys = {'ensemble_long': ['s1'] * 3, 'bin_center': [0.5, 1.5, 2.5]}
    >>> monomers = pd.DataFrame({**keys, 'rPhiMon-scale': [0.8, 0.4, 0.0]})
    >>> crowders = pd.DataFrame({**keys, 'rPhiCrd-scale': [0.2, 0.6, 1.0]})
    >>> total = sum_rule_profiles({'Mon': monomers, 'Crd': crowders},
    ...                           'rPhi', 'cubic', 'r')
    >>> total[['bin_center', 'rPhiSum-scale', 'rPhiSum-norm']]
       bin_center  rPhiSum-scale  rPhiSum-norm
    0         0.5            1.0           1.0
    1         1.5            1.0           1.0
    2         2.5            1.0           1.0
    >>> total['rPhiSum-normalizer'].tolist()
    [1.0, 1.0, 1.0]

    Notes
    -----
    Time and space are :math:`O(n_\\text{species} \\, n_\\text{rows})` for
    the sum, the checks, and the error, plus :math:`O(n_\\text{rows})` per
    state point for the normalization. See :func:`space_profiles` for the
    references. The species' errors are correlated within a whole by
    excluded volume, which the quadrature sum ignores; the per-bin
    estimator is an open question.

    .. versionadded:: 0.6.0
    .. versionchanged:: 0.7.0
       Added ``<prop>Sum-norm-sem``.
    """
    if not profiles:
        raise ValueError("'profiles' is empty.")
    try:
        default = SUM_NORMALIZERS[(geometry, direction)]
    except KeyError:
        raise ValueError(
            f"No density profile along {direction!r} in a {geometry} "
            f"geometry; the (geometry, direction) pairs are "
            f"{sorted(SUM_NORMALIZERS)}."
        ) from None
    if normalizer is None:
        normalizer = default
    else:
        logger.info(
            "Normalizing the %s sum with %r instead of the %s %s default "
            "%r.", prop, normalizer, geometry, direction, default,
        )
    keys = [*by, "bin_center"]
    sum_prop = f"{prop}Sum"

    # The key columns of the first species fix the row order; every other
    # species must match them row for row.
    labels = list(profiles)
    first = profiles[labels[0]]
    total = first[keys].copy()
    total[f"{sum_prop}-scale"] = 0.0
    # The species' normalized errors add in quadrature when every species
    # has one; the column is left out otherwise.
    with_error = all(
        f"{prop}{label}-norm-sem" in profiles[label].columns
        for label in labels
    )
    if with_error:
        total[f"{sum_prop}-norm-sem"] = 0.0
    for label in labels:
        frame = profiles[label]
        if len(frame) != len(first) or not (
            frame[list(by)].to_numpy().tolist()
            == first[list(by)].to_numpy().tolist()
            and np.allclose(frame["bin_center"], first["bin_center"])
        ):
            raise ValueError(
                f"The {label!r} profiles do not match the {labels[0]!r} "
                f"profiles on {keys}."
            )
        total[f"{sum_prop}-scale"] += frame[f"{prop}{label}-scale"].to_numpy()
        if with_error:
            total[f"{sum_prop}-norm-sem"] += (
                frame[f"{prop}{label}-norm-sem"].to_numpy() ** 2
            )
    if with_error:
        total[f"{sum_prop}-norm-sem"] = np.sqrt(total[f"{sum_prop}-norm-sem"])

    # Normalize each state point's sum by its own reference value
    if by:
        groups = [group for _, group in total.groupby(list(by), sort=False)]
    else:
        groups = [total]
    normalized = [_normalize(sum_prop, group, normalizer) for group in groups]
    return pd.concat(normalized, axis=0).reset_index(drop=True)
