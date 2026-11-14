"""\
==========================================================
:mod:`polyphys.manage.parser`
==========================================================

This :mod:`~polyphys.manage.parser` module turns the name of a simulation
artifact into typed attributes. A name such as
``N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000ens1`` is a
run of keyword-value pairs; each project has its own keywords, and each level
of the lineage ``segment``, ``whole``, ``ensemble_long``, ``ensemble``, and
``space`` uses its own subset of them, in its own order.

A parser class is a declaration. It states the project's geometry, polymer
topology, and particle groups; one :class:`KeywordSpec` per filename keyword;
the keyword pattern of each lineage; fixed values; and the particle species
whose bulk number density and volume fraction the parser computes.
:class:`ParserBase` does the parsing, builds the parent lineage names, and
computes the bulk values.

Contents
========

- :class:`KeywordSpec` and :class:`Species`, the two building blocks of a
  declaration.
- :class:`ParserBase`, the base class every parser derives from.

The module ships no study's parser: a study declares its own, as the
"Studies" page shows.

.. versionchanged:: 1.0.0
   The parsers of nine studies left this module. ``TransFociCyl``,
   ``TransFociCub``, ``HnsCub``, and ``HnsCyl`` went to the ``projects/``
   folder of the PolyPhys repository, which is not installed; those of
   five unpublished studies were removed. Importing any of the four from
   here raises :exc:`ImportError` that says where it went.

.. versionchanged:: 1.1.0
   The four parsers are now in the companion repository
   `polyphys-lab <https://github.com/amirhs1/polyphys-lab>`_, as
   ``polyphys_lab.<study>.parser``, and the :exc:`ImportError` points
   there.

Dependencies
============
- Python standard library: `os`, `re`, `collections`, `dataclasses`,
  `decimal`, `pathlib`, `typing`.
- :mod:`polyphys.manage.utils`: :func:`~polyphys.manage.utils.invalid_keyword`
  and the number-density and volume-fraction functions.
- :mod:`polyphys.manage.types`: domain type aliases such as
  :py:obj:`~.LineageT`.

Examples
========
>>> from polyphys.tests.synthetic import ChainCyl
>>> artifact = ChainCyl('N200D10.0ac2.0nc0', 'ensemble', 'all')
>>> print(artifact.dcrowd)
2.0

>>> import polyphys.manage.parser as parser
>>> parser.HnsCub  # doctest: +ELLIPSIS
Traceback (most recent call last):
...
ImportError: HnsCub moved out of PolyPhys in 1.0.0, to polyphys_lab.hns_cub...

Notes
=====
The parsed and computed attributes are set on the instance at construction
with :func:`setattr`, so their names come from the declaration rather than
from the class body; :attr:`ParserBase.attributes` lists the ones a given
lineage sets.
"""
import os
import re
from collections import OrderedDict
from dataclasses import dataclass
from decimal import Context, Decimal
from pathlib import Path
from typing import ClassVar, Optional, Self
from .utils import (
    invalid_keyword,
    number_density_cube,
    number_density_cylinder,
    volume_fraction_cube,
    volume_fraction_cylinder
)
from .types import LineageT, TopologyT, GroupT, GeometryT, Pathish

__all__ = [
    'KeywordSpec',
    'ParserBase',
    'Species',
]

# The parsers of the published studies that left this module in 1.0.0,
# each with its package in polyphys-lab, polyphys_lab.<folder>.
_LAB_URL = 'https://github.com/amirhs1/polyphys-lab'
_MOVED_STUDIES: dict[str, str] = {
    'HnsCub': 'hns_cub',
    'HnsCyl': 'hns_cyl',
    'TransFociCub': 'trans_foci_cub',
    'TransFociCyl': 'trans_foci_cyl',
}

# A keyword is a run of letters (or hyphens); its value is the text up to the
# next keyword.
_TOKEN_RE = re.compile(r"([a-zA-Z\-]+)")

# A file extension, as ``ispath=True`` strips it: a dot and a run of letters
# and digits with at least one letter, such as ``.csv`` or ``.h5``, so that the
# fractional part of a trailing value, such as ``.5`` in ``ac1.5``, stays in
# the name.
_EXTENSION_RE = re.compile(r"\.[a-zA-Z0-9]*[a-zA-Z][a-zA-Z0-9]*")

# The exact decimal arithmetic of a keyword conversion; 50 digits hold the
# product of any filename value with a scale and an offset without rounding.
_DECIMAL_CONTEXT = Context(prec=50)

# The lineages whose names carry the box size, and so the bulk values.
_BULK_LINEAGES: tuple[LineageT, ...] = ('segment', 'whole', 'ensemble_long')

# The attributes that size the simulation box in each geometry.
_BOX_ATTRIBUTES: dict[GeometryT, tuple[str, ...]] = {
    'cubic': ('lcube',),
    'cylindrical': ('lcyl', 'dcyl'),
}


def _format_value(value: int | float) -> str:
    """
    Return the text of a keyword value in a short lineage name.

    A float that is a whole number prints as an integer; any other float
    prints as :func:`repr`, Python's shortest form that reads back as the
    same float, so two different values never print alike. An `int` prints
    unchanged. The cost is :math:`O(1)` per value.

    Parameters
    ----------
    value : int or float
        A parsed keyword value.

    Returns
    -------
    str
        The value as it appears after its keyword in an ``ensemble`` or
        ``space`` name.

    Examples
    --------
    >>> _format_value(2.0), _format_value(0.0), _format_value(40.0)
    ('2', '0', '40')
    >>> _format_value(1.5), _format_value(1.000001), _format_value(200)
    ('1.5', '1.000001', '200')
    """
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return repr(value)


@dataclass(frozen=True)
class KeywordSpec:
    """
    How a parser reads the value of one filename keyword.

    The value :math:`x` that follows the keyword in a name is stored in the
    attribute `attr` as :math:`y = s x + c`, with scale :math:`s` and offset
    :math:`c`; the default :math:`s = 1`, :math:`c = 0` stores the value as
    written.

    Parameters
    ----------
    attr : str
        Name of the parser attribute that stores the value.
    kind : {int, float}
        Type of the stored value.
    scale : float, default 1.0
        Scale :math:`s` of the conversion. Only a `float` keyword may set it.
    offset : float, default 0.0
        Offset :math:`c` of the conversion. Only a `float` keyword may set
        it.

    Raises
    ------
    TypeError
        If `kind` is neither `int` nor `float`, or if an `int` keyword sets
        `scale` or `offset`.

    Examples
    --------
    The keyword ``N`` is a number of monomers:

    >>> KeywordSpec('nmon', int).parse('200')
    200

    The keyword ``r`` is the radius of the cylinder on which the centers of
    the wall beads of size 1 lie, so the cylinder size is :math:`D = 2r - 1`:

    >>> spec = KeywordSpec('dcyl', float, scale=2.0, offset=-1.0)
    >>> spec.parse('10.5')
    20.0
    """
    attr: str
    kind: type[int] | type[float]
    scale: float = 1.0
    offset: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in (int, float):
            raise TypeError(
                f"KeywordSpec {self.attr!r}: 'kind' must be int or float, "
                f"not {self.kind!r}.")
        if self.kind is int and not self.is_identity:
            raise TypeError(
                f"KeywordSpec {self.attr!r}: only a float keyword can have a "
                "'scale' or an 'offset'.")

    @property
    def is_identity(self) -> bool:
        """
        Whether the conversion stores the value as written.
        """
        return self.scale == 1.0 and self.offset == 0.0

    def parse(self, token: str) -> int | float:
        """
        Convert the text that follows a keyword into the stored value.

        Parameters
        ----------
        token : str
            The text between the keyword and the next keyword, such as
            ``'10.5'`` or ``'1.'``.

        Returns
        -------
        int or float
            The value as an `int` for an `int` keyword, else
            :math:`s x + c`. A conversion that is not the identity is done in
            exact decimal arithmetic on the text of `token` and of
            :func:`repr` of :math:`s` and :math:`c`, then rounded once to the
            nearest float, so ``'4.1'`` with :math:`s = 2`, :math:`c = -1`
            gives ``7.2`` rather than the binary result
            ``7.199999999999999``. The cost is :math:`O(1)` per value.

        Raises
        ------
        ValueError
            If `token` is not a number, or if the keyword is an `int` keyword
            and `token` is not a whole number, such as ``'1.5'``.

        Examples
        --------
        >>> KeywordSpec('ncrowd', int).parse('1000.')
        1000
        >>> KeywordSpec('lcube', float, scale=2.0).parse('25.0')
        50.0
        >>> KeywordSpec('dcyl', float, scale=2.0, offset=-1.0).parse('4.1')
        7.2
        """
        try:
            value = float(token)
        except ValueError:
            raise ValueError(f"{token!r} is not a number") from None
        if self.kind is int:
            if not value.is_integer():
                raise ValueError(f"{token!r} is not an integer")
            return int(value)
        if self.is_identity:
            return value
        context = _DECIMAL_CONTEXT
        exact = context.add(
            context.multiply(Decimal(token), Decimal(repr(self.scale))),
            Decimal(repr(self.offset)))
        return float(exact)


@dataclass(frozen=True)
class Species:
    """
    A particle species whose bulk number density and volume fraction a parser
    computes.

    At the ``segment``, ``whole``, and ``ensemble_long`` lineages, a parser
    sets ``rho_bulk_<suffix>``, the number of particles of the species per
    unit volume of the box, and ``phi_bulk_<suffix>``, the fraction of the
    box volume they fill.

    Parameters
    ----------
    suffix : str
        Suffix of the two attribute names, such as ``'m'`` for monomers or
        ``'c'`` for crowders.
    number : str
        Name of the attribute that holds the number of particles.
    size : str
        Name of the attribute that holds the particle diameter.

    Examples
    --------
    >>> crowders = Species('c', 'ncrowd', 'dcrowd')
    >>> crowders.suffix
    'c'
    """
    suffix: str
    number: str
    size: str


class ParserBase:
    """
    Base parser class for extracting information from filenames or file paths
    in a structured project.

    A subclass declares a project's filename grammar and physics through the
    class attributes listed under "Declaration"; this class parses the name,
    builds the names of the parent lineages, and computes the bulk number
    densities and volume fractions. A declaration is checked when the
    subclass is created; a class whose `_patterns` is ``None`` is abstract and
    cannot be instantiated. :meth:`from_name` constructs a parser from a name
    token instead of a path, with ``filepath`` and ``ext`` set to ``None``.

    Parameters
    ----------
    artifact : Pathish
        The artifact to be parsed: a name, or, with ``ispath=True``, a
        relative or absolute file path (string or Path).
    lineage : :py:obj:`~.LineageT`
        The lineage of the artifact, specifying the hierarchical level within
        the project. Must be one of: ``'segment'``, ``'whole'``,
        ``'ensemble_long'``, ``'ensemble'``, or ``'space'``.
    group : :py:obj:`~.GroupT`
        The particle group the artifact describes; one of the subclass's
        `groups`.
    ispath : bool, default False
        Whether `artifact` is a file path. A path is expanded (``~`` and
        environment variables) and resolved, and its stem is parsed.

    Attributes
    ----------
    filepath : pathlib.Path or None
        The resolved absolute directory that holds `artifact` if it is a
        file path, or ``None`` if it is a name.
    filename : str
        The filename extracted from `artifact`, or `artifact` itself.
    ext : str or None
        The file extension, including the leading dot, if `artifact` is a
        file path: :attr:`pathlib.PurePath.suffix`, which is ``''`` for a
        file with no extension. ``None`` if `artifact` is a name.
    name : str
        The lineage name: `filename` without its extension, cut before
        ``.<group>`` (``segment`` and ``whole``) and before the first
        hyphen.
    group : :py:obj:`~.GroupT`
        Particle group type in the project.
    dump_keyword : str
        The attribute that holds the dump period, in time steps, of the
        `group` trajectory: ``dump_keywords[group]``.
    lineage : :py:obj:`~.LineageT`
        The lineage of the artifact.
    project_name : str
        The name of the parser subclass.
    lineage_genealogy : list of :py:obj:`~.LineageT`
        The artifact's `lineage` and its parent lineages.
    lineage_attributes : list of str
        The attributes parsed from `name`, in pattern order.
    physical_attributes : list of str
        The fixed and computed attributes listed for the artifact's
        `lineage`.
    attributes : list of str
        `lineage_attributes` followed by `physical_attributes`.
    segment : str
        The artifact's name at the ``segment`` lineage if it is in
        `lineage_genealogy`, else ``'N/A'``; cut from `name`, with no
        topology suffix.
    whole : str
        As `segment`, at the ``whole`` lineage.
    ensemble_long : str
        As `segment`, at the ``ensemble_long`` lineage.
    ensemble : str
        The artifact's name at the ``ensemble`` lineage if it is in
        `lineage_genealogy`, else ``'N/A'``; built from the parsed values,
        so a whole-number float prints without ``.0``.
    space : str
        As `ensemble`, at the ``space`` lineage.

    Class Attributes
    ----------------
    lineages : list of :py:obj:`~.LineageT`
        List of valid lineage types.
    genealogy : dict of :py:obj:`~.LineageT` to list of :py:obj:`~.LineageT`
        Each lineage mapped to itself and its parent lineages.
    geometry : :py:obj:`~.GeometryT`
        Specifies the geometry of the system.
    topology : :py:obj:`~.TopologyT`
        Specifies how particles are connected (or not) in the system.
    groups : list of :py:obj:`~.GroupT`
        Allowed particle groups in the system.
    dump_keywords : dict of :py:obj:`~.GroupT` to str
        Each group mapped to the attribute that holds the dump period, in
        time steps, of its trajectory.
    genealogy_attributes : dict of :py:obj:`~.LineageT` to OrderedDict
        Each lineage mapped to an OrderedDict of attribute names and their
        filename keywords, in pattern order; built from `_keywords` and
        `_patterns`.
    project_attributes : dict of :py:obj:`~.LineageT` to list of str
        Each lineage mapped to the attributes that are fixed or computed
        rather than parsed.

    Declaration
    -----------
    _geometry : :py:obj:`~.GeometryT`
        ``'cubic'`` or ``'cylindrical'``. A cubic project with species
        declares the attribute ``lcube``; a cylindrical one ``lcyl`` and
        ``dcyl``.
    _topology : :py:obj:`~.TopologyT`
        The polymer topology.
    _groups : list of :py:obj:`~.GroupT`
        The particle groups.
    _dump_keywords : dict of :py:obj:`~.GroupT` to str
        As `dump_keywords` above, one entry per group. Each attribute is set
        at ``segment``, ``whole``, and ``ensemble_long``, by the pattern or
        by `_defaults`.
    _keywords : dict of str to :class:`KeywordSpec`
        Each filename keyword mapped to how its value is read.
    _patterns : dict of :py:obj:`~.LineageT` to tuple of str
        Each lineage mapped to its keywords in name order.
    _project_attributes : dict of :py:obj:`~.LineageT` to list of str
        As `project_attributes` above.
    _defaults : dict of str to float, optional
        Attributes set at every lineage before parsing, such as a fixed
        monomer size.
    _species : tuple of :class:`Species`, optional
        The species whose bulk values are computed.
    _totals : dict of str to tuple of str, optional
        A suffix mapped to species suffixes; ``rho_bulk_<suffix>`` and
        ``phi_bulk_<suffix>`` are the sums of theirs.
    _masses : dict of str to str, optional
        A mass attribute mapped to a size attribute; the mass is the size
        cubed, :math:`m = d^3`, for a uniform mass density.

    Raises
    ------
    TypeError
        On instantiation of a class with no `_patterns`, such as
        `ParserBase` itself; on creation of a subclass whose declaration is
        incomplete or inconsistent.
    ValueError
        If `artifact` is empty; if `lineage` or `group` is not valid; if
        `name` lacks a keyword of the lineage's pattern; or if a keyword's
        value is not a number, or not a whole number for an `int` keyword.

    Notes
    -----
    Parsing splits `name` once into keyword-value pairs and keeps the first
    value of each keyword, so a name may carry keywords the pattern does not
    use. For a name of length :math:`L` and a pattern of :math:`k` keywords,
    parsing costs :math:`O(L + k)` time and :math:`O(L)` memory.

    The bulk values are computed only at the ``segment``, ``whole``, and
    ``ensemble_long`` lineages, whose names carry the box size.

    A value that is a function of parsed values but fits none of the tables
    above, such as a complementary percentage, is set by overriding
    ``_derived_attributes``, which runs at every lineage after the bulk
    values.

    Examples
    --------
    A parser for a hypothetical study of a linear chain among crowders in a
    periodic cube, whose ``bug`` and ``all`` trajectories are dumped every
    1000 and 5000 time steps:

    >>> class ChainCub(ParserBase):
    ...     _geometry = 'cubic'
    ...     _topology = 'linear'
    ...     _groups = ['bug', 'all']
    ...     _dump_keywords = {'bug': 'bdump', 'all': 'adump'}
    ...     _keywords = {
    ...         'N': KeywordSpec('nmon', int),
    ...         'ac': KeywordSpec('dcrowd', float),
    ...         'nc': KeywordSpec('ncrowd', int),
    ...         'l': KeywordSpec('lcube', float, scale=2.0),
    ...         'ens': KeywordSpec('ensemble_id', int),
    ...         'j': KeywordSpec('segment_id', int),
    ...     }
    ...     _patterns = {
    ...         'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
    ...         'whole': ('N', 'ac', 'nc', 'l', 'ens'),
    ...         'ensemble_long': ('N', 'ac', 'nc', 'l'),
    ...         'ensemble': ('N', 'ac', 'nc'),
    ...         'space': ('N', 'ac'),
    ...     }
    ...     _defaults = {'dmon': 1.0, 'bdump': 1000, 'adump': 5000}
    ...     _species = (Species('m', 'nmon', 'dmon'),
    ...                 Species('c', 'ncrowd', 'dcrowd'))
    ...     _project_attributes = {
    ...         'segment': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
    ...                     'phi_bulk_c', 'rho_bulk_c'],
    ...         'whole': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
    ...                   'phi_bulk_c', 'rho_bulk_c'],
    ...         'ensemble_long': ['dmon', 'phi_bulk_m', 'rho_bulk_m',
    ...                           'phi_bulk_c', 'rho_bulk_c'],
    ...         'ensemble': ['dmon'],
    ...         'space': ['dmon'],
    ...     }
    >>> artifact = ChainCub('N100ac2.0nc50l10.0ens1.bug', 'whole', 'bug')
    >>> artifact.lcube
    20.0
    >>> artifact.ensemble_long, artifact.ensemble, artifact.space
    ('N100ac2.0nc50l10.0', 'N100ac2nc50', 'N100ac2')
    >>> artifact.filepath is None, artifact.ext is None
    (True, True)
    >>> artifact.dump_keyword, getattr(artifact, artifact.dump_keyword)
    ('bdump', 1000)
    """
    _lineages: ClassVar[list[LineageT]] = \
        ['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
    _genealogy: ClassVar[dict[LineageT, list[LineageT]]] = {
        'segment': ['segment', 'whole', 'ensemble_long', 'ensemble', 'space'],
        'whole': ['whole', 'ensemble_long', 'ensemble', 'space'],
        'ensemble_long': ['ensemble_long', 'ensemble', 'space'],
        'ensemble': ['ensemble', 'space'],
        'space': ['space'],
    }
    # Every concrete subclass sets these; __init_subclass__ checks them.
    _geometry: ClassVar[GeometryT]
    _topology: ClassVar[TopologyT]
    _groups: ClassVar[list[GroupT]]
    _dump_keywords: ClassVar[dict[GroupT, str]]
    _keywords: ClassVar[Optional[dict[str, KeywordSpec]]] = None
    _patterns: ClassVar[Optional[dict[LineageT, tuple[str, ...]]]] = None
    _defaults: ClassVar[dict[str, float]] = {}
    _species: ClassVar[tuple[Species, ...]] = ()
    _totals: ClassVar[dict[str, tuple[str, ...]]] = {}
    _masses: ClassVar[dict[str, str]] = {}
    _genealogy_attributes: ClassVar[dict[LineageT, OrderedDict[str, str]]]
    _project_attributes: ClassVar[dict[LineageT, list[str]]]

    def __init_subclass__(cls) -> None:
        super().__init_subclass__()
        if cls._patterns is None:
            return
        keywords, patterns = cls._check_declaration()
        cls._genealogy_attributes = {
            lineage: OrderedDict(
                (keywords[keyword].attr, keyword)
                for keyword in patterns[lineage])
            for lineage in cls._lineages
        }

    @classmethod
    def _check_declaration(
        cls
    ) -> tuple[dict[str, KeywordSpec], dict[LineageT, tuple[str, ...]]]:
        """
        Check a subclass's declaration and return its keyword and pattern
        tables.

        Returns
        -------
        tuple of (dict, dict)
            `_keywords` and `_patterns`.

        Raises
        ------
        TypeError
            If a required class attribute is unset; if a lineage has no
            pattern or no project attributes; if a group has no dump
            keyword, or a dump keyword names no group; if a pattern keyword
            has no :class:`KeywordSpec`; if an attribute appears twice in
            one pattern; if the bulk recipe or a dump keyword reads an
            attribute a ``segment``, ``whole``, or ``ensemble_long`` name
            does not set; or if a total names an undeclared species.
        """
        name = cls.__name__
        keywords = cls._keywords
        patterns = cls._patterns
        geometry: Optional[GeometryT] = getattr(cls, '_geometry', None)
        if (geometry is None or getattr(cls, '_topology', None) is None
                or getattr(cls, '_groups', None) is None
                or getattr(cls, '_dump_keywords', None) is None
                or keywords is None or patterns is None
                or getattr(cls, '_project_attributes', None) is None):
            raise TypeError(
                f"{name} must declare '_geometry', '_topology', '_groups', "
                "'_dump_keywords', '_keywords', '_patterns', and "
                "'_project_attributes'.")
        for table in ('_patterns', '_project_attributes'):
            if set(getattr(cls, table)) != set(cls._lineages):
                raise TypeError(
                    f"{name}.{table} must have one entry per lineage: "
                    f"{cls._lineages}.")
        if set(cls._dump_keywords) != set(cls._groups):
            raise TypeError(
                f"{name}._dump_keywords must have one entry per group: "
                f"{cls._groups}.")
        for lineage, pattern in patterns.items():
            for keyword in pattern:
                if keyword not in keywords:
                    raise TypeError(
                        f"{name}: keyword {keyword!r} in the {lineage!r} "
                        "pattern has no KeywordSpec.")
            attrs = [keywords[keyword].attr for keyword in pattern]
            for attr in attrs:
                if attrs.count(attr) > 1:
                    raise TypeError(
                        f"{name}: attribute {attr!r} appears more than once "
                        f"in the {lineage!r} pattern.")
        suffixes = {species.suffix for species in cls._species}
        for total, parts in cls._totals.items():
            for part in parts:
                if part not in suffixes:
                    raise TypeError(
                        f"{name}: total {total!r} sums {part!r}, which is "
                        "not a declared species.")
        needed: list[str] = list(cls._masses.values())
        if cls._species:
            needed.extend(_BOX_ATTRIBUTES[geometry])
        for species in cls._species:
            needed.extend((species.number, species.size))
        for lineage in _BULK_LINEAGES:
            available = set(cls._defaults).union(
                keywords[keyword].attr for keyword in patterns[lineage])
            for attr in needed:
                if attr not in available:
                    raise TypeError(
                        f"{name}: the bulk values need {attr!r}, which the "
                        f"{lineage!r} pattern does not set.")
            for group, attr in cls._dump_keywords.items():
                if attr not in available:
                    raise TypeError(
                        f"{name}: the dump keyword {attr!r} of group "
                        f"{group!r} is not set by the {lineage!r} pattern "
                        "or '_defaults'.")
        return keywords, patterns

    def __init__(
        self,
        artifact: Pathish,
        lineage: LineageT,
        group: GroupT,
        ispath: bool = False
    ) -> None:
        if self._patterns is None:
            raise TypeError(
                f"{type(self).__name__} declares no lineage patterns; "
                "instantiate a subclass that declares a study.")
        artifact_str = os.fspath(artifact)
        if artifact_str == '':
            raise ValueError("'artifact' cannot be an empty string.")

        invalid_keyword(lineage, self.lineages)
        invalid_keyword(group, self.groups)

        self._ispath = ispath
        if ispath:
            # Normalize: expand env/~ and convert to absolute path
            p = Path(
                os.path.expandvars(os.path.expanduser(artifact_str))
            ).resolve()
            self._filepath: Optional[Path] = p.parent
            self._filename: str = p.name
            ext = p.suffix if _EXTENSION_RE.fullmatch(p.suffix) else ''
            self._name: str = p.name[:len(p.name) - len(ext)]
            self._ext: Optional[str] = ext
        else:
            self._filepath = None
            self._filename = artifact_str
            self._name = artifact_str
            self._ext = None

        self._project_name = self.__class__.__name__
        self._lineage: LineageT = lineage
        self._group: GroupT = group

        # Update the self._name by dropping 'group' and/or any other extra
        # beyond the template defined above.
        self._find_name()
        self._lineage_genealogy: list[LineageT] = self._genealogy[lineage]
        self._lineage_attributes = \
            list(self.genealogy_attributes[lineage].keys())
        self._physical_attributes = self.project_attributes[lineage]
        self._attributes = \
            self._lineage_attributes + self._physical_attributes

        self._set_defaults()
        self._parse_name()
        self._set_parents()
        if self.lineage in _BULK_LINEAGES:
            self._dependent_attributes()
        self._derived_attributes()

    def __str__(self) -> str:
        """
        Provide a formatted summary of the parser instance.
        """
        observation = (
            f"Artifact:\n"
            f"    Name: '{self.filename}',\n"
            f"    Geometry: '{self.geometry}',\n"
            f"    Group: '{self.group}',\n"
            f"    Lineage: '{self.lineage}',\n"
            f"    Topology: '{self.topology}',\n"
            f"    Project: '{self.project_name}'"
        )
        return observation

    def __repr__(self) -> str:
        return (
            f"Artifact(filename='{self.filename}', "
            f"geometry='{self.geometry}', "
            f"group='{self.group}', lineage='{self.lineage}', "
            f"topology='{self.topology}', project='{self.project_name}')"
        )

    @classmethod
    def from_name(
        cls,
        name: str,
        lineage: LineageT,
        group: GroupT
    ) -> Self:
        """
        Construct a parser from a *name token* (not a path).

        Sets `filepath` and `ext` to ``None`` and keeps the token as both
        `filename` and `name` before `_find_name()` normalization.

        Examples
        --------
        >>> from polyphys.tests.synthetic import ChainCyl
        >>> parser = ChainCyl.from_name('N200D10.0ac2.0', 'space', 'bug')
        >>> parser.filepath is None, parser.ext is None
        (True, True)
        """
        return cls(name, lineage, group, ispath=False)

    @classmethod
    def pattern(cls, lineage: LineageT) -> tuple[str, ...]:
        """
        Return the filename keywords of `lineage`, in name order.

        Parameters
        ----------
        lineage : :py:obj:`~.LineageT`
            One of ``'segment'``, ``'whole'``, ``'ensemble_long'``,
            ``'ensemble'``, or ``'space'``.

        Returns
        -------
        tuple of str
            The keywords, such as ``'N'`` or ``'ac'``, that a name of
            `lineage` carries, in the order they appear in it.

        Raises
        ------
        TypeError
            If the class declares no patterns, such as `ParserBase` itself.
        ValueError
            If `lineage` is not a valid lineage.

        Examples
        --------
        >>> from polyphys.tests.synthetic import ChainCyl
        >>> ChainCyl.pattern('space')
        ('N', 'D', 'ac')
        """
        patterns = cls._declared_patterns()
        invalid_keyword(lineage, cls._lineages)
        return patterns[lineage]

    @classmethod
    def attributes_at(cls, lineage: LineageT) -> frozenset[str]:
        """
        Return the names of the attributes that a parser of `lineage` sets.

        The set holds the attributes parsed from the lineage's pattern, the
        `_defaults`, the five lineage names (``segment`` to ``space``), and
        the lineage's `_project_attributes`. At ``segment``, ``whole``, and
        ``ensemble_long`` it also holds the bulk number density
        ``rho_bulk_<suffix>`` and volume fraction ``phi_bulk_<suffix>`` of
        each species and each total, and the declared masses. An attribute
        that an override of ``_derived_attributes`` sets counts only if
        `_project_attributes` lists it.

        The set is built from the declaration, without parsing a name; the
        cost is :math:`O(k + n)` for a pattern of :math:`k` keywords and
        :math:`n` declared defaults, species, totals, masses, and project
        attributes.

        Parameters
        ----------
        lineage : :py:obj:`~.LineageT`
            One of ``'segment'``, ``'whole'``, ``'ensemble_long'``,
            ``'ensemble'``, or ``'space'``.

        Returns
        -------
        frozenset of str
            The attribute names.

        Raises
        ------
        TypeError
            If the class declares no patterns, such as `ParserBase` itself.
        ValueError
            If `lineage` is not a valid lineage.

        Examples
        --------
        The bulk values exist at ``ensemble_long``, not at ``space``:

        >>> from polyphys.tests.synthetic import ChainCyl
        >>> 'phi_bulk_c' in ChainCyl.attributes_at('ensemble_long')
        True
        >>> 'phi_bulk_c' in ChainCyl.attributes_at('space')
        False
        """
        cls._declared_patterns()
        invalid_keyword(lineage, cls._lineages)
        names = set(cls._genealogy_attributes[lineage])
        names.update(cls._defaults, cls._lineages,
                     cls._project_attributes[lineage])
        if lineage in _BULK_LINEAGES:
            suffixes = [species.suffix for species in cls._species]
            suffixes.extend(cls._totals)
            for suffix in suffixes:
                names.update((f"rho_bulk_{suffix}", f"phi_bulk_{suffix}"))
            names.update(cls._masses)
        return frozenset(names)

    @classmethod
    def _declared_patterns(cls) -> dict[LineageT, tuple[str, ...]]:
        """
        Return `_patterns`, or raise `TypeError` if the class declares none.
        """
        if cls._patterns is None:
            raise TypeError(
                f"{cls.__name__} declares no lineage patterns; use a "
                "subclass that declares a study.")
        return cls._patterns

    @property
    def lineages(self) -> list[LineageT]:
        """
        List of all the acceptable lineages.
        """
        return self._lineages

    @property
    def genealogy(self) -> dict[LineageT, list[LineageT]]:
        """
        Dictionary of lineages and their parent lineages.
        """
        return self._genealogy

    @property
    def geometry(self) -> GeometryT:
        """
        System geometry in a molecular dynamics system.
        """
        return self._geometry

    @property
    def groups(self) -> list[GroupT]:
        """
        List of valid group names for the subclass.
        """
        return self._groups

    @property
    def dump_keywords(self) -> dict[GroupT, str]:
        """
        Each group mapped to the attribute that holds the dump period, in
        time steps, of its trajectory.
        """
        return self._dump_keywords

    @property
    def topology(self) -> TopologyT:
        """
        Define the polymer topology for the parser subclass.
        """
        return self._topology

    @property
    def genealogy_attributes(self) -> dict[LineageT, OrderedDict[str, str]]:
        """
        Dictionary of lineage-specific attributes. Each key is a lineage type,
        and each value is an OrderedDict mapping attribute names to their
        filename keywords, in pattern order.
        """
        return self._genealogy_attributes

    @property
    def project_attributes(self) -> dict[LineageT, list[str]]:
        """
        Dictionary of project attributes. Each key is a lineage type, and
        each value is the list of attributes that are fixed or computed
        rather than parsed at that lineage.
        """
        return self._project_attributes

    @property
    def filename(self) -> str:
        """
        Return the filename, either extracted from the path or the name
        itself.
        """
        return self._filename

    @property
    def filepath(self) -> Optional[Path]:
        """
        Return the resolved absolute directory that holds the artifact, or
        ``None`` if the parser was built from a name.

        Examples
        --------
        >>> from polyphys.tests.synthetic import ChainCyl
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as folder:
        ...     path = Path(folder, 'N200D10.0ac2.0.csv')
        ...     path.touch()
        ...     parser = ChainCyl(path, 'space', 'bug', ispath=True)
        ...     parser.filepath == Path(folder).resolve()
        True
        >>> ChainCyl('N200D10.0ac2.0', 'space', 'bug').filepath is None
        True
        """
        return self._filepath

    @property
    def ispath(self) -> bool:
        """
        Return whether the artifact was parsed as a filesystem path.
        """
        return self._ispath

    @property
    def group(self) -> GroupT:
        """
        Return the current group.
        """
        return self._group

    @property
    def dump_keyword(self) -> str:
        """
        The attribute that holds the dump period, in time steps, of the
        `group` trajectory.
        """
        return self._dump_keywords[self._group]

    @property
    def ext(self) -> Optional[str]:
        """
        Return the file extension, including the leading dot, or ``None`` if
        the parser was built from a name.

        The extension is :attr:`pathlib.PurePath.suffix`, the part from the
        last dot, when that part holds a letter, such as ``.csv`` or ``.h5``;
        otherwise it is ``''`` and the part stays in :attr:`name`. A name
        that ends with a value, such as ``ac1.5``, or with a dot, such as
        ``ac1.``, so keeps its last value whatever the Python version.

        Examples
        --------
        >>> from polyphys.tests.synthetic import ChainCyl
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as folder:
        ...     path = Path(folder, 'N200D10ac2')
        ...     path.touch()
        ...     ChainCyl(path, 'space', 'bug', ispath=True).ext
        ''
        >>> with tempfile.TemporaryDirectory() as folder:
        ...     path = Path(folder, 'N200D10ac1.5')
        ...     path.touch()
        ...     parser = ChainCyl(path, 'space', 'bug', ispath=True)
        ...     parser.ext, parser.dcrowd
        ('', 1.5)
        >>> ChainCyl('N200D10.0ac2.0', 'space', 'bug').ext is None
        True
        """
        return self._ext

    @property
    def lineage(self) -> LineageT:
        """
        Return the current lineage.
        """
        return self._lineage

    @property
    def name(self) -> str:
        """
        Return the unique name parsed from the filename.
        """
        return self._name

    @property
    def project_name(self) -> str:
        """
        Return the project (parser class) name,
        """
        return self._project_name

    @property
    def attributes(self) -> list[str]:
        """
        Return lineage-specific and project attributes for an artifact.
        """
        return self._attributes

    @property
    def lineage_genealogy(self) -> list[LineageT]:
        """
        Return the parents of a given `lineage`.
        """
        return self._lineage_genealogy

    @property
    def lineage_attributes(self) -> list[str]:
        """
        Return lineage-specific attributes for an artifact with a given
        `lineage`.
        """
        return self._lineage_attributes

    @property
    def physical_attributes(self) -> list[str]:
        """
        Return project-level attributes for an artifact with a given
        `lineage`.
        """
        return self._physical_attributes

    def _find_name(self) -> None:
        """
        parses the unique lineage_name (the first substring of filename
        and/or the segment keyword middle substring) of a filename.
        """
        if self.lineage in ['segment', 'whole']:
            # a 'segment' lineage only used in 'probe' phase
            # a 'whole' lineage used in 'probe' or 'analyze' phases
            # so its lineage_name is either ended by 'group' keyword or '-'.
            # these two combined below:
            base = self.name.split(f".{self.group}")[0]
            self._name = base.split('-')[0]
        else:  # 'ensemble' or 'space' lineages
            self._name = self.name.split('-')[0]

    def _set_defaults(self) -> None:
        """
        Set the attributes declared in `_defaults`.
        """
        for attr, value in self._defaults.items():
            setattr(self, attr, value)

    def _parse_name(self) -> None:
        """
        Set the attributes of the lineage's pattern from the values in
        `name`.

        Raises
        ------
        ValueError
            If `name` lacks a keyword of the pattern, or a keyword's value is
            not a number, or not a whole number for an `int` keyword.
        """
        keywords = self._keywords
        patterns = self._patterns
        assert keywords is not None and patterns is not None  # __init__
        words = _TOKEN_RE.split(self.name)
        values: dict[str, str] = {}
        for keyword, token in zip(words[1::2], words[2::2]):
            values.setdefault(keyword, token)
        for keyword in patterns[self.lineage]:
            spec = keywords[keyword]
            if keyword not in values:
                raise ValueError(
                    f"keyword {keyword!r} for attribute {spec.attr!r} not "
                    f"found in {self.name!r} (lineage {self.lineage!r}).")
            try:
                value = spec.parse(values[keyword])
            except ValueError as err:
                raise ValueError(
                    f"value of keyword {keyword!r} in {self.name!r} "
                    f"(lineage {self.lineage!r}): {err}.") from err
            setattr(self, spec.attr, value)

    def _set_parents(self) -> None:
        """
        Set parent lineage names for each lineage type, following the
        hierarchy defined in `_genealogy`.

        Notes
        -----
        The `self._genealogy` defines the following parent-child hierarchy:

            - space -> ensemble -> ensemble_long -> whole -> segment

        Each lineage on the left has all the lineages on its right.

        The long-form names (``segment``, ``whole``, and ``ensemble_long``)
        are cut from `name` after the value of the parent's last keyword, so
        they keep the original tokens. `_parse_name` converts some values,
        such as ``r`` to `dcyl` and ``l`` to `lcube`; writing those
        attributes back under the original keywords would convert them again
        at every re-parse. Because the cut ends at a keyword's value, a
        ``segment`` or ``whole`` name never carries the topology suffix:
        ``...ens1.ring`` and ``...ens1.j07.ring`` give the ``whole``
        ``...ens1``, so a simulation has one ``whole`` name whichever file it
        was parsed from.

        The short-form names (``ensemble`` and ``space``) use their own
        keywords, such as ``D``, and are built from the attributes. A float
        value that is a whole number prints as an integer and any other
        float as :func:`repr` (``ac2``, ``epshc0``, ``epshc1.5``), whatever
        the name being parsed wrote, such as ``ac2.0`` or ``ac2``.

        Every parent attribute is read once, so an unparsed keyword raises
        `AttributeError`. The cost is :math:`O(L)` in the length :math:`L`
        of `name`.
        """
        for lineage_name in self.lineages:
            lineage_value = 'N/A'
            if lineage_name in self.genealogy[self.lineage]:
                lineage_attr = self.genealogy_attributes[lineage_name]
                lineage_value = ''.join(
                    f"{attr_short}{_format_value(getattr(self, attr_long))}"
                    for attr_long, attr_short in lineage_attr.items()
                )
                if lineage_name in ['segment', 'whole', 'ensemble_long']:
                    last_keyword = next(reversed(lineage_attr.values()))
                    lineage_value = self._name_through(last_keyword)
            setattr(self, lineage_name, lineage_value)

    def _name_through(self, keyword: str) -> str:
        """
        Return `name` up to and including the value of `keyword`, without a
        trailing dot.

        Parameters
        ----------
        keyword : str
            A short-form keyword, such as ``ens``, that appears in `name`.

        Returns
        -------
        str
            The prefix of `name` that ends with the value of `keyword`.

        Raises
        ------
        ValueError
            If `keyword` is not a keyword of `name`.
        """
        for match in re.finditer(r"([a-zA-Z]+)([^a-zA-Z]*)", self.name):
            if match.group(1) == keyword:
                return self.name[:match.end()].rstrip('.')
        raise ValueError(f"'{keyword}' keyword not found in '{self.name}'.")

    def _bulk_values(self, number: float, size: float) -> tuple[float, float]:
        """
        Return the bulk number density and volume fraction of `number`
        particles of diameter `size` in the periodic box.

        Parameters
        ----------
        number : float
            Number of particles.
        size : float
            Particle diameter.

        Returns
        -------
        tuple of (float, float)
            The number density and the volume fraction, from
            :func:`~polyphys.manage.utils.number_density_cube` and
            :func:`~polyphys.manage.utils.volume_fraction_cube` with the side
            `lcube`, or their cylinder counterparts with `lcyl` and `dcyl`,
            all with ``pbc=True``.
        """
        if self.geometry == 'cubic':
            lcube = getattr(self, 'lcube')
            return (number_density_cube(number, size, lcube, pbc=True),
                    volume_fraction_cube(number, size, lcube, pbc=True))
        lcyl = getattr(self, 'lcyl')
        dcyl = getattr(self, 'dcyl')
        return (number_density_cylinder(number, size, lcyl, dcyl, pbc=True),
                volume_fraction_cylinder(number, size, lcyl, dcyl, pbc=True))

    def _dependent_attributes(self) -> None:
        """
        Compute the bulk values of each species, their declared totals, and
        the declared masses.

        Raises
        ------
        ValueError
            If a species leaves its centers no room in the box, such as a
            species as wide as the cylinder; the message names the species,
            `name`, and the lineage.
        """
        for species in self._species:
            try:
                rho, phi = self._bulk_values(getattr(self, species.number),
                                             getattr(self, species.size))
            except ValueError as err:
                raise ValueError(
                    f"species {species.suffix!r} of {self.name!r} "
                    f"(lineage {self.lineage!r}): {err}") from err
            setattr(self, f"rho_bulk_{species.suffix}", rho)
            setattr(self, f"phi_bulk_{species.suffix}", phi)
        for total, parts in self._totals.items():
            for quantity in ('rho_bulk', 'phi_bulk'):
                setattr(self, f"{quantity}_{total}",
                        sum(getattr(self, f"{quantity}_{part}")
                            for part in parts))
        for mass, size in self._masses.items():
            setattr(self, mass, getattr(self, size) ** 3)

    def _derived_attributes(self) -> None:
        """
        Set the attributes that a subclass derives from parsed values outside
        the declaration tables; the base class sets none.

        It runs at every lineage, after parsing, the parent names, and the
        bulk values, so an override reads only attributes that the lineage's
        pattern, `_defaults`, or the bulk recipe set. An override costs
        :math:`O(1)` per attribute it sets.
        """


def __getattr__(name: str) -> object:
    """
    Say where a study parser that left this module in 1.0.0 went.

    Python calls this (:pep:`562`) only for a name the module does not
    define. An :exc:`ImportError`, not an :exc:`AttributeError`, so that
    ``from polyphys.manage.parser import HnsCub`` shows the message: Python
    replaces an :exc:`AttributeError` raised here with its own.

    Raises
    ------
    ImportError
        If `name` is the parser of a published study that left the
        module.
    AttributeError
        For any other name.
    """
    if name in _MOVED_STUDIES:
        raise ImportError(
            f"{name} moved out of PolyPhys in 1.0.0, to "
            f"polyphys_lab.{_MOVED_STUDIES[name]}.parser in the "
            f"polyphys-lab repository, {_LAB_URL}; copy it into your "
            f"study's code.")
    raise AttributeError(
        f"module 'polyphys.manage.parser' has no attribute {name!r}")
