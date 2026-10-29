"""\
==========================================================
:mod:`polyphys.manage.parser`
==========================================================

This :mod:`~polyphys.manage.parser` module turns the name of a simulation
artifact into typed attributes. A name such as
``N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000ens1`` is a
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
- Cylindrical confinement: :class:`SumRuleCyl`, :class:`TransFociCyl`, and
  :class:`HnsCyl`.
- Cubic box: :class:`TwoMonDepCub`, :class:`TransFociCub`, :class:`HnsCub`,
  :class:`SumRuleCubHeteroRing`, :class:`SumRuleCubHeteroLinear`, and
  :class:`SumRuleCubHeteroRingBiDisperse`.

Dependencies
============
- Python standard library: `os`, `re`, `collections`, `dataclasses`,
  `pathlib`, `typing`.
- :mod:`polyphys.manage.utils`: :func:`~polyphys.manage.utils.invalid_keyword`
  and the number-density and volume-fraction functions.
- :mod:`polyphys.manage.types`: domain type aliases such as
  :py:data:`~polyphys.manage.types.LineageT`.

Examples
========
>>> artifact = SumRuleCyl('N200D10.0ac2.0nc0', 'ensemble', 'all')
>>> print(artifact.dcrowd)
2.0

>>> artifact = HnsCyl('N200D20.0nh16ac1.0epshc1.0nc0', 'space', 'nucleoid')
>>> print(artifact.eps_hc)
1.0

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

# A keyword is a run of letters (or hyphens); its value is the text up to the
# next keyword.
_TOKEN_RE = re.compile(r"([a-zA-Z\-]+)")

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

    The keyword ``r`` is the radius of the cylinder on which the centres of
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
            ``scale * float(token) + offset``.

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
        return self.scale * value + self.offset


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
    densities and volume fractions.

    Parameters
    ----------
    artifact : Pathish
        The artifact to be parsed: a name, or, with ``ispath=True``, a
        relative or absolute file path (string or Path).
    lineage : :py:data:`LineageT`
        The lineage of the artifact, specifying the hierarchical level within
        the project. Must be one of: ``'segment'``, ``'whole'``,
        ``'ensemble_long'``, ``'ensemble'``, or ``'space'``.
    group : :py:data:`GroupT`
        The particle group the artifact describes; one of the subclass's
        `groups`.
    ispath : bool, default False
        Whether `artifact` is a file path. A path is expanded (``~`` and
        environment variables) and resolved, and its stem is parsed.

    Alternative Constructors
    ------------------------
    .. method:: from_name(name, lineage, group)
       Construct a parser directly from a *name token* (not a path).
       This bypasses filesystem resolution. In this mode, ``filepath`` and
       ``ext`` are ``None``, and the given token is treated as both
       the ``filename`` and ``name`` prior to :meth:`_find_name`
       normalization.

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
    group : :py:data:`GroupT`
        Particle group type in the project.
    lineage : :py:data:`LineageT`
        The lineage of the artifact.
    project_name : str
        The name of the parser subclass.
    lineage_genealogy : list of :py:data:`LineageT`
        The artifact's `lineage` and its parent lineages.
    lineage_attributes : list of str
        The attributes parsed from `name`, in pattern order.
    physical_attributes : list of str
        The fixed and computed attributes listed for the artifact's
        `lineage`.
    attributes : list of str
        `lineage_attributes` followed by `physical_attributes`.
    segment, whole, ensemble_long, ensemble, space : str
        The artifact's name at each lineage of `lineage_genealogy`, and
        ``'N/A'`` at the others. ``segment``, ``whole``, and
        ``ensemble_long`` are cut from `name` and carry no topology suffix;
        ``ensemble`` and ``space`` print a whole-number float without
        ``.0``. :meth:`_set_parents` gives both rules.

    Class Attributes
    ----------------
    lineages : list of :py:data:`LineageT`
        List of valid lineage types.
    genealogy : dict of (:py:data:`LineageT` -> list of :py:data:`LineageT`)
        Each lineage mapped to itself and its parent lineages.
    geometry : :py:data:`GeometryT`
        Specifies the geometry of the system.
    topology : :py:data:`TopologyT`
        Specifies how particles are connected (or not) in the system.
    groups : list of :py:data:`GroupT`
        Allowed particle groups in the system.
    genealogy_attributes : dict of (\
:py:data:`LineageT` -> OrderedDict[str, str])
        Each lineage mapped to an OrderedDict of attribute names and their
        filename keywords, in pattern order; built from `_keywords` and
        `_patterns`.
    project_attributes : dict of (:py:data:`LineageT` -> list of str)
        Each lineage mapped to the attributes that are fixed or computed
        rather than parsed.

    Declaration
    -----------
    A subclass sets these class attributes:

    _geometry : :py:data:`GeometryT`
        ``'cubic'`` or ``'cylindrical'``. A cubic project with species
        declares the attribute ``lcube``; a cylindrical one ``lcyl`` and
        ``dcyl``.
    _topology : :py:data:`TopologyT`
        The polymer topology.
    _groups : list of :py:data:`GroupT`
        The particle groups.
    _keywords : dict of (str -> :class:`KeywordSpec`)
        Each filename keyword mapped to how its value is read.
    _patterns : dict of (:py:data:`LineageT` -> tuple of str)
        Each lineage mapped to its keywords in name order.
    _project_attributes : dict of (:py:data:`LineageT` -> list of str)
        As `project_attributes` above.
    _defaults : dict of (str -> float), optional
        Attributes set at every lineage before parsing, such as a fixed
        monomer size.
    _species : tuple of :class:`Species`, optional
        The species whose bulk values are computed.
    _totals : dict of (str -> tuple of str), optional
        A suffix mapped to species suffixes; ``rho_bulk_<suffix>`` and
        ``phi_bulk_<suffix>`` are the sums of theirs.
    _masses : dict of (str -> str), optional
        A mass attribute mapped to a size attribute; the mass is the size
        cubed, :math:`m = d^3`, for a uniform mass density.

    A declaration is checked when the subclass is created; a class whose
    `_patterns` is ``None`` is abstract and cannot be instantiated.

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
    :meth:`_derived_attributes`, which runs at every lineage after the bulk
    values.

    Examples
    --------
    A parser for a hypothetical study of a linear chain among crowders in a
    periodic cube:

    >>> class ChainCub(ParserBase):
    ...     _geometry = 'cubic'
    ...     _topology = 'linear'
    ...     _groups = ['bug', 'all']
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
    ...     _defaults = {'dmon': 1.0}
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
            pattern or no project attributes; if a pattern keyword has no
            :class:`KeywordSpec`; if an attribute appears twice in one
            pattern; if the bulk recipe reads an attribute a ``segment``,
            ``whole``, or ``ensemble_long`` name does not set; or if a total
            names an undeclared species.
        """
        name = cls.__name__
        keywords = cls._keywords
        patterns = cls._patterns
        geometry: Optional[GeometryT] = getattr(cls, '_geometry', None)
        if (geometry is None or getattr(cls, '_topology', None) is None
                or getattr(cls, '_groups', None) is None
                or keywords is None or patterns is None
                or getattr(cls, '_project_attributes', None) is None):
            raise TypeError(
                f"{name} must declare '_geometry', '_topology', '_groups', "
                "'_keywords', '_patterns', and '_project_attributes'.")
        for table in ('_patterns', '_project_attributes'):
            if set(getattr(cls, table)) != set(cls._lineages):
                raise TypeError(
                    f"{name}.{table} must have one entry per lineage: "
                    f"{cls._lineages}.")
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
                "instantiate a concrete parser such as SumRuleCyl.")
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
            self._name: str = p.stem
            self._ext: Optional[str] = p.suffix
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
        >>> parser = SumRuleCyl.from_name('N200D10.0ac2.0', 'space', 'bug')
        >>> parser.filepath is None, parser.ext is None
        (True, True)
        """
        return cls(name, lineage, group, ispath=False)

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
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as folder:
        ...     path = Path(folder, 'N200D10.0ac2.0.csv')
        ...     path.touch()
        ...     parser = SumRuleCyl(path, 'space', 'bug', ispath=True)
        ...     parser.filepath == Path(folder).resolve()
        True
        >>> SumRuleCyl('N200D10.0ac2.0', 'space', 'bug').filepath is None
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
    def ext(self) -> Optional[str]:
        """
        Return the file extension, including the leading dot, or ``None`` if
        the parser was built from a name.

        The extension is :attr:`pathlib.PurePath.suffix`: only the part after
        the last dot, and ``''`` for a file with no extension.

        Examples
        --------
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as folder:
        ...     path = Path(folder, 'N200D10ac2')
        ...     path.touch()
        ...     SumRuleCyl(path, 'space', 'bug', ispath=True).ext
        ''
        >>> SumRuleCyl('N200D10.0ac2.0', 'space', 'bug').ext is None
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
        """
        for species in self._species:
            rho, phi = self._bulk_values(getattr(self, species.number),
                                         getattr(self, species.size))
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


class TwoMonDepCub(ParserBase):
    """
    Parse an artifact name of the *TwoMonDep* project: two monomers held at
    a fixed surface-to-surface distance among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: am#nm#ac#nc#hl#sd#dt#bdump#adump#tdump#ens#.j#
    - ``whole``: am#nm#ac#nc#hl#sd#dt#bdump#adump#tdump#ens#
    - ``ensemble_long``: am#nm#ac#nc#hl#sd#dt#bdump#adump#tdump#
    - ``ensemble``: nm#am#ac#nc#sd#
    - ``space``: nm#am#ac#nc#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the two monomers, ``'all'`` for every
        particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon : float
        Monomer diameter; keyword ``am``.
    nmon : int
        Number of monomers; keyword ``nm``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    lcube : float
        Side of the cube, :math:`L = 2 l_h`, from the half-side :math:`l_h`
        (keyword ``hl``).
    d_sur : float
        Surface-to-surface distance between the two monomers; keyword
        ``sd``.
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    tdump : int
        Output period of the LAMMPS thermo values, in time steps; keyword
        ``tdump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of monomers.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    Examples
    --------
    >>> artifact = TwoMonDepCub('nm2am5.0ac1.0nc1000', 'space', 'bug')
    >>> print(artifact.nmon)
    2
    """
    _geometry = 'cubic'
    _topology = 'atomic'
    _groups = ['bug', 'all']
    _keywords = {
        'am': KeywordSpec('dmon', float),
        'nm': KeywordSpec('nmon', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'hl': KeywordSpec('lcube', float, scale=2.0),
        'sd': KeywordSpec('d_sur', float),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'tdump': KeywordSpec('tdump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('am', 'nm', 'ac', 'nc', 'hl', 'sd', 'dt', 'bdump',
                    'adump', 'tdump', 'ens', 'j'),
        'whole': ('am', 'nm', 'ac', 'nc', 'hl', 'sd', 'dt', 'bdump', 'adump',
                  'tdump', 'ens'),
        'ensemble_long': ('am', 'nm', 'ac', 'nc', 'hl', 'sd', 'dt', 'bdump',
                          'adump', 'tdump'),
        'ensemble': ('nm', 'am', 'ac', 'nc', 'sd'),
        'space': ('nm', 'am', 'ac', 'nc'),
    }
    _species = (Species('m', 'nmon', 'dmon'), Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                          'rho_bulk_c'],
        'ensemble': [],
        'space': []
    }


class SumRuleCyl(ParserBase):
    """
    Parse an artifact name of the *SumRuleCyl* project: a linear chain among
    crowders in a cylinder that is periodic along its axis.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
    - ``whole``: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
    - ``ensemble_long``: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
    - ``ensemble``: N#D#ac#nc#
    - ``space``: N#D#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the chain, ``'all'`` for every
        particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon : int
        Monomer diameter, fixed at 1; not in the name.
    nmon : int
        Number of monomers; keyword ``N``.
    dcrowd : float
        Crowder diameter; keyword ``sig`` in the long names, ``ac`` in the
        short ones.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    lcyl : float
        Length of the cylinder, :math:`Z = 2 l_z`, from the half-length
        :math:`l_z` (keyword ``lz``).
    dcyl : float
        Size of the cylinder, :math:`D = 2r - 1`, from the radius :math:`r`
        (keyword ``r``) in the long names; :math:`D` itself (keyword ``D``)
        in the short ones.
    epsilon : float
        Strength of the wall-particle LJ interaction; keyword ``epsilon``.
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of monomers.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    The cylindrical wall is a layer of beads of size :math:`a_w = 1` whose
    centres lie on the cylinder of radius :math:`r` set in LAMMPS (keyword
    ``r``), so the cylinder size is :math:`D = 2r - a_w`, stored as `dcyl`.
    The keyword ``lz`` is the half-length of the cylinder, so its length is
    :math:`Z = 2 l_z`, stored as `lcyl`. The centres of particles of
    species :math:`\\alpha` with size :math:`a_\\alpha` are confined to a
    cylinder of size :math:`D_{\\mathrm{eff},\\alpha} = D - a_\\alpha`, and
    the bulk number density and volume fraction of that species both use
    the effective volume :math:`V_{\\mathrm{eff},\\alpha} = \\pi Z
    D_{\\mathrm{eff},\\alpha}^2 / 4`.

    Examples
    --------
    >>> artifact = SumRuleCyl('N200D10.0ac2.0nc0', 'ensemble', 'all')
    >>> print(artifact.dcrowd)
    2.0
    """
    _geometry = 'cylindrical'
    _topology = 'linear'
    _groups = ['bug', 'all']
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'epsilon': KeywordSpec('epsilon', float),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'sig': KeywordSpec('dcrowd', float),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'epsilon', 'r', 'lz', 'sig', 'nc', 'dt', 'bdump',
                    'adump', 'ens', 'j'),
        'whole': ('N', 'epsilon', 'r', 'lz', 'sig', 'nc', 'dt', 'bdump',
                  'adump', 'ens'),
        'ensemble_long': ('N', 'epsilon', 'r', 'lz', 'sig', 'nc', 'dt',
                          'bdump', 'adump'),
        'ensemble': ('N', 'D', 'ac', 'nc'),
        'space': ('N', 'D', 'ac'),
    }
    _defaults = {'dmon': 1}
    _species = (Species('m', 'nmon', 'dmon'), Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                    'rho_bulk_c'],
        'whole': ['dmon',  'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                  'rho_bulk_c'],
        'ensemble_long': ['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                          'rho_bulk_c'],
        'ensemble': ['dmon'],
        'space': ['dmon']
    }


class SumRuleCubHeteroRing(ParserBase):
    """
    Parse an artifact name of the *SumRuleCubHeteroRing* project: a ring of
    small and large monomers among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
    - ``whole``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#ac#nc#
    - ``space``: ns#nl#al#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the ring, ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer; keyword ``ml``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    mcrowd : float
        Crowder mass, :math:`m_c = d_c^3`, at ``segment``, ``whole``, and
        ``ensemble_long``; ``-1``, a placeholder, at ``ensemble`` and
        ``space``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    The mass density is uniform across all species, so a species whose mass
    is not in the name has the mass :math:`m_i = d_i^3`, where :math:`d_i` is
    its diameter.

    Examples
    --------
    >>> artifact = SumRuleCubHeteroRing('ns400nl5al5.0ac1.0', 'space', 'bug')
    >>> print(artifact.dcrowd)
    1.0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                    'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                  'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt',
                          'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'ac'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1, 'mcrowd': -1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c', 'ncrowd', 'dcrowd'))
    _totals = {'m': ('m_small', 'm_large')}
    _masses = {'mcrowd': 'dcrowd'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                    'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }


class SumRuleCubHeteroLinear(ParserBase):
    """
    Parse an artifact name of the *SumRuleCubHeteroLinear* project: a linear
    chain of small and large monomers among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.linear
    - ``whole``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.linear
    - ``ensemble_long``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#ac#nc#
    - ``space``: ns#nl#al#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the chain, ``'all'`` for every
        particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer; keyword ``ml``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    mcrowd : float
        Crowder mass, :math:`m_c = d_c^3`, at ``segment``, ``whole``, and
        ``ensemble_long``; ``-1``, a placeholder, at ``ensemble`` and
        ``space``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    The mass density is uniform across all species, so a species whose mass
    is not in the name has the mass :math:`m_i = d_i^3`, where :math:`d_i` is
    its diameter.

    Examples
    --------
    >>> artifact = SumRuleCubHeteroLinear('ns800nl5al6.0ac3.0', 'space', 'all')
    >>> print(artifact.dmon_large)
    6.0
    """
    _geometry = 'cubic'
    _topology = 'linear'
    _groups = ['bug', 'all']
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                    'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                  'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt',
                          'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'ac'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1, 'mcrowd': -1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c', 'ncrowd', 'dcrowd'))
    _totals = {'m': ('m_small', 'm_large')}
    _masses = {'mcrowd': 'dcrowd'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                    'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }


class SumRuleCubHeteroRingBiDisperse(ParserBase):
    """
    Parse an artifact name of the *SumRuleCubHeteroRingBiDisperse* project: a
    ring of small and large monomers among crowders of two sizes in a
    periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: al#nl#ns#ac#ab#pctc#nc#nb#l#dt#bdump#adump#ens#.j#.ring
    - ``whole``: al#nl#ns#ac#ab#pctc#nc#nb#l#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: al#nl#ns#ac#ab#pctc#nc#nb#l#dt#bdump#adump#
    - ``ensemble``: al#nl#ns#ac#ab#pctc#nc#nb#
    - ``space``: al#nl#ns#ac#ab#pctc#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the ring, ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer, :math:`m_l = d_l^3`; only at ``segment``,
        ``whole``, and ``ensemble_long``.
    dcrowd_one, dcrowd_two : float
        Diameters of the type-one and type-two crowders; keywords ``ac`` and
        ``ab``.
    ncrowd_one, ncrowd_two : int
        Numbers of type-one and type-two crowders; keywords ``nc`` and
        ``nb``.
    mcrowd_one, mcrowd_two : float
        Masses of the type-one and type-two crowders, :math:`m = d^3`; only
        at ``segment``, ``whole``, and ``ensemble_long``.
    pct_one : float
        Percentage of the total crowder volume fraction
        :math:`\\phi_c = \\phi_{c,1} + \\phi_{c,2}` that the type-one crowders
        carry, :math:`100\\,\\phi_{c,1}/\\phi_c`; keyword ``pctc``.
    pct_two : float
        Percentage that the type-two crowders carry, ``100 - pct_one``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c_one, phi_bulk_c_one : float
        Bulk number density and volume fraction of type-one crowders.
    rho_bulk_c_two, phi_bulk_c_two : float
        Bulk number density and volume fraction of type-two crowders.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of all crowders: the sums of
        the type-one and type-two values.

    Notes
    -----
    The bulk values and the computed masses exist only at ``segment``,
    ``whole``, and ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    A space is a line of state points at fixed `pct_one` along which the
    total crowder volume fraction :math:`\\phi_c` varies; both crowder
    numbers change along it, so the ``space`` name drops ``nc`` and ``nb``.

    The mass density is uniform across all species, so a species whose mass
    is not in the name has the mass :math:`m_i = d_i^3`, where :math:`d_i` is
    its diameter.

    Every pair of particle types interacts through a Lennard-Jones potential
    of strength :math:`\\epsilon = 1`, the energy unit, and there is no wall:
    the project's LAMMPS input sets :math:`\\epsilon = 1` for all ten pairs of
    its four particle types.

    Examples
    --------
    >>> artifact = SumRuleCubHeteroRingBiDisperse(
    ...     'al4.0nl5ns400ac1.0ab5.0pctc40.0', 'space', 'bug')
    >>> print(artifact.dcrowd_two, artifact.pct_one, artifact.pct_two)
    5.0 40.0 60.0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd_one', float),
        'ab': KeywordSpec('dcrowd_two', float),
        'pctc': KeywordSpec('pct_one', float),
        'nc': KeywordSpec('ncrowd_one', int),
        'nb': KeywordSpec('ncrowd_two', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ns', 'ac', 'ab', 'pctc', 'nc', 'nb', 'l',
                    'dt', 'bdump', 'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ns', 'ac', 'ab', 'pctc', 'nc', 'nb', 'l', 'dt',
                  'bdump', 'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ns', 'ac', 'ab', 'pctc', 'nc', 'nb',
                          'l', 'dt', 'bdump', 'adump'),
        'ensemble': ('al', 'nl', 'ns', 'ac', 'ab', 'pctc', 'nc', 'nb'),
        'space': ('al', 'nl', 'ns', 'ac', 'ab', 'pctc'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c_one', 'ncrowd_one', 'dcrowd_one'),
                Species('c_two', 'ncrowd_two', 'dcrowd_two'))
    _totals = {'m': ('m_small', 'm_large'), 'c': ('c_one', 'c_two')}
    _masses = {'mmon_large': 'dmon_large', 'mcrowd_one': 'dcrowd_one',
               'mcrowd_two': 'dcrowd_two'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mmon_large', 'mcrowd_one',
                    'mcrowd_two', 'pct_two', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large',
                    'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                    'phi_bulk_c_one', 'rho_bulk_c_one', 'phi_bulk_c_two',
                    'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mmon_large', 'mcrowd_one',
                  'mcrowd_two', 'pct_two', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c_one',
                  'rho_bulk_c_one', 'phi_bulk_c_two', 'rho_bulk_c_two',
                  'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small', 'mmon_large',
                          'mcrowd_one', 'mcrowd_two', 'pct_two',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large',
                          'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c_one',
                          'rho_bulk_c_one', 'phi_bulk_c_two',
                          'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'pct_two'],
        'space': ['dmon_small', 'mmon_small', 'pct_two'],
    }

    def _derived_attributes(self) -> None:
        """
        Set `pct_two`, the percentage of the total crowder volume fraction
        that the type-two crowders carry: ``100 - pct_one``.
        """
        self.pct_two: float = 100 - getattr(self, 'pct_one')


class TransFociCyl(ParserBase):
    """
    Parse an artifact name of the *TransFociCyl* project: a ring of small and
    large monomers among crowders in a cylinder that is periodic along its
    axis.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.j#.ring
    - ``whole``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#D#ac#nc#
    - ``space``: ns#nl#al#D#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the ring, ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer; keyword ``ml``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    mcrowd : float
        Crowder mass, :math:`m_c = d_c^3`, at ``segment``, ``whole``, and
        ``ensemble_long``; ``-1``, a placeholder, at ``ensemble`` and
        ``space``.
    lcyl : float
        Length of the cylinder, :math:`Z = 2 l_z`, from the half-length
        :math:`l_z` (keyword ``lz``).
    dcyl : float
        Size of the cylinder, :math:`D = 2r - 1`, from the radius :math:`r`
        (keyword ``r``) in the long names; :math:`D` itself (keyword ``D``)
        in the short ones.
    epsilon_small : float
        Strength of the wall-small-monomer LJ interaction; keyword ``epss``.
    epsilon_large : float
        Strength of the wall-large-monomer LJ interaction; keyword ``epsl``.
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    - The bulk values exist only at ``segment``, ``whole``, and
      ``ensemble_long``. The lineage names are documented in
      :class:`ParserBase`.

    - The mass density is uniform across all species, so a species whose
      mass is not in the name has the mass :math:`m_i = d_i^3`, where
      :math:`d_i` is its diameter.

    - The cylindrical wall is a layer of beads of size :math:`a_w = 1` whose
      centres lie on the cylinder of radius :math:`r` set in LAMMPS (keyword
      ``r``), so the cylinder size is :math:`D = 2r - a_w`, stored as `dcyl`.
      The keyword ``lz`` is the half-length of the cylinder, so its length is
      :math:`Z = 2 l_z`, stored as `lcyl`. The centres of particles of
      species :math:`\\alpha` with size :math:`a_\\alpha` are confined to a
      cylinder of size :math:`D_{\\mathrm{eff},\\alpha} = D - a_\\alpha`, and
      the bulk number density and volume fraction of that species both use
      the effective volume :math:`V_{\\mathrm{eff},\\alpha} = \\pi Z
      D_{\\mathrm{eff},\\alpha}^2 / 4`.

    Examples
    --------
    >>> artifact = TransFociCyl(
    ... 'ns500nl5al3.0D20.0ac2.0nc0',
    ... 'ensemble',
    ... 'bug'
    ... )
    >>> print(artifact.dcyl)
    20.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _keywords = {
        'epss': KeywordSpec('epsilon_small', float),
        'epsl': KeywordSpec('epsilon_large', float),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac', 'nc',
                    'lz', 'dt', 'bdump', 'adump', 'ens', 'j'),
        'whole': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac', 'nc',
                  'lz', 'dt', 'bdump', 'adump', 'ens'),
        'ensemble_long': ('epss', 'epsl', 'r', 'al', 'nl', 'ml', 'ns', 'ac',
                          'nc', 'lz', 'dt', 'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'D', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'D', 'ac'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1, 'mcrowd': -1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c', 'ncrowd', 'dcrowd'))
    _totals = {'m': ('m_small', 'm_large')}
    _masses = {'mcrowd': 'dcrowd'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                    'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }


class TransFociCub(ParserBase):
    """
    Parse an artifact name of the *TransFociCub* project: a ring of small and
    large monomers among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
    - ``whole``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
    - ``ensemble_long``: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#
    - ``ensemble``: ns#nl#al#ac#nc#
    - ``space``: ns#nl#al#ac#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'bug', 'all'}
        Particle group: ``'bug'`` for the ring, ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon_small : int
        Diameter of a small monomer, fixed at 1; not in the name.
    nmon_small : int
        Number of small monomers; keyword ``ns``.
    mmon_small : int
        Mass of a small monomer, fixed at 1; not in the name.
    dmon_large : float
        Diameter of a large monomer; keyword ``al``.
    nmon_large : int
        Number of large monomers; keyword ``nl``.
    mmon_large : float
        Mass of a large monomer; keyword ``ml``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    mcrowd : float
        Crowder mass, :math:`m_c = d_c^3`, at ``segment``, ``whole``, and
        ``ensemble_long``; ``-1``, a placeholder, at ``ensemble`` and
        ``space``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step; keyword ``dt``.
    bdump : int
        Dump period of the ``bug`` trajectory, in time steps; keyword
        ``bdump``.
    adump : int
        Dump period of the ``all`` trajectory, in time steps; keyword
        ``adump``.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m_small, phi_bulk_m_small : float
        Bulk number density and volume fraction of small monomers.
    rho_bulk_m_large, phi_bulk_m_large : float
        Bulk number density and volume fraction of large monomers.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of all monomers: the sums of
        the small and large values.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    The mass density is uniform across all species, so a species whose mass
    is not in the name has the mass :math:`m_i = d_i^3`, where :math:`d_i` is
    its diameter.

    Examples
    --------
    >>> artifact = TransFociCub(
    ... 'al5.0nl5ml125.0ns400ac1.0nc0l100.0dt0.005bdump2000adump5000.ring',
    ... 'ensemble_long',
    ... 'all'
    ... )
    >>> print(artifact.ncrowd)
    0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']
    _keywords = {
        'al': KeywordSpec('dmon_large', float),
        'nl': KeywordSpec('nmon_large', int),
        'ml': KeywordSpec('mmon_large', float),
        'ns': KeywordSpec('nmon_small', int),
        'ac': KeywordSpec('dcrowd', float),
        'nc': KeywordSpec('ncrowd', int),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'dt': KeywordSpec('dt', float),
        'bdump': KeywordSpec('bdump', int),
        'adump': KeywordSpec('adump', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                    'adump', 'ens', 'j'),
        'whole': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt', 'bdump',
                  'adump', 'ens'),
        'ensemble_long': ('al', 'nl', 'ml', 'ns', 'ac', 'nc', 'l', 'dt',
                          'bdump', 'adump'),
        'ensemble': ('ns', 'nl', 'al', 'ac', 'nc'),
        'space': ('ns', 'nl', 'al', 'ac'),
    }
    _defaults = {'dmon_small': 1, 'mmon_small': 1, 'mcrowd': -1}
    _species = (Species('m_small', 'nmon_small', 'dmon_small'),
                Species('m_large', 'nmon_large', 'dmon_large'),
                Species('c', 'ncrowd', 'dcrowd'))
    _totals = {'m': ('m_small', 'm_large')}
    _masses = {'mcrowd': 'dcrowd'}
    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                    'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'rho_bulk_m_large',
                  'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'rho_bulk_m_large', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }


class HnsCub(ParserBase):
    """
    Parse an artifact name of the *HnsCub* project: a ring polymer with H-NS
    proteins among crowders in a periodic cube.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: N#kbmm#nh#ac#l#epshc#nc#ens#.j#.ring
    - ``whole``: N#kbmm#nh#ac#l#epshc#nc#ens#.ring
    - ``ensemble_long``: N#kbmm#nh#ac#l#epshc#nc#
    - ``ensemble``: N#kbmm#nh#ac#epshc#nc#
    - ``space``: N#kbmm#nh#ac#epshc#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'nucleoid', 'all'}
        Particle group: ``'nucleoid'`` or ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon : int
        Monomer diameter, fixed at 1; not in the name.
    nmon : int
        Number of monomers; keyword ``N``.
    dhns : int
        Diameter of an H-NS protein, fixed at 1; not in the name.
    dhns_patch : float
        Diameter of the patch at each pole of an H-NS protein, fixed at
        0.178; not in the name.
    nhns : int
        Number of H-NS proteins; keyword ``nh``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    bend_mm : float
        Bending rigidity of the DNA monomers; keyword ``kbmm``.
    eps_hm : int
        Strength of the attractive LJ interaction between the H-NS poles and
        the monomers, fixed at 29; not in the name.
    eps_hc : float
        Strength of the attractive LJ interaction between the H-NS cores and
        the crowders; keyword ``epshc``.
    lcube : float
        Side of the cube, :math:`L = 2 l`, from the half-side :math:`l`
        (keyword ``l``).
    dt : float
        Simulation time step, fixed at 0.005; not in the name.
    ndump : int
        Dump period of the ``nucleoid`` trajectory, fixed at 5000 time steps;
        not in the name.
    adump : int
        Dump period of the ``all`` trajectory, fixed at 10000 time steps; not
        in the name.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of monomers.
    rho_bulk_hns, phi_bulk_hns : float
        Bulk number density and volume fraction of H-NS proteins.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    Examples
    --------
    >>> artifact = HnsCub(
    ... 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
    ... 'whole',
    ... 'nucleoid'
    ... )
    >>> print(artifact.nhns)
    8
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'kbmm': KeywordSpec('bend_mm', float),
        'nh': KeywordSpec('nhns', int),
        'ac': KeywordSpec('dcrowd', float),
        'l': KeywordSpec('lcube', float, scale=2.0),
        'epshc': KeywordSpec('eps_hc', float),
        'nc': KeywordSpec('ncrowd', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc', 'ens', 'j'),
        'whole': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc', 'ens'),
        'ensemble_long': ('N', 'kbmm', 'nh', 'ac', 'l', 'epshc', 'nc'),
        'ensemble': ('N', 'kbmm', 'nh', 'ac', 'epshc', 'nc'),
        'space': ('N', 'kbmm', 'nh', 'ac', 'epshc'),
    }
    _defaults = {'dmon': 1, 'dhns': 1, 'dhns_patch': 0.178, 'dt': 0.005,
                 'ndump': 5000, 'adump': 10000, 'eps_hm': 29}
    _species = (Species('m', 'nmon', 'dmon'), Species('hns', 'nhns', 'dhns'),
                Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                    'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                    'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'whole': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                  'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt', 'ndump',
                  'adump', 'eps_hm'],
        'ensemble_long': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                          'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'ensemble': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'space': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']
    }


class HnsCyl(ParserBase):
    """
    Parse an artifact name of the *HnsCyl* project: a ring polymer with H-NS
    proteins among crowders in a cylinder that is periodic along its axis.

    Each lineage has its own name pattern; ``#`` stands for a value:

    - ``segment``: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.j#.ring
    - ``whole``: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.ring
    - ``ensemble_long``: N#kbmm#r#nh#ac#lz#epshc#nc#
    - ``ensemble``: N#D#nh#ac#epshc#nc#
    - ``space``: N#D#nh#ac#epshc#

    Parameters
    ----------
    artifact : Pathish
        Name to be parsed, or a file path with ``ispath=True``.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Lineage of the artifact.
    group : {'nucleoid', 'all'}
        Particle group: ``'nucleoid'`` or ``'all'`` for every particle.
    ispath : bool, default False
        Whether `artifact` is a file path.

    Attributes
    ----------
    dmon : int
        Monomer diameter, fixed at 1; not in the name.
    nmon : int
        Number of monomers; keyword ``N``.
    dhns : int
        Diameter of an H-NS protein, fixed at 1; not in the name.
    dhns_patch : float
        Diameter of the patch at each pole of an H-NS protein, fixed at
        0.178; not in the name.
    nhns : int
        Number of H-NS proteins; keyword ``nh``.
    dcrowd : float
        Crowder diameter; keyword ``ac``.
    ncrowd : int
        Number of crowders; keyword ``nc``.
    bend_mm : float
        Bending rigidity of the DNA monomers; keyword ``kbmm``.
    eps_hm : int
        Strength of the attractive LJ interaction between the H-NS poles and
        the monomers, fixed at 29; not in the name.
    eps_hc : float
        Strength of the attractive LJ interaction between the H-NS cores and
        the crowders; keyword ``epshc``.
    lcyl : float
        Length of the cylinder, :math:`Z = 2 l_z`, from the half-length
        :math:`l_z` (keyword ``lz``).
    dcyl : float
        Size of the cylinder, :math:`D = 2r - 1`, from the radius :math:`r`
        (keyword ``r``) in the long names; :math:`D` itself (keyword ``D``)
        in the short ones.
    dt : float
        Simulation time step, fixed at 0.005; not in the name.
    ndump : int
        Dump period of the ``nucleoid`` trajectory, fixed at 5000 time steps;
        not in the name.
    adump : int
        Dump period of the ``all`` trajectory, fixed at 10000 time steps; not
        in the name.
    ensemble_id : int
        Index of the whole in its ensemble; keyword ``ens``.
    segment_id : int
        Index of the segment in its whole, zero-padded in the name; keyword
        ``j``.
    rho_bulk_m, phi_bulk_m : float
        Bulk number density and volume fraction of monomers.
    rho_bulk_hns, phi_bulk_hns : float
        Bulk number density and volume fraction of H-NS proteins.
    rho_bulk_c, phi_bulk_c : float
        Bulk number density and volume fraction of crowders.

    Notes
    -----
    The bulk values exist only at ``segment``, ``whole``, and
    ``ensemble_long``. The lineage names are documented in
    :class:`ParserBase`.

    The cylindrical wall is a layer of beads of size :math:`a_w = 1` whose
    centres lie on the cylinder of radius :math:`r` set in LAMMPS (keyword
    ``r``), so the cylinder size is :math:`D = 2r - a_w`, stored as `dcyl`.
    The keyword ``lz`` is the half-length of the cylinder, so its length is
    :math:`Z = 2 l_z`, stored as `lcyl`. The centres of particles of
    species :math:`\\alpha` with size :math:`a_\\alpha` are confined to a
    cylinder of size :math:`D_{\\mathrm{eff},\\alpha} = D - a_\\alpha`, and
    the bulk number density and volume fraction of that species both use
    the effective volume :math:`V_{\\mathrm{eff},\\alpha} = \\pi Z
    D_{\\mathrm{eff},\\alpha}^2 / 4`.

    Examples
    --------
    >>> artifact = HnsCyl(
    ... 'N200D20.0nh16ac1.0epshc1.0nc0',
    ... 'space',
    ... 'nucleoid'
    ... )
    >>> print(artifact.eps_hc)
    1.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']
    _keywords = {
        'N': KeywordSpec('nmon', int),
        'kbmm': KeywordSpec('bend_mm', float),
        'r': KeywordSpec('dcyl', float, scale=2.0, offset=-1.0),
        'D': KeywordSpec('dcyl', float),
        'nh': KeywordSpec('nhns', int),
        'ac': KeywordSpec('dcrowd', float),
        'lz': KeywordSpec('lcyl', float, scale=2.0),
        'epshc': KeywordSpec('eps_hc', float),
        'nc': KeywordSpec('ncrowd', int),
        'ens': KeywordSpec('ensemble_id', int),
        'j': KeywordSpec('segment_id', int),
    }
    _patterns = {
        'segment': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc', 'ens',
                    'j'),
        'whole': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc', 'ens'),
        'ensemble_long': ('N', 'kbmm', 'r', 'nh', 'ac', 'lz', 'epshc', 'nc'),
        'ensemble': ('N', 'D', 'nh', 'ac', 'epshc', 'nc'),
        'space': ('N', 'D', 'nh', 'ac', 'epshc'),
    }
    _defaults = {'dmon': 1, 'dhns': 1, 'dhns_patch': 0.178, 'dt': 0.005,
                 'ndump': 5000, 'adump': 10000, 'eps_hm': 29}
    _species = (Species('m', 'nmon', 'dmon'), Species('hns', 'nhns', 'dhns'),
                Species('c', 'ncrowd', 'dcrowd'))
    _project_attributes = {
        'segment': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                    'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                    'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'whole': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                  'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt', 'ndump',
                  'adump', 'eps_hm'],
        'ensemble_long': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m',
                          'phi_bulk_c', 'rho_bulk_c', 'phi_bulk_hns',
                          'rho_bulk_hns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'ensemble': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm'],
        'space': ['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']
    }
