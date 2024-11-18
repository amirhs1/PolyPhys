"""\
==========================================================
:mod:`polyphys.manage.parser`
==========================================================

This :mod:`~polyphys.manage.parser` module defines a suite of base and
specialized parser classes for extracting structured, lineage-based information
from filenames and file paths in molecular simulation projects. Each parser
subclass represents a unique project type and geometry (e.g., cubic,
cylindrical), allowing lineage-specific parsing of artifact attributes from
filenames.

The lineage attribute hierarchy includes the following levels: 'segment',
'whole', 'ensemble_long', 'ensemble', and 'space'. Each level represents a
different scope within the project, ranging from individual simulation segments
to complete simulation spaces. Different classes support parsing of project
attributes such as particle dimensions, densities, and simulation box
dimensions. These classes employ various methods to dynamically interpret and
calculate attributes based on parsed lineage and project requirements.

Classes
=======
.. autoclass:: ParserBase
   :members:
   :undoc-members:
.. autoclass:: TwoMonDep
.. autoclass:: SumRuleCyl
.. autoclass:: SumRuleCubHeteroRing
.. autoclass:: SumRuleCubHeteroLinear
.. autoclass:: TransFociCyl
.. autoclass:: TransFociCub
.. autoclass:: HnsCub
.. autoclass:: HnsCyl

Dependencies
============
- `os`: For handling file paths.
- `re`: For regular expressions in filename parsing.
- `typing`: For type hinting.
- `abc`: For defining abstract base classes.
- `collections`: For ordered dictionary functionality.
- `.utilizer`: Utility functions, e.g., `invalid_keyword`.
- `..analyze.measurer`: Measurement functions such as `number_density_cube` and
  `volume_fraction_cube`.

Usage
=====
Classes in this module are instantiated with artifact filenames, lineages, and
group types, which they parse into a rich set of attributes. These classes
facilitate programmatic access to file-based information for complex
projects.

Examples
========
>>> artifact = SumRuleCyl("N200D10.0ac2.0nc0", 'ensemble', 'all')
>>> print(artifact.dcrowd)
2.0

>>> artifact = HnsCyl("N200D20.0nh16ac1.0epshc1.0nc0", 'space', 'nucleoid')
>>> print(artifact.eps_hc)
1.0

Notes
=====
These classes use dynamic attribute setting based on parsed filename data.
To handle attributes that may or may not exist at runtime, the `getattr`
function is commonly used, ensuring robustness to varied filename patterns.
"""
import os
import re
from typing import Dict, List, Literal, ClassVar, Optional
from abc import ABC, abstractmethod
from collections import OrderedDict
from .utilizer import invalid_keyword
from ..analyze.measurer import (
    number_density_cube,
    number_density_cylinder,
    volume_fraction_cube,
    volume_fraction_cylinder
)
from .typer import LineageT, TopologyT, GroupT, GeometryT


class ParserBase(ABC):
    """
    Base parser class for extracting information from filenames or file paths
    in a structured project. Designed to enforce lineage, geometry, and group
    conventions across subclasses for specific project types.

    Parameters
    ----------
    artifact : str
        Artifact that is parsed for extracting information. Can be a filename
        or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        The lineage of the name, specifying the hierarchical level within the
        project.
    group : str
        The particle group type in the project. Used for specific group-based
        parsing.

    Attributes
    ----------
    filepath : str
        The filepath if `artifact` is a filepath, otherwise "N/A".
    filename : str
        The filename extracted from `artifact` if it is a filepath, otherwise
        `artifact` itself.
    group : {'bug', 'nucleoid', 'all'}
        Particle group type in the project.
    name : str
        The unique name derived from `filename` based on `lineage` and `group`
        conventions.
    project_name : str
        The name of the project class (subclass of `ParserBase`), automatically
        assigned as the subclass's name.
    lineage_genealogy: List[str]
        List of parent lieages for an `artifact` with a given `lineage`.
    lineage_attributes: List[str]
        List of parsed/dynamically-defined attributes specific to an `artifact`
        with a given `lineage`,
    physical_attributes: List[str]
        List of computed/system attributes specific to an `artifact` with a
        given `lineage`, including computed attributes.
    attributes: List[str]
        List of attributes specific to an `artifact` with a given `lineage`,
        including parsed and computed attributes.

    Class Attributes
    ----------------
    lineages : list of str
        List of valid lineage types.
    genealogy : dict of lists
        Dictionary defining the hierarchical relationship between lineages.
        Each lineage points to its parent lineages in a hierarchy.

    Abstract Class Properties
    -------------------------
    geometry : str
        Specifies geometry of the system
    topology : str
        Specifies how particles are connected (how) or not in the system.
    groups : List[str]
        Specifies particle groups in the system
    genealogy_attributes : Dict[str, OrderedDict[str, str]]
        Dictionary defining lineage-specific attributes for each lineage type.
    project_attributes : Dict[str, List[str]]
        Dictionary defining project-level attributes that remain constant but
        are not extractable from a filename.

    Methods
    -------
    _find_name() -> None
        Parses and sets the unique name based on `lineage` and `group`.
    _set_parents() -> None
        Sets pattern names for each `lineage` based on `_genealogy`.
    _initiate_attributes() -> None
        Defines and initializes subclass-specific attributes. (Abstract method)
    _parse_name() -> None
        Parses lineage-specific attributes based on the filename.
        (Abstract method)
    _bulk_attributes() -> None
        Computes physical attributes for the current lineage based on primary
        attributes. (Abstract method)
    """
    _lineages: ClassVar[List[LineageT]] = \
        ['segment', 'whole', 'ensemble_long', 'ensemble', 'space']
    _genealogy: ClassVar[Dict[LineageT, List[LineageT]]] = {
        'segment': ['segment', 'whole', 'ensemble_long', 'ensemble', 'space'],
        'whole': ['whole', 'ensemble_long', 'ensemble', 'space'],
        'ensemble_long': ['ensemble_long', 'ensemble', 'space'],
        'ensemble': ['ensemble', 'space'],
        'space': ['space'],
    }
    _geometry: ClassVar[Optional[GeometryT]] = None
    _topology: ClassVar[Optional[TopologyT]] = None
    _groups: ClassVar[Optional[List[GroupT]]] = None
    _genealogy_attributes: ClassVar[Dict[LineageT, OrderedDict[str, str]]] = \
        None
    _project_attributes: ClassVar[Dict[str, List[str]]] = None

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: GroupT
    ) -> None:
        self._filepath, self._filename = (
            ("N/A", artifact) if '/' not in artifact and '\\' not in artifact
            else os.path.split(artifact)
        )
        invalid_keyword(lineage, self._lineages)
        invalid_keyword(group, self._groups)
        self._lineage = lineage
        self._group = group
        self._project_name = self.__class__.__name__
        self._lineage_genealogy: List[LineageT] = self._genealogy[lineage]
        self._lineage_attributes = \
            list(self._genealogy_attributes[lineage].keys())
        self._physical_attributes = self._project_attributes[lineage]
        self._attributes = \
            self._lineage_attributes + self._physical_attributes
        self._find_name()

    def __str__(self) -> str:
        """
        Provides a formatted summary of the parser instance.
        """
        observation = (
            f"Arifact:\n"
            f"    Name: '{self.filename}',\n"
            f"    Geometry: '{self._geometry}',\n"
            f"    Group: '{self._group}',\n"
            f"    Lineage: '{self._lineage}',\n"
            f"    Topology: '{self._topology}'"
        )
        return observation

    def __repr__(self) -> str:
        return (
            f"Artifact('{self.filename}' in geometry"
            f" '{self._geometry}' from group '{self._group}' with"
            f" lineage '{self._lineage}' and"
            f" topology '{self._topology}')"
        )

    @property
    def geometry(self) -> GeometryT:
        """
        System geometry in a molecular dynamics system.
        """
        if self._geometry is None:
            raise AttributeError("'_geometry' has not been initialized.")
        return self._geometry

    @property
    def groups(self) -> List[GroupT]:
        """
        List of valid group names for the subclass.
        """
        if self._groups is None:
            raise AttributeError("'_groups' has not been initialized.")
        return self._groups

    @property
    def topology(self) -> TopologyT:
        """
        Defines the polymer topology for the parser subclass.
        """
        if self._topology is None:
            raise AttributeError("'_topology' has not been initialized.")
        return self._topology

    @property
    def genealogy_attributes(self) -> Dict[LineageT, OrderedDict[str, str]]:
        """
        Dictionary of lineage-specific attributes. Each key is a lineage type,
        and each value is an OrderedDict mapping attribute names to their
        short-form representations.
        """
        if self._topology is None:
            raise AttributeError(
                "'_genealogy_attributes' has not been initialized.")
        return self._genealogy_attributes

    @property
    def project_attributes(self) -> Dict[str, List[str]]:
        """
        Dictionary of project attributes. Each key is a lineage type,
        and each value is an OrderedDict mapping attribute names to their
        short-form representations.
        """
        if self._topology is None:
            raise AttributeError(
                "'_project_attributes' has not been initialized.")
        return self._project_attributes

    @property
    def filename(self) -> str:
        """
        Returns the filename, either extracted from the path or the name
        itself.
        """
        return self._filename

    @property
    def filepath(self) -> str:
        """
        Returns the full filepath or 'N/A' if not a valid path.
        """
        return self._filepath

    @property
    def group(self) -> GroupT:
        """
        Returns the current group.
        """
        return self._group

    @property
    def lineage(self) -> LineageT:
        """
        Returns the current lineage.
        """
        return self._lineage

    @property
    def name(self) -> str:
        """
        Returns the unique name parsed from the filename.
        """
        return self._name

    @property
    def project_name(self) -> str:
        """
        Returns the project (parser class) name,
        """
        return self._project_name

    @property
    def attributes(self) -> List[str]:
        """
        Returns lineage-specific andp project attributes for an artifact.
        """
        return self._attributes

    @property
    def lineage_genealogy(self) -> List[LineageT]:
        """
        Returns the parents of a given `lineage`.
        """
        return self._lineage_genealogy

    @property
    def lineage_attributes(self) -> List[str]:
        """
        Returns lineage-specific attributes for an artifact with a given
        `lineage`.
        """
        return self._lineage_attributes

    @property
    def physical_attributes(self) -> List[str]:
        """
        Returns project-level attributes for an artifact with a given
        `lineage`.
        """
        return self._physical_attributes

    def _find_name(self) -> None:
        """
        Parses and sets the unique `lineage_name` from the filename
        based on the `lineage` and `group`.

        Notes
        -----
        - For 'segment' and 'whole' lineages, names typically end with the
          group keyword or a hyphen.
        - For 'ensemble_long', 'ensemble', and 'space', names are derived
          from the first substring in the filename.
        """
        if self._lineage in ['segment', 'whole']:
            self._name = \
                self.filename.split("." + self._group)[0].split("-")[0]
        else:
            self._name = self._filename.split("-")[0]

    @abstractmethod
    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes. Lineage attributes are
        set dynamically via `_parse_name` method.
        """

    def _set_parents(self) -> None:
        """
        Sets parent lineage names for each lineage type, following the
        hierarchy defined in `_genealogy`.

        Notes
        -----
        The `self._genealogy` defines the following parent-child hierarchy:

            - space -> ensemble -> ensemble_long -> whole -> segment

        Each lineage on the left has all the lineages on its right.
        """
        for lineage_name in self._lineages:
            lineage_value = "N/A"
            if lineage_name in self._genealogy[self._lineage]:
                lineage_value = ""
                lineage_attr = self._genealogy_attributes[lineage_name]
                for attr_long, attr_short in lineage_attr.items():
                    lineage_value += \
                            f"{attr_short}{getattr(self, attr_long)}"
            setattr(self, lineage_name, lineage_value)

    @abstractmethod
    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """

    @abstractmethod
    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """


class TwoMonDepCub(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *TwoMonDep* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#.j#
      One of multiple chunks of a complete artifact.
    - `whole`: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: nm#am#ac#nc#sd#
      Short name for an 'ensemble' artifact.
    - `space`: nm#am#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon : float
        Size (diameter) of a monomer. Its associated keyword is 'am'.
    nmon : int
        Number of monomers. Its associated keyword is 'N'.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    lcube : float
        Length of the simulation box, inferred from 'hl' keyword
        (half-length of hthe simulation box).
    d_sur : float
        Surface-to-surface distance between two monomers fixed in space.
        Its associated keyword is 'sd'.
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    tdump : int
        Frequency by which 'themo' variables are written in a 'lammps'
        log file. Its associated keyword is 'tdump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = (
    ..."nm2am5.0ac1.0"
    ...'space',
    ...'bug'
    ... )
    >>> print(artifact.nmon)
    2
    """
    _geometry = 'cubic'
    _topology = 'atomic'
    _groups = ['bug', 'all']

    _genealogy_attributes = {
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#.j#
        'segment': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'hl', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump', 'ensemble_id': 'ens',
            'segment_id': 'j'}
            ),
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#
        'whole': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'hl', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump', 'ensemble_id': 'ens'}
            ),
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump# :
        'ensemble_long': OrderedDict({
            'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
            'lcube': 'hl', 'd_sur': 'sd', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'tdump': 'tdump'}
            ),
        # Pattern: nm#am#ac#nc#sd :
        'ensemble': OrderedDict(
            {'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'd_sur': 'sd'}
             ),
        # pttern: nm#am#ac# :
        'space': OrderedDict(
            {'dmon': 'am', 'nmon': 'nm', 'dcrowd': 'ac', 'ncrowd': 'nc'}
            )
    }

    _project_attributes = {
        'segment': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                          'rho_bulk_c'],
        'ensemble': [],
        'space': []
        }

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon', 'lcube', 'dcrowd', 'dt', 'd_sur']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'hl':
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m = number_density_cube(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m = volume_fraction_cube(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_c = number_density_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_c = volume_fraction_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )


class SumRuleCyl(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *SumRuleCyl* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
      One of multiple chunks of a complete artifact.
    - `whole`: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: N#D#ac#nc#
      Short name for an 'ensemble' artifact.
    - `space`: N#D#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon : float
        Size (diameter) of a monomer. Its associated keyword is 'am'.
    nmon : int
        Number of monomers. Its associated keyword is 'N'.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    lcyl : float
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis).
    dcyl : float
        Size (or diameter) of the cylindrical confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement).
    epsilon: float
        Wall-particle LJ interaction strength. Its associated keyword is
        'epsilon' keyword.
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    The cylindrical wall is implemented in LAMMPS by using wall-forming
    particles of size 1.0. Thus, the actual size of the cylinder size
    (diameter), :math:`D`, is :math:`D=2r-1.0`,  :math:`r` is the radius of
    the cylindrical region defined in LAMMPS.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = SumRuleCyl(
    ..."N200D10.0ac2.0nc0"
    ...'ensemble',
    ...'all'
    ... )
    >>> print(artifact.dcrowd)
    2.0
    """
    _geometry = 'cylindrical'
    _topology = 'linear'
    _groups = ['bug', 'all']

    _genealogy_attributes = {
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
        'segment':  OrderedDict(
            {'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
             'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
             'adump': 'adump', 'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
        'whole':  OrderedDict(
            {'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
             'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
             'adump': 'adump', 'ensemble_id': 'ens'}),
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
        'ensemble_long':  OrderedDict(
            {'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
             'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
             'adump': 'adump'}),
        # Pattern: N#D#ac#nc#
        'ensemble':  OrderedDict(
            {'nmon': 'N', 'dcyl': 'D', 'dcrowd': 'ac', 'ncrowd': 'nc'}),
        # Pattern: N#D#ac#
        'space':  OrderedDict({'nmon': 'N', 'dcyl': 'D', 'dcrowd': 'ac'})
    }
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

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon: float = 1
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon', 'dcyl', 'lcyl', 'epsilon', 'dcrowd', 'dt']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'lz':
                    # Cylinder full length from its half-length
                    setattr(self, attr, 2 * getattr(self, attr))
                if keyword == 'r':
                    # Cylinder size is twice its radius, correcting for
                    # wall-forming particles with size 1
                    setattr(self, attr, 2 * getattr(self, attr) - 1.0)
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m = number_density_cylinder(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_m = volume_fraction_cylinder(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.rho_bulk_c = number_density_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_c = volume_fraction_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )


class SumRuleCubHeteroRing(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *SumRuleCubHeteroRing* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
      One of multiple chunks of a complete artifact.
    - `whole`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: ns#nl#al#ac#nc#
      Short name for an 'ensemble' artifact.
    - `space`: ns#nl#al#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon_small: float
        Size (diameter) of a monomer
    nmon_small: int
        number of small monomers. Its associated keyword is 'ns'.
    mmon_small: float
        Mass of a small monomer
    dmon_large: float
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int
        number of large monomers. Its associated keyword is 'nl'.
    mmon_large: float, default np.nan
        Mass of a large monomer. Its associated keyword is 'ml'.
    nmon: int
        Total number of monomers.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder.
    lcube : float
        Length of the simulation box, inferred from 'l' keyword
        (half-length of hthe simulation box).
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m_small : float
        Bulk number density fraction of small monomers.
    phi_bulk_m_small : float
        Bulk volume fraction of small monomers
    rho_bulk_m_large : float
        Bulk number density fraction of large monomers.
    phi_bulk_m_large : float
        Bulk volume fraction of large monomers
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    The mass density is uniform across all species. For any species whose mass
    is not explicitly parsed, the mass is defined as :math:`m_{i} = d_{i}^3`,
    where :math:`d_{i}` represents the species' diameter.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = SumRuleCubHeteroRing(
    ..."ns400nl5am5.0ac1.0"
    ...'space',
    ...'bug'
    ... )
    >>> print(artifact.dcrowd)
    1.0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']

    _genealogy_attributes = {
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
        'segment': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
        'whole': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcube': 'l',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
        'ensemble_long': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump'}),
        # Pattern: ns#nl#al#ac#nc#
        'ensemble': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac', 'ncrowd': 'nc'}),
        # Pattern: ns#nl#al#ac#
        'space': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac'})
    }

    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                  'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'phi_bulk_m_large', 'rho_bulk_m',
                          'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.mmon_small: float = self.dmon_small**3
        self.mcrowd: float = -1
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m_small: float = -1
            self.rho_bulk_m_small: float = -1
            self.phi_bulk_m_large: float = -1
            self.rho_bulk_m_large: float = -1
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon_large', 'lcube', 'mmon_large', 'dcrowd', 'dt']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'l':
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m_small = number_density_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_small = volume_fraction_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m_large = number_density_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_large = volume_fraction_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m = self.rho_bulk_m_small + self.rho_bulk_m_large
        self.phi_bulk_m = self.phi_bulk_m_small + self.phi_bulk_m_large

        self.mcrowd = getattr(self, 'dcrowd') ** 3
        self.rho_bulk_c = number_density_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_c = volume_fraction_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )


class SumRuleCubHeteroLinear(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *SumRuleCubHeteroRing* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.linear
      One of multiple chunks of a complete artifact.
    - `whole`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.linear
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.linear
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: ns#nl#al#ac#nc#
      Short name for an 'ensemble' artifact.
    - `space`: ns#nl#al#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon_small: float
        Size (diameter) of a monomer
    nmon_small: int
        number of small monomers. Its associated keyword is 'ns'.
    mmon_small: float
        Mass of a small monomer
    dmon_large: float
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int
        number of large monomers. Its associated keyword is 'nl'.
    mmon_large: float, default np.nan
        Mass of a large monomer. Its associated keyword is 'ml'.
    nmon: int
        Total number of monomers.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder.
    lcube : float
        Length of the simulation box, inferred from 'l' keyword
        (half-length of hthe simulation box).
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m_small : float
        Bulk number density fraction of small monomers.
    phi_bulk_m_small : float
        Bulk volume fraction of small monomers
    rho_bulk_m_large : float
        Bulk number density fraction of large monomers.
    phi_bulk_m_large : float
        Bulk volume fraction of large monomers
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    The mass density is uniform across all species. For any species whose mass
    is not explicitly parsed, the mass is defined as :math:`m_{i} = d_{i}^3`,
    where :math:`d_{i}` represents the species' diameter.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = SumRuleCubHeteroLinear(
    ..."ns800nl5am6.0ac3.0"
    ...'space',
    ...'all'
    ... )
    >>> print(artifact.dmon_large)
    6.0
    """
    _geometry = 'cubic'
    _topology = 'linear'
    _groups = ['bug', 'all']

    _genealogy_attributes = {
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
        'segment': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
        'whole': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcube': 'l',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
        'ensemble_long': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump'}),
        # Pattern: ns#nl#al#ac#nc#
        'ensemble': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac', 'ncrowd': 'nc'}),
        # Pattern: ns#nl#al#ac#
        'space': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac'})
    }

    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                  'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'phi_bulk_m_large', 'rho_bulk_m',
                          'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.mmon_small: float = self.dmon_small**3
        self.mcrowd: float = -1
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m_small: float = -1
            self.rho_bulk_m_small: float = -1
            self.phi_bulk_m_large: float = -1
            self.rho_bulk_m_large: float = -1
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon_large', 'lcube', 'mmon_large', 'dcrowd', 'dt']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'l':
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m_small = number_density_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_small = volume_fraction_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m_large = number_density_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_large = volume_fraction_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m = self.rho_bulk_m_small + self.rho_bulk_m_large
        self.phi_bulk_m = self.phi_bulk_m_small + self.phi_bulk_m_large

        self.mcrowd = getattr(self, 'dcrowd') ** 3
        self.rho_bulk_c = number_density_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_c = volume_fraction_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )


class TransFociCyl(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *TransFociCyl* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.j#.ring
      One of multiple chunks of a complete artifact.
    - `whole`: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.ring
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#.ring
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: ns#nl#al#D#ac#nc#
      Short name for an 'ensemble' artifact.
    - `space`: ns#nl#al#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bugg', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon_small: float
        Size (diameter) of a monomer
    nmon_small: int
        number of small monomers. Its associated keyword is 'ns'.
    mmon_small: float
        Mass of a small monomer
    dmon_large: float
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int
        number of large monomers. Its associated keyword is 'nl'.
    mmon_large: float, default np.nan
        Mass of a large monomer. Its associated keyword is 'ml'.
    nmon: int
        Total number of monomers.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder.
    lcyl : float
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis).
    dcyl : float
        Size (or diameter) of the cylindrical confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement).
    epsilon_s: float
        Wall-small-monomer LJ interaction strength. Its associated keyword is
        'epss' keyword.
    epsilon_l: float
        Wall-large-monomer LJ interaction strength. Its associated keyword is
        'espl' keyword.
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m_small : float
        Bulk number density fraction of small monomers.
    phi_bulk_m_small : float
        Bulk volume fraction of small monomers
    rho_bulk_m_large : float
        Bulk number density fraction of large monomers.
    phi_bulk_m_large : float
        Bulk volume fraction of large monomers
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    - The mass density is uniform across all species. For any species whose
      mass is not explicitly parsed, the mass is defined as
      :math:`m_{i} = d_{i}^3`, where :math:`d_{i}` represents the species'
      diameter.

    - The cylindrical wall is implemented in LAMMPS by using wall-forming
      particles of size 1.0. Thus, the actual size of the cylinder size
      (diameter), :math:`D`, is :math:`D=2r-1.0`,  :math:`r` is the radius of
      the cylindrical region defined in LAMMPS.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = TransFociCyl(
    ..."ns500nl5al3.0D20.0ac2.0nc0"
    ...'ensemble',
    ...'bug'
    ... )
    >>> print(artifact.dcyl)
    20.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['bug', 'all']

    _genealogy_attributes: Dict[str, OrderedDict[str, str]] = {
        # Pattern: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.j#.ring
        'segment': OrderedDict(
            {'epsilon_small': 'epss', 'epsilon_large': 'epsl', 'dcyl': 'r',
             'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcyl': 'lz',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.ring
        'whole': OrderedDict(
            {'epsilon_small': 'epss', 'epsilon_large': 'epsl', 'dcyl': 'r',
             'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcyl': 'lz',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens'}),
        # Pattern: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#.ring
        'ensemble_long': OrderedDict(
            {'epsilon_small': 'epss', 'epsilon_large': 'epsl', 'dcyl': 'r',
             'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcyl': 'lz',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump'}),
        # Pattern: ns#nl#al#D#ac#nc#
        'ensemble': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcyl': 'D', 'dcrowd': 'ac', 'ncrowd': 'nc'}),
        # Pattern: ns#nl#al#D#ac#
        'space': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac'})
    }

    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                  'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'phi_bulk_m_large', 'rho_bulk_m',
                          'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.mmon_small: float = self.dmon_small**3
        self.mcrowd: float = -1
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m_small: float = -1
            self.rho_bulk_m_small: float = -1
            self.phi_bulk_m_large: float = -1
            self.rho_bulk_m_large: float = -1
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon_large', 'dcyl', 'lcyl', 'epsilon_small',
                       'epsilon_large', 'mmon_large', 'dcrowd', 'dt']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'lz':
                    # Cylinder full length from its half-length
                    setattr(self, attr, 2 * getattr(self, attr))
                if keyword == 'r':
                    # Cylinder size is twice its radius, correcting for
                    # wall-forming particles with size 1
                    setattr(self, attr, 2 * getattr(self, attr) - 1.0)
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m_small = number_density_cylinder(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_m_small = volume_fraction_cylinder(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.rho_bulk_m_large = number_density_cylinder(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_m_large = volume_fraction_cylinder(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.rho_bulk_m = self.rho_bulk_m_small + self.rho_bulk_m_large
        self.phi_bulk_m = self.phi_bulk_m_small + self.phi_bulk_m_large

        self.mcrowd = getattr(self, 'dcrowd') ** 3
        self.rho_bulk_c = number_density_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_c = volume_fraction_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )


class TransFociCub(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *TransFociCub* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
      One of multiple chunks of a complete artifact.
    - `whole`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: ns#nl#al#ac#nc#
      Short name for an 'ensemble' artifact.
    - `space`: ns#nl#al#ac#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon_small: float
        Size (diameter) of a monomer
    nmon_small: int
        number of small monomers. Its associated keyword is 'ns'.
    mmon_small: float
        Mass of a small monomer
    dmon_large: float
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int
        number of large monomers. Its associated keyword is 'nl'.
    mmon_large: float, default np.nan
        Mass of a large monomer. Its associated keyword is 'ml'.
    nmon: int
        Total number of monomers.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder.
    lcube : float
        Length of the simulation box, inferred from 'l' keyword
        (half-length of hthe simulation box).
    dt : float
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m_small : float
        Bulk number density fraction of small monomers.
    phi_bulk_m_small : float
        Bulk volume fraction of small monomers
    rho_bulk_m_large : float
        Bulk number density fraction of large monomers.
    phi_bulk_m_large : float
        Bulk volume fraction of large monomers
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    The mass density is uniform across all species. For any species whose mass
    is not explicitly parsed, the mass is defined as :math:`m_{i} = d_{i}^3`,
    where :math:`d_{i}` represents the species' diameter.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = TransFociCub(
    ..."al5.0nl5ml125.0ns400ac1.0nc0l100.0dt0.005bdump2000adump5000.ring"
    ...'ensemble_long',
    ...'all'
    ... )
    >>> print(artifact.nc)
    0
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['bug', 'all']

    _genealogy_attributes = {
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
        'segment': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
        'whole': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc', 'lcube': 'l',
             'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump',
             'ensemble_id': 'ens'}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
        'ensemble_long': OrderedDict(
            {'dmon_large': 'al', 'nmon_large': 'nl', 'mmon_large': 'ml',
             'nmon_small': 'ns', 'dcrowd': 'ac', 'ncrowd': 'nc',
             'lcube': 'l', 'dt': 'dt', 'bdump': 'bdump', 'adump': 'adump'}),
        # Pattern: ns#nl#al#ac#nc#
        'ensemble': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac', 'ncrowd': 'nc'}),
        # Pattern: ns#nl#al#ac#
        'space': OrderedDict(
            {'nmon_small': 'ns', 'nmon_large': 'nl', 'dmon_large': 'al',
             'dcrowd': 'ac'})
    }

    _project_attributes = {
        'segment': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                    'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                    'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'whole': ['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                  'rho_bulk_m_small', 'phi_bulk_m_large', 'phi_bulk_m_large',
                  'rho_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble_long': ['dmon_small', 'mmon_small',  'mcrowd',
                          'phi_bulk_m_small', 'rho_bulk_m_small',
                          'phi_bulk_m_large', 'phi_bulk_m_large', 'rho_bulk_m',
                          'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c'],
        'ensemble': ['dmon_small', 'mmon_small', 'mcrowd'],
        'space': ['dmon_small', 'mmon_small', 'mcrowd']
    }

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['bug', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.mmon_small: float = self.dmon_small**3
        self.mcrowd: float = -1
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m_small: float = -1
            self.rho_bulk_m_small: float = -1
            self.phi_bulk_m_large: float = -1
            self.rho_bulk_m_large: float = -1
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['dmon_large', 'lcube', 'mmon_large', 'dcrowd', 'dt']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'l':
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m_small = number_density_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_small = volume_fraction_cube(
            getattr(self, 'nmon_small'),
            getattr(self, 'dmon_small'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m_large = number_density_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m_large = volume_fraction_cube(
            getattr(self, 'nmon_large'),
            getattr(self, 'dmon_large'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_m = self.rho_bulk_m_small + self.rho_bulk_m_large
        self.phi_bulk_m = self.phi_bulk_m_small + self.phi_bulk_m_large

        self.mcrowd = getattr(self, 'dcrowd') ** 3
        self.rho_bulk_c = number_density_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_c = volume_fraction_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )


class HnsCub(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *HnsCub* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: N#kbmm#nh#ac#l#epshc#nc#ens#.#j#.ring
      One of multiple chunks of a complete artifact.
    - `whole`: N#kbmm#nh#ac#l#epshc#nc#ens#.ring
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: N#kbmm#nh#ac#l#epshc#nc#.ring
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: N#kbmm#nh#ac#epshc#nc#
      Short name for an 'ensemble' artifact.
    - `space`: N#kbmm#nh#ac#epshc#ring
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'nucleoid', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon: float
        Size (diameter) of a monomer
    nmon: int
        number of small monomers. Its associated keyword is 'ns'.
    dhns: float, default 1.0
        Size (diameter) of a H-NS protein.
    dhns_patch: float, default 0.178
        Size (diameter) of a H-NS protein patch at its pole.
    nhns: int
        Number of H-NS protein. Its associated keyword is 'nh'.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    bend_mm: float
        Bending rigidity of DNA monomers. Its associated keyword is 'kbmm'.
    eps_hm: float, default 29.0
        The strength of attractive LJ interaction between hns poles and
        monomers.
    eps_hc: float
        The strength of attractive LJ interaction between H-NS cores and
        crowders. Its associated keyword is 'epshc'.
    lcube : float
        Length of the simulation box, inferred from 'l' keyword
        (half-length of hthe simulation box).
    dt : float, default 0.005
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int, default 5000
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int, default 10000
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_hns : float
        Bulk number density fraction of H-NS proteins.
    phi_bulk_hns : float
        Bulk volume fraction of H-NS proteins.
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = HnsCub(
    ..."N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring"
    ...'whole',
    ...'nucleoid'
    ... )
    >>> print(artifact.nhns)
    8
    """
    _geometry = 'cubic'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']

    _genealogy_attributes = {
        # Pattern: N#kbmm#nh#ac#l#epshc#nc#ens#.#j#.ring
        'segment': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'nhns': 'nh', 'dcrowd': 'ac',
             'lcube': 'l', 'eps_hc': 'epshc', 'ncrowd': 'nc',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: N#kbmm#nh#ac#l#epshc#nc#ens#.ring
        'whole': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'nhns': 'nh', 'dcrowd': 'ac',
             'lcube': 'l', 'eps_hc': 'epshc', 'ncrowd': 'nc',
             'ensemble_id': 'ens'}),
        # Pattern: N#kbmm#nh#ac#l#epshc#nc#.ring
        'ensemble_long': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'nhns': 'nh', 'dcrowd': 'ac',
             'lcube': 'l', 'eps_hc': 'epshc', 'ncrowd': 'nc'}),
        # Pattern: N#kbmm#nh#ac#epshc#nc#
        'ensemble': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'nhns': 'nh', 'dcrowd': 'ac',
             'eps_hc': 'epshc', 'ncrowd': 'nc'}),
        # Pattern: N#kbmm#nh#ac#epshc#ring
        'space': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'nhns': 'nh', 'dcrowd': 'ac',
             'eps_hc': 'epshc'})
    }

    _project_attributes = {
        'segment': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m'
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

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['nucleoid', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.dhns: float = 1
        self.dhns_patch: float = 0.178
        self.dt: float = 0.005
        self.bdump: int = 5000
        self.adump: int = 10000
        self.eps_hm: float = 29
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_hns: float = -1
            self.rho_bulk_hns: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['kbmm', 'lcube', 'dcrowd', 'epshc']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'l':
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m = number_density_cube(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_m = volume_fraction_cube(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_hns = number_density_cube(
            getattr(self, 'nhns'),
            getattr(self, 'dhns'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_hns = volume_fraction_cube(
            getattr(self, 'nhns'),
            getattr(self, 'dhns'),
            getattr(self, 'lcube')
        )
        self.rho_bulk_c = number_density_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )
        self.phi_bulk_c = volume_fraction_cube(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcube')
        )


class HnsCyl(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *HnsCyl* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system attributes:

    - `segment`: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.j#.ring
      One of multiple chunks of a complete artifact.
    - `whole`: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.ring
      A complete artifact. It may be a collection of segments.
    - `ensemble_long`: N#kbmm#r#nh#ac#lz#epshc#nc#.ring
      Detailed name for an 'ensemble' artifact.
    - `ensemble`: N#D#nh#ac#epshc#nc#
      Short name for an 'ensemble' artifact.
    - `space`: N#D#nh#ac#epshc#nc#
      A 'space' artifact.

    For the above four lineages, the short names (eywords) are physical
    attributes where their values (shown by "#" sign) are float or integer
    number. See `genealogy_attributes` below for long names of attribues.

    Other than attributes inhertied from the parent class `ParserBase`, this
    class dynamically defines new attributes based on the list of physical
    attributes of a given `lineage` as define in the `genealogy_attributes`
    class attribute.

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'nucleoid', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    dmon: float
        Size (diameter) of a monomer
    nmon: int
        number of small monomers. Its associated keyword is 'ns'.
    dhns: float, default 1.0
        Size (diameter) of a H-NS protein.
    dhns_patch: float, default 0.178
        Size (diameter) of a H-NS protein patch at its pole.
    nhns: int
        Number of H-NS protein. Its associated keyword is 'nh'.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
    bend_mm: float
        Bending rigidity of DNA monomers. Its associated keyword is 'kbmm'.
    eps_hm: float, default 29.0
        The strength of attractive LJ interaction between hns poles and
        monomers.
    eps_hc: float
        The strength of attractive LJ interaction between H-NS cores and
        crowders. Its associated keyword is 'epshc'.
    lcyl : float
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis).
    dcyl : float
        Size (or diameter) of the cylindrical confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement).
    dt : float, default 0.005
        Simulation timestep. Its associated keyword is 'dt'.
    bdump : int, default 5000
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump : int, default 10000
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id : int
        The ensemble number of a 'whole' artifact in an ensemble. Its
        associated keyword is 'ens'.
    segment_id : int
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a artifact file. Its associated keyword is 'j'.
    rho_bulk_m : float
        Bulk number density fraction of monomers.
    phi_bulk_m : float
        Bulk volume fraction of monomers
    rho_bulk_hns : float
        Bulk number density fraction of H-NS proteins.
    phi_bulk_hns : float
        Bulk volume fraction of H-NS proteins.
    rho_bulk_c : float
        Bulk number density fraction of crowders
    phi_bulk_c : float
        Bulk volume fraction of crowders
    space : str
        A space's name.
    ensemble : str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long : str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole : str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment : str, "N/A"
        A segment's name if applicable, otherwise "N/A"

    Notes
    -----
    The cylindrical wall is implemented in LAMMPS by using wall-forming
    particles of size 1.0. Thus, the actual size of the cylinder size
    (diameter), :math:`D`, is :math:`D=2r-1.0`,  :math:`r` is the radius of
    the cylindrical region defined in LAMMPS.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = HnsCyl(
    ..."N200D20.0nh16ac1.0epshc1.0nc0"
    ...'space',
    ...'nucleoid'
    ... )
    >>> print(artifact.eps_hc)
    1.0
    """
    _geometry = 'cylindrical'
    _topology = 'ring'
    _groups = ['nucleoid', 'all']

    _genealogy_attributes = {
        # Pattern: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.j#.ring
        'segment': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'dcyl': 'r', 'nhns': 'nh',
             'dcrowd': 'ac', 'lcyl': 'lz', 'eps_hc': 'epshc', 'ncrowd': 'nc',
             'ensemble_id': 'ens', 'segment_id': 'j'}),
        # Pattern: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.ring
        'whole': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'dcyl': 'r', 'nhns': 'nh',
             'dcrowd': 'ac', 'lcyl': 'lz', 'eps_hc': 'epshc', 'ncrowd': 'nc',
             'ensemble_id': 'ens'}),
        # Pattern: N#kbmm#r#nh#ac#lz#epshc#nc#.ring
        'ensemble_long': OrderedDict(
            {'nmon': 'N', 'bend_mm': 'kbmm', 'dcyl': 'r', 'nhns': 'nh',
             'dcrowd': 'ac', 'lcyl': 'lz', 'eps_hc': 'epshc', 'ncrowd': 'nc'}),
        # Pattern: N#D#nh#ac#epshc#nc#
        'ensemble': OrderedDict(
            {'nmon': 'N', 'dcyl': 'D', 'nhns': 'nh', 'dcrowd': 'ac',
             'eps_hc': 'epshc', 'ncrowd': 'nc'}),
        # Pattern: N#D#nh#ac#epshc#nc#
        'space': OrderedDict(
            {'nmon': 'N', 'dcyl': 'D', 'nhns': 'nh', 'dcrowd': 'ac',
             'eps_hc': 'epshc'})
    }

    _project_attributes = {
        'segment': ['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m'
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

    def __init__(
        self,
        artifact: str,
        lineage: LineageT,
        group: Literal['nucleoid', 'all']
    ) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon_small: float = 1
        self.dhns: float = 1
        self.dhns_patch: float = 0.178
        self.dt: float = 0.005
        self.bdump: int = 5000
        self.adump: int = 10000
        self.eps_hm: float = 29
        if self._lineage in ['segment', 'whole', 'ensemble_long']:
            self.phi_bulk_m: float = -1
            self.rho_bulk_m: float = -1
            self.phi_bulk_hns: float = -1
            self.rho_bulk_hns: float = -1
            self.phi_bulk_c: float = -1
            self.rho_bulk_c: float = -1

    def _parse_name(self) -> None:
        """
        Parses lineage attributes from the `name` attribute, assigning them
        dynamically as class attributes.

        Notes
        -----
        Lineage attributes are macroscopic physical attributes of the systems.
        They are added to the class dynamically as new class attribute upon
        class instantiation.
        """
        name_strs = re.compile(r"([a-zA-Z\-]+)")
        words = name_strs.split(self._name)
        attrs_float = ['bend_mm', 'dcyl', 'lcyl', 'dcrowd', 'eps_hc']
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == 'lz':
                    # Cylinder full length from its half-length
                    setattr(self, attr, 2 * getattr(self, attr))
                if keyword == 'r':
                    # Cylinder size is twice its radius, correcting for
                    # wall-forming particles with size 1
                    setattr(self, attr, 2 * getattr(self, attr) - 1.0)
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        self.rho_bulk_m = number_density_cylinder(
            getattr(self, 'nmon'),
            getattr(self, 'dmon'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_m = volume_fraction_cylinder(
            getattr(self, 'nhns'),
            getattr(self, 'dhns'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.rho_bulk_hns = number_density_cylinder(
            getattr(self, 'nhns'),
            getattr(self, 'dhns'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_hns = volume_fraction_cylinder(
            getattr(self, 'nhns'),
            getattr(self, 'dhns'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.rho_bulk_c = number_density_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
        self.phi_bulk_c = volume_fraction_cylinder(
            getattr(self, 'ncrowd'),
            getattr(self, 'dcrowd'),
            getattr(self, 'lcyl'),
            getattr(self, 'dcyl')
        )
