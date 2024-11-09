import os
import re
import warnings
from typing import TypeVar, IO, Tuple, Dict, List, Union
from abc import ABC, abstractmethod
from dataclasses import dataclass
from collections import OrderedDict
import numpy as np
import pandas as pd
from .utilizer import invalid_keyword, openany_context, InputT


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

    Class Attributes
    ----------------
    _lineages : list of str
        List of valid lineage types.
    _genealogy : dict of lists
        Dictionary defining the hierarchical relationship between lineages.
        Each lineage points to its parent lineages in a hierarchy.

    Abstract Properties
    -------------------
    _geometry : str
        Specifies geometry of the system
    _topology : str
        Specifies how particles are connected (how) or not in the system.
    _groups : List[str]
        Specifies particle groups in the system
    _genealogy_attributes : Dict[str, OrderedDict[str, str]]
        Dictionary defining lineage-specific attributes for each lineage type.
    _project_attributes : Dict[str, List[str]]
        Dictionary defining project-level attributes that remain constant but
        are not extractable from a filename.
    _attributes : List[str]
        List of all the attributs for an `artifact` with a given `lineage`,
        containing both lineage and project attriubes.

    Methods
    -------
    _find_name() -> None
        Parses and sets the unique name based on `lineage` and `group`.
    _set_parents() -> None
        Sets pattern names for each `lineage` based on `_genealogy`.
    _initiate_attributes() -> None
        Defines and initializes subclass-specific attributes. (Abstract method)
    _parse_lineage_name() -> None
        Parses lineage-specific attributes based on the filename.
        (Abstract method)
    _bulk_attributes() -> None
        Computes physical attributes for the current lineage based on primary
        attributes. (Abstract method)
    """
    _lineages = ["segment", "whole", "ensemble_long", "ensemble", "space"]
    _genealogy = {
        "segment": ["segment", "whole", "ensemble_long", "ensemble", "space"],
        "whole": ["whole", "ensemble_long", "ensemble", "space"],
        "ensemble_long": ["ensemble_long", "ensemble", "space"],
        "ensemble": ["ensemble", "space"],
        "space": ["space"],
    }

    @property
    @abstractmethod
    def _geometry(self) -> str:
        """
        Defines the system geometry for the parser subclass.
        """

    @property
    @abstractmethod
    def _groups(self) -> List[str]:
        """
        List of valid group names for the subclass.
        """

    @property
    @abstractmethod
    def _topology(self) -> str:
        """
        Defines the polymer topology for the parser subclass.
        """

    @property
    @abstractmethod
    def _genealogy_attributes(self) -> Dict[str, OrderedDict[str, str]]:
        """
        Dictionary of lineage-specific attributes. Each key is a lineage type,
        and each value is an OrderedDict mapping attribute names to their
        short-form representations.
        """

    @property
    @abstractmethod
    def _project_attributes(self) -> List[str]:
        """
        Dictionary of project attributes. Each key is a lineage type,
        and each value is an OrderedDict mapping attribute names to their
        short-form representations.
        """

    def __init__(
        self,
        artifact: str,
        lineage: str,
        group: str
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
        self._attributes = (
            list(self._genealogy_attributes[self._lineage].keys())
            + self._project_attributes[self._lineage])
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
    def lineage(self) -> str:
        """
        Gets the current lineage.
        """
        return self._lineage

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
        if self._lineage in ["segment", "whole"]:
            self._name = \
                self.filename.split("." + self._group)[0].split("-")[0]
        else:
            self._name = self._filename.split("-")[0]

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
    def attributes(self) -> Dict[str, OrderedDict[str, str]]:
        """
        Returns lineage-specific attributes for instance `lineage`.
        """
        return self._attributes[self._lineage]

    @property
    def genealogy_attributes(self) -> Dict[str, OrderedDict[str, str]]:
        """
        Returns lineage-specific attributes for the all the lineages.
        """
        return self._genealogy_attributes

    @property
    def project_attributes(self) -> Dict[str, List[str]]:
        """
        Returns project-level attributes.
        """
        return self._project_attributes

    @abstractmethod
    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.
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


class TwoMonDep(ParserBase):
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

    Other than attributes inhertied from the parent class `ParserBase` as
    explained below, this class dynamically defines new attributes based on the
    list of physical attributes of a given `lineage` as define in the
    `genealogy_attributes` class attribute.

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
    attributes: List[str]
        List of attributes specific to an `artifact` with a given `lineage`,
        including parsed and computed attributes.
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
    rho_m_bulk : float
        Bulk number density fraction of monomers.
    phi_m_bulk : float
        Bulk volume fraction of monomers
    rho_c_bulk : float
        Bulk number density fraction of crowders
    phi_c_bulk : float
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

    Class Attributes
    ----------------
    _geometry : str
        Specifies geometry of the system
    _topology : str
        Specifies how particles are connected (how) or not in the system.
    _groups : List[str]
        Specifies particle groups in the system
    _genealogy_attributes : Dict[str, Dict[str, str]]
        Maps `lineage` names to attribute keywords for parsing.
    _project_attributes : Dict[str, List[Optional[str]]]
        Specifies additional physical attributes for each `lineage`.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = (
    ..."nm2am5.0ac1.0"
    ..."space",
    ..."bug"
    ... )
    Artifact('nm2am5.0ac1.0' in geometry 'cubic' from group 'bug' with lineage 'space' and topology 'atomic')
    """
    _geometry = "cubic"
    _topology = "atomic"
    _groups = ["bug", "all"]

    _genealogy_attributes: Dict[str, OrderedDict[str, str]] = {
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#.j#
        "segment": OrderedDict({
            "nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump", "ensemble_id": "ens",
            "segment_id": "j"}
            ),
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#
        "whole": OrderedDict({
            "nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump", "ensemble_id": "ens"}
            ),
        # Pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump# :
        "ensemble_long": OrderedDict({
            "dmon": "am", "nmon": "nm", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump"}
            ),
        # Pattern: nm#am#ac#nc#sd :
        "ensemble": OrderedDict(
            {"nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
             "d_sur": "sd"}
             ),
        # pttern: nm#am#ac# :
        "space": OrderedDict(
            {"nmon": "nm", "dmon": "am",  "dcrowd": "ac"}
            )
    }

    _project_attributes: Dict[str, List[str]] = {
        "segment": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "whole": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "ensemble_long": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk",
                          "rho_c_bulk"],
        "ensemble": [],
        "space": []
        }

    def __init__(self, artifact: str, lineage: str, group: str) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        self._lineage_genealogy = super()._genealogy[lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        if self._lineage in ["segment", "whole", "ensemble_long"]:
            self.phi_m_bulk: float = -1
            self.rho_m_bulk: float = -1
            self.phi_c_bulk: float = -1
            self.rho_c_bulk: float = -1

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
        attrs_float = ["dmon", "lcube", "dcrowd", "dt", "d_sur"]
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == "hl":
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        vol_cell = getattr(self, 'lcube') ** 3
        vol_mon = np.pi * getattr(self, 'dmon') ** 3 / 6
        self.rho_m_bulk = getattr(self, 'nmon') / vol_cell
        self.phi_m_bulk = self.rho_m_bulk * vol_mon
        vol_crowd = np.pi * getattr(self, 'dcrowd') ** 3 / 6
        self.rho_c_bulk = getattr(self, 'ncrowd') / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


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

    Other than attributes inhertied from the parent class `ParserBase` as
    explained below, this class dynamically defines new attributes based on the
    list of physical attributes of a given `lineage` as define in the
    `genealogy_attributes` class attribute.

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
    attributes: List[str]
        List of attributes specific to an `artifact` with a given `lineage`,
        including parsed and computed attributes.
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
    rho_m_bulk : float
        Bulk number density fraction of monomers.
    phi_m_bulk : float
        Bulk volume fraction of monomers
    rho_c_bulk : float
        Bulk number density fraction of crowders
    phi_c_bulk : float
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

    Class Attributes
    ----------------
    _geometry : str
        Specifies geometry of the system
    _topology : str
        Specifies how particles are connected (how) or not in the system.
    _groups : List[str]
        Specifies particle groups in the system
    _genealogy_attributes : Dict[str, Dict[str, str]]
        Maps `lineage` names to attribute keywords for parsing.
    _project_attributes : Dict[str, List[Optional[str]]]
        Specifies additional physical attributes for each `lineage`.

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
    ..."ensemble",
    ..."all"
    ... )
    ... print(artifact.dcrowd)
    2.0
    """
    _geometry = "cylindrical"
    _topology = "linear"
    _groups = ["bug", "all"]

    _genealogy_attributes: Dict[str, OrderedDict[str, str]] = {
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
        "segment":  OrderedDict(
            {"nmon": "N", "epsilon": "epsilon", "dcyl": "r", "lcyl": "lz",
             "dcrowd": "sig", "ncrowd": "nc", "dt": "dt", "bdump": "bdump",
             "adump": "adump", "ensemble_id": "ens", "segment_id": "j"}),
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
        "whole":  OrderedDict(
            {"nmon": "N", "epsilon": "epsilon", "dcyl": "r", "lcyl": "lz",
             "dcrowd": "sig", "ncrowd": "nc", "dt": "dt", "bdump": "bdump",
             "adump": "adump", "ensemble_id": "ens"}),
        # Pattern: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
        "ensemble_long":  OrderedDict(
            {"nmon": "N", "epsilon": "epsilon", "dcyl": "r", "lcyl": "lz",
             "dcrowd": "sig", "ncrowd": "nc", "dt": "dt", "bdump": "bdump",
             "adump": "adump"}),
        # Pattern: N#D#ac#nc#
        "ensemble":  OrderedDict(
            {"nmon": "N", "dcyl": "D", "dcrowd": "ac", "ncrowd": "nc"}),
        # Pattern: N#D#ac#
        "space":  OrderedDict({"nmon": "N", "dcyl": "D", "dcrowd": "ac"})
    }
    _project_attributes: Dict[str, List[str]] = {
        "segment": ["dmon", "phi_m_bulk", "rho_m_bulk", "phi_c_bulk",
                    "rho_c_bulk"],
        "whole": ["dmon",  "phi_m_bulk", "rho_m_bulk", "phi_c_bulk",
                  "rho_c_bulk"],
        "ensemble_long": ["dmon", "phi_m_bulk", "rho_m_bulk", "phi_c_bulk",
                          "rho_c_bulk"],
        "ensemble": ["dmon"],
        "space": ["dmon"]
    }

    def __init__(self, artifact: str, lineage: str, group: str) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        self._lineage_genealogy = super()._genealogy[self.lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._dependant_attributes()

    def _initiate_attributes(self) -> None:
        """
        Defines and initiates the project attributes.

        Notes
        -----
        The negative initial values are unphysical.
        """
        self.dmon: float = 1
        if self._lineage in ["segment", "whole", "ensemble_long"]:
            self.phi_m_bulk: float = -1
            self.rho_m_bulk: float = -1
            self.phi_c_bulk: float = -1
            self.rho_c_bulk: float = -1

    def _parse_lineage_name(self) -> None:
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
        attrs_float = ["dmon", "dcyl", "lcyl", "epsilon", "dcrowd", "dt"]
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == "lz":
                    # Cylinder full length from its half-length
                    setattr(self, attr, 2 * getattr(self, attr))
                if keyword == "r":
                    # Cylinder size is twice its radius, correcting for
                    # wall-forming particles with size 1
                    setattr(self, attr, 2 * getattr(self, attr) - 1.0)
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _dependant_attributes(self) -> None:
        """
        Calculates system attributes based on parsed values.
        """
        vol_cell_m = (np.pi
                      * (getattr(self, 'dcyl') - getattr(self, 'dmon'))**2
                      * getattr(self, 'lcyl') / 4.0)
        vol_cell_c = (np.pi
                      * (getattr(self, 'dcyl') - getattr(self, 'dcrowd'))**2
                      * getattr(self, 'lcyl') / 4.0)
        vol_mon = np.pi * getattr(self, 'dmon') ** 3 / 6
        self.rho_m_bulk = getattr(self, 'nmon') / vol_cell_m
        self.phi_m_bulk = self.rho_m_bulk * vol_mon
        vol_crowd = np.pi * getattr(self, 'dcrowd') ** 3 / 6
        self.rho_c_bulk = getattr(self, 'ncrowd') / vol_cell_c
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class SumRuleCubHeteroRing(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *SumRuleCyl* project, utilizing specific filename patterns.

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

    Other than attributes inhertied from the parent class `ParserBase` as
    explained below, this class dynamically defines new attributes based on the
    list of physical attributes of a given `lineage` as define in the
    `genealogy_attributes` class attribute.

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
    attributes: List[str]
        List of attributes specific to an `artifact` with a given `lineage`,
        including parsed and computed attributes.
    dmon_small: float
        Size (diameter) of a monomer
    nmon_small: int
        number of small monomers. Its associated keyword is 'ns'.
    dmon_large: float
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int
        number of large monomers. Its associated keyword is 'nl'.
    nmon: int
        Total number of monomers.
    dcrowd: float
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd : int
        Number of crowders. Its associated keyword is 'nc'.
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
    rho_m_bulk : float
        Bulk number density fraction of monomers.
    phi_m_bulk : float
        Bulk volume fraction of monomers
    rho_c_bulk : float
        Bulk number density fraction of crowders
    phi_c_bulk : float, default np.na
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

    Class Attributes
    ----------------
    _geometry : str
        Fixed to 'cylindrical' for TwoMonDep project.
    _topology : str
        Specifies polymer topology; e.g., 'atomic' or 'linear'.
    _groups : List[str]
        Possible particle groups for the `TwoMonDep` project.
    _genealogy_attributes : Dict[str, Dict[str, str]]
        Maps `lineage` names to attribute keywords for parsing.
    _project_attributes : Dict[str, List[Optional[str]]]
        Specifies additional physical attributes for each `lineage`.

    Examples
    --------
    Creating a instance to parse a filename with specified lineage and group.

    >>> artifact = SumRuleCubHeteroRing(
    ..."ns400nl5am5.0ac1.0"
    ..."space",
    ..."bug"
    ... )
    ... print(artifact.dcrowd)
    1.0
    """
    _geometry = "cubic"
    _topology = "ring"
    _groups = ["bug", "all"]

    _genealogy_attributes: Dict[str, OrderedDict[str, str]] = {
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
        "segment": OrderedDict(
            {"dmon_large": "al", "nmon_large": "nl", "mmon_large": "ml",
             "nmon_small": "ns", "dcrowd": "ac", "ncrowd": "nc",
             "lcube": "l", "dt": "dt", "bdump": "bdump", "adump": "adump",
             "ensemble_id": "ens", "segment_id": "j"}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
        "whole": OrderedDict(
            {"dmon_large": "al", "nmon_large": "nl", "mmon_large": "ml",
             "nmon_small": "ns", "dcrowd": "ac", "ncrowd": "nc", "lcube": "l",
             "dt": "dt", "bdump": "bdump", "adump": "adump",
             "ensemble_id": "ens"}),
        # Pattern: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
        "ensemble_long": OrderedDict(
            {"dmon_large": "al", "nmon_large": "nl", "mmon_large": "ml",
             "nmon_small": "ns", "dcrowd": "ac", "ncrowd": "nc",
             "lcube": "l", "dt": "dt", "bdump": "bdump", "adump": "adump"}),
        # Pattern: ns#nl#al#ac#nc#
        "ensemble": OrderedDict(
            {"nmon_small": "ns", "nmon_large": "nl", "dmon_large": "al",
             "dcrowd": "ac", "ncrowd": "nc"}),
        # Pattern: ns#nl#al#ac#
        "space": OrderedDict(
            {"nmon_large": "nl", "nmon_small": "ns", "dmon_large": "al",
             "dcrowd": "ac"})
    }

    _project_attributes: Dict[str, List[str]] = {
        "segment": ["dmon_small", "mmon_small", "mcrowd", "phi_m_bulk",
                    "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "whole": ["dmon_small", "mmon_small", "mcrowd", "phi_m_bulk",
                  "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "ensemble_long": ["dmon_small", "mmon_small", "mcrowd", "phi_m_bulk",
                          "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "ensemble": ["dmon_small", "mmon_small", "mcrowd"],
        "space": ["dmon_small", "mmon_small", "mcrowd"]
    }

    def __init__(self, artifact: str, lineage: str, group: str,) -> None:
        super().__init__(artifact, lineage, group)
        self._initiate_attributes()
        self._parse_name()
        self._set_parents()
        self._lineage_genealogy: List[str] = super()._genealogy[self.lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

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
        if self._lineage in ["segment", "whole", "ensemble_long"]:
            self.phi_m_bulk: float = -1
            self.rho_m_bulk: float = -1
            self.phi_c_bulk: float = -1
            self.rho_c_bulk: float = -1

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
        attrs_float = ["dmon_large", "lcube", "mmon_large", "dcrowd", "dt"]
        for attr, keyword in self._genealogy_attributes[self._lineage].items():
            try:
                val = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(val) if attr in attrs_float else int(float(val)))
                if keyword == "l":
                    # Cube full side from its half-side
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _bulk_attributes(self) -> None:
        """
        computes some physical attributes of a lineage based on its
        primary attributes.
        """
        self.mcrowd = getattr(self, "dcrowd") ** 3
        vol_cell = getattr(self, 'lcube') ** 3
        vol_mon_s = np.pi * getattr(self, 'dmon_small') ** 3 / 6
        vol_mon_l = np.pi * getattr(self, 'dmon_large') ** 3 / 6
        self.rho_m_bulk = getattr(self, 'nmon') / vol_cell
        self.phi_m_bulk = (
            vol_mon_s * getattr(self, 'nmon_small')
            + vol_mon_l * getattr(self, 'nmon_large')
        ) / vol_cell
        vol_crowd = np.pi * getattr(self, 'dcrowd') ** 3 / 6
        self.rho_c_bulk = getattr(self, 'ncrowd') / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class SumRuleCubHeteroLinear(ParserBase):
    """
    parses a `name` (which can be a filename or filepath based on the value
    of `ispath` argument) to extract information about a project's file
    based on a pattern pre-defined by the `lineage` in the
    'cubic' geometry for a given `group` in the project.

    In the geometry 'cubic', these patterns are used to parse a `name`:

        segment: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.linear
            One of multiple chunks of a complete simulation or measurement.
        whole: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.linear
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.linear
            Long name of an ensemble.
        ensemble: ns#nl#al#ac#nc#
            A collection of 'wholes' (complete simulations) that differs only
            in their initial conditions (e.g., random number seed).
        space: ns#nl#al#ac#
            A collection of ensembles with a unique set of all the input
            parameters except number of crowders (nc).

    In the above lineages, the keywords are attributes where their values
    (shown by "#" sign) are float or integer number. If a lineage does not
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:

        whole: 'j'
        ensemble_long: 'ens', and 'j'
        ensemble: 'ens', 'j', 'l', 'dt', 'bdump', 'adump', and 'ml'
        space: 'ens' , 'j', 'l', 'dt', 'bdump', 'adump', 'ml', and 'nc'

    There are some difference between the keywords of physical attributes and
    their associated attributes in the `TransFuci` class. Below, the these two
    types of attributes are explained.

    The mass density is uniform and is the same for all the species so,
    m = d ** 3 is used to define mass for species whose masses are not parsed.

    Parameters
    ----------
    name: str
        Name that is parsed for extracting information.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}
        Type of the lineage of the name.
    geometry: 'cylindrical'
        Shape of the simulation box.
    group: {'bug', 'all'}
        Type of the particle group. 'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    topology:
        The topology of the polymer.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    _pathname: str, default "N/A"
        Equal to `name` if `name` is a filepath, otherwise "N/A".
    _filename: str
        Name of a the file referred to by `name` if `name` is a filepath,
        otherwise the `name` itself.
    _lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}, default whole
        Type of the lineage of the name.
    _geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    _group: {'bug', 'all'}
        Type of the particle group.  'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    _ispath: bool, default True
        Whether the name is a filepath or a simple name.
    _topology: str, default 'linear'
        The topology of the polymer.
    lineage_name: str,
        The unique name of type extracted from self.fullname
    dmon_small: float, default 1
        Size (diameter) of a monomer
    nmon_small: int, np.nan
        number of small monomers. Its associated keyword is 'ns'.
    mmon_small: float, default dmon_small**3
        Mass of a small monomer
    dmon_large: float, np.nan
        Size (diameter) of a large monomer. Its associated keyword is 'al'.
    nmon_large: int, np.nan
        number of large monomers. Its associated keyword is 'nl'.
    mmon_large: float, default np.nan
        Mass of a large monomer. Its associated keyword is 'ml'.
    nmon: int, np.nan
        Total number of monomers.
    acrowd: int, np.nan
        Size of crowders. Its associated keyword is 'ac'.
    ncrowd: int, np.nan
        number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder.
    eps_others: float, default 1
        Unit of the LJ interaction strength.
    lcube: float, np.nan
        Side of the cubic simulation box.
    dt: float, np.nan
        Simulation timestep. Its associated keyword is 'dt'.
    bdump: int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump: int, default np.nan
        Frequency by which 'all' configurations are dumped in a 'segment'
        trajectory file. Its associated keyword is 'adump'.
    ensemble_id: int, np.nan
        The ensemble number of a 'whole' simulation in an ensemble. Its
        associated keyword is 'ens'.
    segment_id: int, np.nan
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a whole file. Its associated keyword is 'j'.
    space: str
        A space's name
    ensemble: str, "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long: str, "N/A"
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole: str, "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment: str. "N/A"
        A segment's name if applicable, otherwise "N/A"
    rho_m_bulk: float, default np.nan
        Bulk number density fraction of monomers
    phi_m_bulk: float, default np.nan
        Bulk volume fraction of monomers
    rho_c_bulk: float, default np.nan
        Bulk number density fraction of crowders
    phi_c_bulk: float, default np.nan
        Bulk volume fraction of crowders

    Class Attributes
    ----------------
    _groups: list of str
        Possible groups of the `SumRule` project.
    _lineage_attributes: dict of dict
        a dictionary of `lineage` names. For each `lineage`, a dictionary
        maps the keywords of physical attributes in that lineage to their
        corresponding attributes in this class.
    _project_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created. This
        attributes are either used in the simulation/experiment/run but not use
        the `name`, or are created within this class.
    """
    _groups = ["bug", "all"]
    _lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # lcube twice of l
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # lcube twice of l
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # lcube twice of l
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
        },
        "ensemble": {
            "dmon_large": "al",
            "nmon_large": "nl",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
        },
        "space": {
            "dmon_large": "al",
            "nmon_large": "nl",
            "nmon_small": "ns",
            "dcrowd": "ac",
        },
    }
    _project_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon_small",
            "mmon_small",
            "mcrowd",
            "eps_others",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "whole": [
            "dmon_small",
            "mmon_small",
            "mcrowd",
            "eps_others",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble_long": [
            "dmon_small",
            "mmon_small",
            "mcrowd",
            "eps_others",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble": ["dmon_small", "mmon_small", "mcrowd", "eps_others"],
        "space": ["dmon_small", "mmon_small", "mcrowd", "eps_others"],
    }
    _geometry_error = "'TransFociCub' is used for the 'cubic' geometry."

    def __init__(
        self,
        name: str,
        lineage: str,
        geometry: str,
        group: str,
        topology: str,
        ispath: bool = True,
    ) -> None:
        invalid_keyword(geometry, ["cubic"], self._geometry_error)
        invalid_keyword(group, self._groups)
        super().__init__(name, lineage, geometry, group, topology, ispath)
        self._initiate_attributes()
        self._parse_lineage_name()
        self._set_parents()
        self.attributes = (
            list(self._lineage_attributes[self.lineage].keys())
            + self._project_attributes[self.lineage]
        )
        self.genealogy: List[str] = super()._genealogy[self.lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

    def _initiate_attributes(self) -> None:
        """
        defines and initiates the class attributes based on the physical
        attributes defined for the project.

        The negative initial values are unphysical.
        """
        # group attributes
        self.dmon_small: float = 1
        self.nmon_small: int = -1
        self.mmon_small: float = self.dmon_small**3
        self.nmon_large: int = -1
        self.dmon_large: float = -1
        self.mmon_large: float = -1
        self.nmon: int = -1
        self.phi_m_bulk: float = -1
        self.rho_m_bulk: float = -1
        self.dcrowd: float = -1
        self.ncrowd: int = -1
        self.mcrowd: float = -1
        self.phi_c_bulk: float = -1
        self.rho_c_bulk: float = -1
        # system attributes
        self.ensemble_id: int = -1
        self.segment_id: int = -1
        self.dt: float = 0
        self.bdump: int = -1
        self.adump: int = -1
        # cubic attributes
        self.eps_others: float = 1
        self.lcube: float = -1

    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.lineage_name)
        attributes_float = ["dmon_large",
                            "lcube", "mmon_large", "dcrowd", "dt"]
        for attr_n, attr_kw in self._lineage_attributes[self.lineage].items():
            try:
                attr_value = words[words.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                if attr_kw == "l":
                    attr_value = 2 * attr_value
                setattr(self, attr_n, attr_value)
            except ValueError:
                print(
                    f"'{attr_kw}'"
                    " attribute keyword is not in "
                    f"'{self.lineage_name}'"
                    " lineage name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not."
                )
        setattr(self, "nmon", self.nmon_large + self.nmon_small)
        self.mcrowd = self.dcrowd**3
        warnings.warn("Mass is scaled with sized", UserWarning)

    def _set_parents(self) -> None:
        """
        set to parent names for a lineage_name based on its lineage.

        The following map is used for setting relationships:

            'segment': A child of 'whole' lineage.
            'whole': A child of 'ensemble' lineage.
            'ensemble': A child of 'space' lineage.
            'space': The root of other lineages.

        It is assumed that 'nc' is the last attribute short-key in a
        lineage_name of types: 'ensemble', 'ensemble_long', 'whole', 'segment'.
        """
        convention_warning = (
            "It is assumed that 'nc' is the last attribute"
            + " short-key in a lineage_name of types:"
            + " 'ensemble', 'ensemble_long', 'whole', 'segment'."
        )
        space_name = (
            "ns"
            + str(self.nmon_small)
            + "nl"
            + str(self.nmon_large)
            + "al"
            + str(self.dmon_large)
            + "ac"
            + str(self.dcrowd)
        )
        ensemble_name = space_name + "nc" + str(self.ncrowd)
        if self.lineage == "space":
            self.space = self.lineage_name
            self.ensemble = "N/A"
            self.ensemble_long = "N/A"
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == "ensemble":
            self.space = self.lineage_name.split("nc")[0]
            warnings.warn(convention_warning, UserWarning)
            self.ensemble = self.lineage_name
            self.ensemble_long = "N/A"
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == "ensemble_long":
            self.space = space_name
            self.ensemble = ensemble_name
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == "whole":
            self.space = space_name
            self.ensemble = ensemble_name
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name.split("ens")[0]
            self.whole = self.lineage_name
            self.segment = "N/A"
        else:
            self.space = space_name
            self.ensemble = ensemble_name
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name.split("ens")[0]
            self.whole = self.lineage_name.split(".j")[0]
            self.segment = self.lineage_name

    def _bulk_attributes(self) -> None:
        """
        computes some physical attributes of a lineage based on its
        primary attributes.
        """
        vol_cell = self.lcube**3
        vol_mon_s = np.pi * self.dmon_small**3 / 6
        vol_mon_l = np.pi * self.dmon_large**3 / 6
        self.rho_m_bulk = self.nmon / vol_cell
        self.phi_m_bulk = (
            vol_mon_s * self.nmon_small + vol_mon_l * self.nmon_large
        ) / vol_cell
        vol_crowd = np.pi * self.dcrowd**3 / 6
        self.rho_c_bulk = self.ncrowd / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd