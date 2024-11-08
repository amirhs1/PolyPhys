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
        String defining the system geometry, such as 'cylindrical' or 'cubic'.
    _groups : List[Union[str, None]]
        List of valid group names for the subclass.
    _topology : str
        String defining the polymer topology, such as 'linear'.
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

    @property
    @abstractmethod
    def _attributes(self) -> List[str]:
        """
        List of all the attributs for an `artifact` with a given `lineage`,
        containing both lineage and project attriubes.
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
        self._lineage = lineage
        self._group = group
        self._project_name = self.__class__.__name__
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
        Parses the lineage-specific attributes based on the filename.
        """

    @abstractmethod
    def _bulk_attributes(self) -> None:
        """
        Computes physical attributes for the current lineage based on primary
        attributes.
        """


class TwoMonDep(ParserBase):
    """
    Extracts structured information about an artifact from its name in the
    *TwoMonDep* project, utilizing specific filename patterns.

    Each lineage level has a unique naming pattern used to parse key physical
    and system parameters:

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

    Other than the attribute inhertied from the parent class `ParserBase` as
    explained below

    Parameters
    ----------
    artifact : str
        Name to be parsed, either a filename or filepath.
    lineage : {'segment', 'whole', 'ensemble_long', 'ensemble', 'space'}
        Type of the lineage of the name.
    group : {'bug', 'nucleoid', 'all'}
        Particle group type, with `bug` representing a single polymer.

    Attributes
    ----------
    attributes: List[str]
        List of attributes specific to *TwoMonDep* project,
        including parsed and computed attributes. Details include:
        List of class attribute set dynamically within the class. Details
        include:

        - dmon : float
            Size (diameter) of a monomer. Its associated keyword is 'am'.
        - nmon : int
            Number of monomers. Its associated keyword is 'N'.
        - dcrowd: float
            Size (diameter) of a crowder. Its associated keyword is 'ac'.
        - ncrowd : int
            Number of crowders. Its associated keyword is 'nc'.
        - lcube : float
            Length of the simulation box, inferred from 'hl' keyword
            (half-length of hthe simulation box).
        - d_sur : float
            Surface-to-surface distance between two monomers fixed in space.
            Its associated keyword is 'sd'.
        - dt : float
            Simulation timestep. Its associated keyword is 'dt'.
        - bdump : int
            Frequency by which 'bug' configurations are dumped in a 'bug'
            trajectory file. Its associated keyword is 'bdump'.
        - adump : int
            Frequency by which 'all' configurations are dumped in a 'segment'
            trajectory file. Its associated keyword is 'adump'.
        - tdump : int
            Frequency by which 'themo' variables are written in a 'lammps'
            log file. Its associated keyword is 'tdump'.
        - ensemble_id : int
            The ensemble number of a 'whole' artifact in an ensemble. Its
            associated keyword is 'ens'.
        - segment_id : int, np.nan
            The 'segment_id' keyword starts with 'j', ends with a 'padded'
            number such as '05' or '14', showing the succession of segments
            in a artifact file. Its associated keyword is 'j'.
        - space : str
            A space's name.
        - ensemble : str, "N/A"
            An ensemble's name if applicable, otherwise "N/A"
        - ensemble_long : str, "N/A"
            The name of ensemble derived from 'whole' name if applicable,
            otherwise "N/A"
        - whole : str, "N/A"
            A whole's name if applicable, otherwise "N/A"
        - segmen : str, "N/A"
            A segment's name if applicable, otherwise "N/A"
        - rho_m_bulk : float, default np.nan
            Bulk number density fraction of monomers
        - phi_m_bulk : float, default np.nan
            Bulk volume fraction of monomers
        - rho_c_bulk : float, default np.nan
            Bulk number density fraction of crowders
        - phi_c_bulk : float, default np.nan
            Bulk volume fraction of crowders

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
    """
    _geometry = "cubic"
    _topology = "atomic"
    _groups = ["bug", "all"]

    _genealogy_attributes: Dict[str, Dict[str, str]] = {
        # pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#.j#
        "segment": OrderedDict({
            "nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump", "ensemble_id": "ens",
            "segment_id": "j"}
            ),
        # pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump#ens#
        "whole": OrderedDict({
            "nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump", "ensemble_id": "ens"}
            ),
        # pattern: am#nm#ac#nc#hl#sd#dt#bdump#adump$tdump# :
        "ensemble_long": OrderedDict({
            "dmon": "am", "nmon": "nm", "dcrowd": "ac", "ncrowd": "nc",
            "lcube": "hl", "d_sur": "sd", "dt": "dt", "bdump": "bdump",
            "adump": "adump", "tdump": "tdump"}
            ),
        # pattern: nm#am#ac#nc#sd :
        "ensemble": OrderedDict(
            {"nmon": "nm", "dmon": "am", "dcrowd": "ac", "ncrowd": "nc",
             "d_sur": "sd"}
             ),
        # pttern: nm#am#ac# :
        "space": OrderedDict(
            {"nmon": "nm", "dmon": "am",  "dcrowd": "ac"}
            )
    }

    _project_attributes = {
        "segment": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "whole": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk"],
        "ensemble_long": ["phi_m_bulk", "rho_m_bulk", "phi_c_bulk",
                          "rho_c_bulk"],
        "ensemble": [],
        "space": []
        }

    def __init__(self, artifact: str, lineage: str, group: str) -> None:
        invalid_keyword(group, self._groups)
        super().__init__(artifact, lineage, group)
        self._parse_name()
        self._set_parents()
        self._attributes = (
            list(self._genealogy_attributes[self.lineage].keys())
            + self._project_attributes[self.lineage]
        )
        self._lineage_genealogy = super()._genealogy[lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

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
                value = words[words.index(keyword) + 1]
                setattr(self,
                        attr,
                        float(value) if attr in attrs_float else int(float(value)))
                if attr == "lcube":
                    setattr(self, attr, 2 * getattr(self, attr))
            except ValueError:
                print(f"'{keyword}' attribute not found in '{self._name}'")

    def _bulk_attributes(self) -> None:
        """
        Calculates system attributes and bulk properties based on parsed
        values.
        """
        vol_cell = getattr(self, 'lcube') ** 3
        vol_mon = np.pi * getattr(self, 'dmon') ** 3 / 6
        self.rho_m_bulk = getattr(self, 'nmon') / vol_cell
        self.phi_m_bulk = self.rho_m_bulk * vol_mon
        vol_crowd = np.pi * getattr(self, 'dcrowd') ** 3 / 6
        self.rho_c_bulk = getattr(self, 'ncrowd') / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd
