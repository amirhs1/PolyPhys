from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
import pandas as pd
import os
import re
from typing import TypeVar, IO, Tuple, Dict, List
import warnings
from collections import OrderedDict
from .utilizer import invalid_keyword, openany_context, InputT

TExcludedVolume = TypeVar("TExcludedVolume", bound="ExcludedVolume")
TFreeEnergyVirial = TypeVar("TFreeEnergyVirial", bound="FreeEnergyVirial")


class ParserBase(ABC):
    """
    parses a `name` (which can be a filename or  a filepath based on the
    value of `ispath` argument) to extract information about a
    project's file based on a pattern pre-defined by the `lineage` in a given
    `geometry` for a given `group`.

    Parameters
    ----------
    name: str
        Name that is parsed for extracting information.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}
        Type of the lineage of the name.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.  'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    topology:
        The topology of the polymer.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    filepath: str, default "N/A"
        Equal to `name` if `name` is a filepath, otherwise "N/A".
    filename: str
        Name of a file referred to by `name` if `name` is a filepath,
        otherwise the `name` itself.
    _lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}, default whole
        Type of the lineage of the name.
    _geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    _group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.  'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    _ispath: bool, default True
        Whether a `name` is a filename or a filepath.
    _topology: str, default 'linear'
        The topology of the polymer.
    lineage_name: str
        The unique name of type extracted from self.fullname

    Class Attributes
    ----------------
    _lineage_attributes: dict of None
        A dictionary of `lineage` names. For each `lineage`, a dictionary
        maps the keywords of physical attributes in that `lineage` of a project
        to their corresponding attributes in the project's associated parser
        class.
    _physical_attributes: dict of None
        A dictionary of `lineage` names. For each `lineage`, a list of parser
        attributes that either are defined or have default values for the
        `lineage` in the class.
    _genealogy: dict of lists
        A dictionary of `lineage` names. For each `lineage`, a list is defined
        that contains the parent-like lineage attributes of that `lineage`.
    """
    _lineage_attributes: Dict[str, None] = {
        "segment": None,
        "whole": None,
        "ensemble_long": None,
        "ensemble": None,
        "space": None,
    }
    _physical_attributes: Dict[str, None] = {
        "segment": None,
        "whole": None,
        "ensemble_long": None,
        "ensemble": None,
        "space": None,
    }
    _genealogy: Dict[str, List[str]] = {
        "segment": [
            "lineage_name",
            "segment",
            "whole",
            "ensemble_long",
            "ensemble",
            "space",
        ],
        "whole": [
            "lineage_name",
            "whole",
            "ensemble_long",
            "ensemble",
            "space"
        ],
        "ensemble_long": [
            "lineage_name",
            "ensemble_long",
            "ensemble",
            "space",
        ],
        "ensemble": ["lineage_name", "ensemble", "space"],
        "space": ["lineage_name", "space"],
    }

    def __init__(
        self,
        name: str,
        lineage: str,
        geometry: str,
        group: str,
        topology: str,
        ispath: bool = True,
    ) -> None:
        self.filepath = name
        self.filename = ''
        invalid_keyword(lineage, list(self._lineage_attributes.keys()))
        self._lineage = lineage
        self._geometry = geometry
        self._group = group
        self._topology = topology
        self._ispath = ispath
        if self._ispath is True:
            self.filename, _ = os.path.splitext(self.filepath)
            self.filename: str = self.filename.split("/")[-1]
        else:
            self.filename = name
            self.filepath = "N/A"
        self._find_lineage_name()

    def __str__(self) -> str:
        observation: str = f"""
        Observation:
            Name: '{self.filename}',
            Geometry: '{self._geometry},
            Group: '{self._group}',
            Lineage: '{self._lineage}',
            Polymer topology: '{self._topology}'
        """
        return observation

    def __repr__(self) -> str:
        return (
            f"Observation('{self.filename}' in geometry"
            + f" '{self._geometry}' from group '{self._group}' with"
            + f" lineage '{self._lineage}' and "
            + f" polymer topology '{self._topology}')"
        )

    @property
    def lineage(self) -> str:
        """
        gets the current `lineage`.

        Returns
        -------
        str
            returns the current `lineage`.
        """
        return self._lineage

    def _find_lineage_name(self) -> None:
        """
        parses the unique lineage_name (the first substring of filename
        and/or the segment keyword middle substring) of a filename.
        """
        if self._lineage in ["segment", "whole"]:
            # a 'segment' lineage only used in 'probe' phase
            # a 'whole' lineage used in 'probe' or 'analyze' phases
            # so its lineage_name is either ended by 'group' keyword or "-".
            # these two combined below:
            self.lineage_name: str = \
                self.filename.split("." + self._group)[0].split("-")[0]
        else:  # 'ensemble' or 'space' lineages
            self.lineage_name = self.filename.split("-")[0]

    @property
    def lineage_attributes(self) -> Dict[str, None]:
        """
        gets the `lineage_attributes`.

        Returns
        -------
        str
            returns the `lineage_attributes`.
        """
        return self._lineage_attributes

    @property
    def physical_attributes(self) -> Dict[str, None]:
        """
        gets the `physical_attributes`.

        Returns
        -------
        str
            returns the `physical_attributes`.
        """
        return self._physical_attributes

    @property
    def ispath(self) -> bool:
        """
        gets the current `ispath`.

        Returns
        -------
        str
            returns the current `ispath`.
        """
        return self._ispath

    @property
    def geometry(self) -> str:
        """
        gets the current `geometry`.

        Returns
        -------
        str
            returns the current `geometry`.
        """
        return self._geometry

    @property
    def group(self) -> str:
        """
        gets the current `group`.

        Returns
        -------
        str
            returns the current `group`.
        """
        return self._group

    @property
    def topology(self) -> str:
        """
        gets the current `ispath`.

        Returns
        -------
        str
            returns the current `ispath`.
        """
        return self._topology

    @abstractmethod
    def _initiate_attributes(self) -> None:
        """
        defines and initiates the class attributes based on the physical
        attributes defined for the project.
        """
        pass

    @abstractmethod
    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        pass

    @abstractmethod
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
        pass

    @abstractmethod
    def _bulk_attributes(self) -> None:
        """
        computes some physical attributes of a lineage based on its
        primary attributes.
        """
        pass


class SumRuleCyl(ParserBase):
    name: str
    lineage: str
    geometry: str
    group: str
    topology: str
    ispath: bool = True
    """
    parses a `name` (which can be a filename or a file path based on the
    value of `ispath` argument) to extract information about a project's
    file based on a pattern pre-defined by the `lineage` in the
    'cylindrical' geometry for a given `group` in the project.

    In the geometry 'cylindrical', these patterns are used to parse a `name`:

        segment: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
            One of multiple chunks of a complete simulation or measurement.
        whole: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
            Long name of an ensemble.
        ensemble: N#D#ac#nc#
            A collection of 'wholes' (complete simulations) that differs only
            in their initial conditions (e.g., random number seed).
        space: N#D#ac#
            A collection of ensembles with a unique set of all the input
            parameters except number of crowders (nc).

    In the above lineages, the keywords are attributes where their values
    (shown by "#" sign) are float or integer number. If a lineage does not
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:

        whole: 'j'
        ensemble_long: 'ens', and 'j'
        ensemble: 'lz', 'dt', 'bdump', 'adump', 'ens', and 'j'
        space:  'lz', 'nc', 'dt', 'bdump', 'adump', 'ens' , and 'j'

    There are some difference between the keywords of physical attributes and
    their associated attributes in the `SumRule` class. Below, the these two
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
    group: {'bug', 'nucleoid', 'all'}
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
    _group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.  'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    _ispath: bool, default True
        Whether the name is a filepath or a simple name.
    _topology: str, default 'linear'
        The topology of the polymer.
    lineage_name: str,
        The unique name of the type extracted from self.fullname
    dmon: float, default 1
        Size (diameter) of a monomer
    nmon: int, np.nan
        Number of monomers. Its associated keyword is 'N'.
    mmon: float, default dmon**3
        Mass of a monomer
    dcrowd: float, np.nan
        Size (diameter) of a crowder. Its associated keyword is 'sig' or 'ac'.
    ncrowd: int, np.nan
        Number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder
    dwall: float, default 1
        Wall-forming particles diameter
    epsilon: float, np.nan
        Wall-particle LJ interaction strength. Its associated keyword is
        'epsilon' keyword.
    eps_others: float, default 1
        Other LJ interaction strengths
    dcyl: float, np.nan
        Size (or diameter) of the `cylindrical` or `slit` confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement. Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl/2] inclusive, is the
        domain over which x and y cartesian coordinates are defined in the
        'cylindrical' geometry; however, `dcyl`, ranged from '[-dcyl/2,dcyl.2]
        inclusive, is the domain over which z cartesian coordinate is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the cylinder size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
    lcyl: float, np.nan
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis.
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
    segment: str, "N/A"
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
    _physical_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created. This
        attributes are either used in the simulation/experiment/run but not use
        the `name`, or are created within this class.
    """
    _groups = ["bug", "all"]
    _lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # dcyl twice of r
            "nmon": "N",
            "epsilon": "epsilon",
            "dcyl": "r",
            "lcyl": "lz",
            "dcrowd": "sig",
            "ncrowd": "nc",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # dcyl twice of r
            "nmon": "N",
            "epsilon": "epsilon",
            "dcyl": "r",
            "lcyl": "lz",
            "dcrowd": "sig",
            "ncrowd": "nc",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # dcyl twice of r
            "nmon": "N",
            "epsilon": "epsilon",
            "dcyl": "r",
            "lcyl": "lz",
            "dcrowd": "sig",
            "ncrowd": "nc",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
        },
        "ensemble": {"nmon": "N", "dcyl": "D", "dcrowd": "ac", "ncrowd": "nc"},
        "space": {"nmon": "N", "dcyl": "D", "dcrowd": "ac"},
    }
    _physical_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon",
            "mmon",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "whole": [
            "dmon",
            "mmon",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble_long": [
            "dmon",
            "mmon",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble": ["dmon", "mmon", "eps_others", "mcrowd", "dwall"],
        "space": ["dmon", "mmon", "eps_others", "mcrowd", "dwall"],
    }
    _geometry_error = "'SumRuleCyl' is used for the 'cylindrical' geometry."

    def __init__(
        self,
        name: str,
        lineage: str,
        geometry: str,
        group: str,
        topology: str,
        ispath: bool = True,
    ) -> None:
        invalid_keyword(geometry, ["cylindrical"], self._geometry_error)
        invalid_keyword(group, self._groups)
        super().__init__(name, lineage, geometry, group, topology, ispath)
        self._initiate_attributes()
        self._parse_lineage_name()
        self._set_parents()
        self.attributes = (
            list(self._lineage_attributes[self.lineage].keys())
            + self._physical_attributes[self.lineage]
        )
        self.genealogy: List[str] = super()._genealogy[lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

    def _initiate_attributes(self) -> None:
        """
        defines and initiates the class attributes based on the physical
        attributes defined for the project.

        The negative initial values are unphysical.
        """
        # group attributes
        self.dmon: float = 1
        self.nmon: int = -1
        self.mmon: float = self.dmon**3
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
        # geometry attributes:
        self.dcyl: float = -1
        self.lcyl: float = -1
        self.dwall: float = 1
        self.eps_others: float = 1
        self.epsilon: float = 0

    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.lineage_name)
        attributes_float = ["dmon", "dcyl", "lcyl", "epsilon", "dcrowd", "dt"]
        for attr_n, attr_kw in self._lineage_attributes[self.lineage].items():
            try:
                attr_value = words[words.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                if attr_kw == "lz":
                    attr_value = 2 * attr_value
                if attr_kw == "r":
                    attr_value = 2 * attr_value - self.dwall
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
        self.mcrowd = self.dcrowd**3

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
            self.space = (
                "N" + str(self.nmon) + "D" + str(self.dcyl) +
                "ac" + str(self.dcrowd)
            )
            self.ensemble = self.space + "nc" + str(self.ncrowd)
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == "whole":
            self.space = (
                "N" + str(self.nmon) + "D" + str(self.dcyl) +
                "ac" + str(self.dcrowd)
            )
            self.ensemble = self.space + "nc" + str(self.ncrowd)
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name.split("ens")[0]
            self.whole = self.lineage_name
            self.segment = "N/A"
        else:
            self.space = (
                "N" + str(self.nmon) + "D" + str(self.dcyl) +
                "ac" + str(self.dcrowd)
            )
            self.ensemble = self.space + "nc" + str(self.ncrowd)
            warnings.warn(convention_warning, UserWarning)
            self.ensemble_long = self.lineage_name.split("ens")[0]
            self.whole = self.lineage_name.split(".j")[0]
            self.segment = self.lineage_name

    def _bulk_attributes(self) -> None:
        """
        computes some physical attributes of a lineage based on its
        primary attributes.
        """
        vol_cell = np.pi * self.dcyl**2 * self.lcyl / 4.0
        vol_mon = np.pi * self.dmon**3 / 6
        self.rho_m_bulk = self.nmon / vol_cell
        self.phi_m_bulk = self.rho_m_bulk * vol_mon
        vol_crowd = np.pi * self.dcrowd**3 / 6
        self.rho_c_bulk = self.ncrowd / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class TransFociCyl(ParserBase):
    name: str
    lineage: str
    geometry: str
    group: str
    topology: str
    ispath: bool = True
    """
    parses a `name` (which can be a filename or filepath based on the value
    of `ispath` argument) to extract information about a project's file
    based on a pattern pre-defined by the `lineage` in the
    'cylindrical' geometry for a given `group` in the project.

    In the geometry 'cylindrical', these patterns are used to parse a `name`:

        segment: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.j#.ring
            One of multiple chunks of a complete simulation or measurement.
        whole: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#ens#.ring
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: epss#epsl#r#al#nl#ml#ns#ac#nc#lz#dt#bdump#adump#.ring
            Long name of an ensemble.
        ensemble: ns#nl#al#D#ac#nc#
            A collection of 'wholes' (complete simulations) that differs only
            in their initial conditions (e.g., random number seed).
        space: ns#nl#al#D#ac#
            A collection of ensembles with a unique set of all the input
            parameters except number of crowders (nc).

    In the above lineages, the keywords are attributes where their values
    (shown by "#" sign) are float or integer number. If a lineage does not
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:

        whole: 'j'
        ensemble_long: 'ens', and 'j'
        ensemble: 'ens', 'j', 'lz', 'dt', 'bdump', 'adump', and 'ml'
        space: 'ens' , 'j', 'lz', 'dt', 'bdump', 'adump', 'ml', and 'nc'

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
    dcrowd: float, np.nan
        Size (diameter) of a crowder. Its associated keyword is 'ac'.
    ncrowd: int, np.nan
        number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder
    dwall: float, default 1
        Wall-forming particles diameter
    epsilon_s: float, np.nan
        Wall-small-monomer LJ interaction strength. Its associated keyword is
        'epss' keyword.
    epsilon_l: float, np.nan
        Wall-large-monomer LJ interaction strength. Its associated keyword is
        'espl' keyword.
    eps_others: float, default 1
        Other LJ interaction strengths
    dcyl: float, np.nan
        Size (or diameter) of the `cylindrical` or `slit` confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement. Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl/2] inclusive, is the
        domain over which x and y cartesian coordinates are defined in the
        'cylindrical' geometry; however, `dcyl`, ranged from '[-dcyl/2,dcyl.2]
        inclusive, is the domain over which z cartesian coordinate is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the cylinder size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
    lcyl: float, np.nan
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis.
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
    segment: str, "N/A"
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
    _physical_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created. This
        attributes are either used in the simulation/experiment/run but not use
        the `name`, or are created within this class.
    """
    _groups = ["bug", "all"]
    _lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # dcyl twice of r
            "epsilon_small": "epss",
            "epsilon_large": "epss",
            "dcyl": "r",
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # dcyl twice of r
            "epsilon_small": "epss",
            "epsilon_large": "epss",
            "dcyl": "r",
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # dcyl twice of r
            "epsilon_small": "epss",
            "epsilon_large": "epss",
            "dcyl": "r",
            "dmon_large": "al",
            "nmon_large": "nl",
            "mmon_large": "ml",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
            "dt": "dt",
            "bdump": "bdump",
            "adump": "adump",
        },
        "ensemble": {
            "dcyl": "D",
            "dmon_large": "al",
            "nmon_large": "nl",
            "nmon_small": "ns",
            "dcrowd": "ac",
            "ncrowd": "nc",
        },
        "space": {
            "dcyl": "D",
            "dmon_large": "al",
            "nmon_large": "nl",
            "nmon_small": "ns",
            "dcrowd": "ac",
        },
    }
    _physical_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon_small",
            "mmon_small",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "whole": [
            "dmon_small",
            "mmon_small",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble_long": [
            "dmon_small",
            "mmon_small",
            "eps_others",
            "mcrowd",
            "dwall",
            "phi_m_bulk",
            "rho_m_bulk",
            "phi_c_bulk",
            "rho_c_bulk",
        ],
        "ensemble": [
            "dmon_small", "mmon_small", "eps_others", "mcrowd", "dwall"],
        "space": ["dmon_small", "mmon_small", "eps_others", "mcrowd", "dwall"],
    }
    _geometry_error = "'TransFociCyl' is used for the 'cylindrical' geometry."

    def __init__(
        self,
        name: str,
        lineage: str,
        geometry: str,
        group: str,
        topology: str,
        ispath: bool = True,
    ) -> None:
        invalid_keyword(geometry, ["cylindrical"], self._geometry_error)
        invalid_keyword(group, self._groups)
        super().__init__(name, lineage, geometry, group, topology, ispath)
        self._initiate_attributes()
        self._parse_lineage_name()
        self._set_parents()
        self.attributes = (
            list(self._lineage_attributes[self.lineage].keys())
            + self._physical_attributes[self.lineage]
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
        # geometry attributes
        self.lcyl: float = -1
        self.dwall: float = 1
        self.eps_others: float = 1
        self.dcyl: float = -1
        self.epsilon_small: float = -1
        self.epsilon_large: float = -1

    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.lineage_name)
        attributes_float = ["dmon_large", "dcyl", "lcyl", "epsilon_small"
                            "epsilon_large", "mmon_large", "dcrowd", "dt"]
        for attr_n, attr_kw in self._lineage_attributes[self.lineage].items():
            try:
                attr_value = words[words.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                if attr_kw == "lz":
                    attr_value = 2 * attr_value
                if attr_kw == "r":
                    attr_value = 2 * attr_value - self.dwall
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
            + "D"
            + str(self.dcyl)
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
        vol_cell = np.pi * self.dcyl**2 * self.lcyl / 4.0
        vol_mon_s = np.pi * self.dmon_small**3 / 6
        vol_mon_l = np.pi * self.dmon_large**3 / 6
        self.rho_m_bulk = self.nmon / vol_cell
        self.phi_m_bulk = (
            vol_mon_s * self.nmon_small + vol_mon_l * self.nmon_large
        ) / vol_cell
        vol_crowd = np.pi * self.dcrowd**3 / 6
        self.rho_c_bulk = self.ncrowd / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class TransFociCub(ParserBase):
    name: str
    lineage: str
    geometry: str
    group: str
    topology: str
    ispath: bool = True
    """
    parses a `name` (which can be a filename or filepath based on the value
    of `ispath` argument) to extract information about a project's file
    based on a pattern pre-defined by the `lineage` in the
    'cubic' geometry for a given `group` in the project.

    In the geometry 'cubic', these patterns are used to parse a `name`:

        segment: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.j#.ring
            One of multiple chunks of a complete simulation or measurement.
        whole: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#ens#.ring
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: al#nl#ml#ns#ac#nc#l#dt#bdump#adump#.ring
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
    _physical_attributes: dict of lists
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
    _physical_attributes: Dict[str, List[str]] = {
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
            + self._physical_attributes[self.lineage]
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


class HnsCub(ParserBase):
    name: str
    lineage: str
    geometry: str
    group: str
    topology: str
    ispath: bool = True
    """
    parses a `name` (which can be a filename or a pathname based on the value
    of `ispath` argument) to extract information about a project's file
    based on a pattern pre-defined by the `lineage` in the
    'cubic' geometry for a given `group` in the project.

    In the geometry 'cubic', these patterns are used to parse a `name`:
        segment: N#kbmm#nh#ac#l#epshc#nc#ens.#j#.ring
            One of multiple chunks of a complete simulation or measurement.
        whole: N#kbmm#nh#ac#l#epshc#nc#ens.ring
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: N#kbmm#nh#ac#l#epshc#nc#.ring
            Long name of an ensemble.
        ensemble: N#kbmm#nh#ac#epshc#nc#
            A collection of 'wholes' (complete simulations) that differs only
            in their initial conditions (e.g., random number seed).
        space: N#kbmm#nh#ac#epshc#
            A collection of ensembles with a unique set of all the input
            parameters except number of crowders (nc).

    In the above lineages, the keywords are attributes where their values
    (shown by "#" sign) are float or integer number. If a lineage does not
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:

        whole: 'j'
        ensemble_long: 'ens', and 'j'
        ensemble: 'ens', 'j', 'l', 'dt', 'ndump', and 'adump'
        space: 'ens' , 'j', 'l', 'dt', 'ndump', 'adump', and 'nc'

    There are some difference between the keywords of physical attributes and
    their associated attributes in the `HnsCub` class. Below, the these two
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
    geometry: 'cubic'
        Shape of the simulation box.
    group: {'nucleoid', 'all'}
        Type of the particle group. 'nucleoid' is used for a single polymer
        with its nucleoid-associated proteins within the nucleoid. 'all' is
        used for all the particles/atoms in the system.
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
    _group: {'nucleoid', 'all'}
        Type of the particle group.  'nucleoid' is used for a single polymer
        with its nucleoid-associated proteins within the nucleoid. 'all' is
        used for all the particles/atoms in the system.
    _ispath: bool, default True
        Whether the name is a filepath or a simple name.
    _topology: str, default 'linear'
        The topology of the polymer.
    lineage_name: str,
        The unique name of type extracted from self.fullname
    dmon: float, default 1
        Size (diameter) of a monomer
    nmon: int, np.nan
        number of monomers. Its associated keyword is 'N'.
    mmon: float, default dmon**3
        Mass of a small monomer
    dhns: float, default 1
        Size (diameter) of a hns protein
    dhns_patch: float, default 0.178
        Size (diameter) of a hns protein patch at its pole.
    nhns: int, np.nan
        number of hns protein. Its associated keyword is 'nh'.
    mhns: float, default dhns**3
        Mass of a hns protein
    dcrowd: float, np.nan
        Size (diameter) of a monomer. Its associated keyword is 'ac'.
    ncrowd: int, np.nan
        number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder
    bend_mm: float, default np.nan
        Bending rigidity of DNA monomers. Its associated keyword is 'kbmm'.
    eps_hm: float, default np.nan
        The strength of attractive LJ interaction between hns poles and
        monomers. Its associated keyword is 'epshm'.
    eps_hc: float, default np.nan
        The strength of attractive LJ interaction between hns cores and
        crowders. Its associated keyword is 'epshc'.
    eps_others: float, default 1
        Unit of the LJ interaction strength
    lcube: float, np.nan
        Side of the cubic simulation box.
    dt: float, np.nan
        Simulation timestep. Its associated keyword is 'dt'.
    ndump: int
        Frequency by which 'bug' configurations are dumped in a 'nucleoid'
        trajectory file. Its associated keyword is 'ndump'.
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
    rho_hns_bulk: float, default np.nan
        Bulk number density fraction of hns proteins
    phi_hns_bulk: float, default np.nan
        Bulk volume fraction of hns proteins
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
    _physical_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created. This
        attributes are either used in the simulation/experiment/run but not use
        the `name`, or are created within this class.
    """
    _groups = ["nucleoid", "all"]
    _lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
        },
        "ensemble": {
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hm": "epshm",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
        },
        "space": {
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
        },
    }
    _physical_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "dt", "ndump", "adump"
            ],
        "whole": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk",  "dt", "ndump", "adump"
            ],
        "ensemble_long": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "dt", "ndump", "adump"
            ],
        "ensemble": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "dt", "ndump", "adump"
            ],
        "space": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "dt", "ndump", "adump"
            ]
    }
    _geometry_error = "'HnsCub' is used for the 'cubic' geometry."

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
            + self._physical_attributes[self.lineage]
        )
        self.genealogy: List[str] = super()._genealogy[self.lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

    def _initiate_attributes(self) -> None:
        """
        defines and initiates the class attributes based on the physical
        attributes defined for the project.

        All negative initial values are unphysical.
        """
        # group attributes
        self.dmon: float = 1
        self.nmon: int = -1
        self.mmon: float = self.dmon**3
        self.phi_m_bulk: float = -1
        self.rho_m_bulk: float = -1
        self.dhns: float = 1
        self.dhns_patch: float = 0.178
        self.nhns: int = -1
        self.mhns: float = self.dhns**3
        self.phi_hns_bulk: float = -1
        self.rho_hns_bulk: float = -1
        self.dcrowd: float = -1
        self.ncrowd: int = -1
        self.mcrowd: float = -1
        self.phi_c_bulk: float = -1
        self.rho_c_bulk: float = -1
        # system attributes
        self.ensemble_id: int = -1
        self.segment_id: int = -1
        self.dt: float = 0.005
        self.ndump: int = 5000
        self.adump: int = 10000
        # cubic attributes
        self.eps_hm: float = 29
        self.bend_mm: float = -1
        self.eps_hc: float = -1
        self.eps_others: float = 1
        self.lcube: float = -1

    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.lineage_name)
        attributes_float = ["kbmm", "lcube", "dcrowd", "epshc"]
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
        self.mcrowd = self.dcrowd**3

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
            "N"
            + str(self.nmon)
            + "kbmm"
            + str(self.bend_mm)
            + "nh"
            + str(self.nhns)
            + "ac"
            + str(self.dcrowd)
            + "epshc"
            + str(self.eps_hc)
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
        vol_mon = np.pi * self.dmon**3 / 6
        self.rho_m_bulk = self.nmon / vol_cell
        self.phi_m_bulk = vol_mon * self.nmon / vol_cell
        vol_hns = np.pi * self.dhns**3 / 6
        self.rho_hns_bulk = self.nhns / vol_cell
        self.phi_hns_bulk = vol_hns * self.nhns / vol_cell
        vol_crowd = np.pi * self.dcrowd**3 / 6
        self.rho_c_bulk = self.ncrowd / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class HnsCyl(ParserBase):
    name: str
    lineage: str
    geometry: str
    group: str
    topology: str
    ispath: bool = True
    """
    parses a `name` (which can be a filename or a pathname based on the value
    of `ispath` argument) to extract information about a project's file
    based on a pattern pre-defined by the `lineage` in the
    'cubic' geometry for a given `group` in the project.

    In the geometry 'cylindrical', these patterns are used to parse a `name`:
        segment: N#kbmm#r#nh#ac#lz#epshc#nc#ens#.j#.ring
            One of multiple chunks of a complete simulation or measurement.
        whole: N#kbmm#r#nh#ac#r#lz#epshc#nc#ens#.ring
            A complete simulation or measurement; a collection of 'segments'.
        ensemble_long: N#kbmm#r#nh#ac#lz#epshc#nc#.ring
            Long name of an ensemble.
        ensemble: N#D#nh#ac#epshc#nc#
            A collection of 'wholes' (complete simulations) that differs only
            in their initial conditions (e.g., random number seed).
        space: N#D#nh#ac#epshc#
            A collection of ensembles with a unique set of all the input
            parameters except number of crowders (nc).

    In the above lineages, the keywords are attributes where their values
    (shown by "#" sign) are float or integer number. If a lineage does not
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:

        whole: 'j'
        ensemble_long: 'ens', and 'j'
        ensemble: 'ens', 'j', 'l',
        space: 'ens' , 'j', 'l', and 'nc'

    There are some difference between the keywords of physical attributes and
    their associated attributes in the `HnsCub` class. Below, the these two
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
    group: {'nucleoid', 'all'}
        Type of the particle group. 'nucleoid' is used for a single polymer
        with its nucleoid-associated proteins within the nucleoid. 'all' is
        used for all the particles/atoms in the system.
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
    _group: {'nucleoid', 'all'}
        Type of the particle group.  'nucleoid' is used for a single polymer
        with its nucleoid-associated proteins within the nucleoid. 'all' is
        used for all the particles/atoms in the system.
    _ispath: bool, default True
        Whether the name is a filepath or a simple name.
    _topology: str, default 'linear'
        The topology of the polymer.
    lineage_name: str,
        The unique name of type extracted from self.fullname
    dmon: float, default 1
        Size (diameter) of a monomer
    nmon: int, np.nan
        number of monomers. Its associated keyword is 'N'.
    mmon: float, default dmon**3
        Mass of a small monomer
    dhns: float, default 1
        Size (diameter) of a hns protein
    dhns_patch: float, default 0.178
        Size (diameter) of a hns protein patch at its pole.
    nhns: int, np.nan
        number of hns protein. Its associated keyword is 'nh'.
    mhns: float, default dhns**3
        Mass of a hns protein
    dcrowd: float, np.nan
        Size (diameter) of a monomer. Its associated keyword is 'ac'.
    ncrowd: int, np.nan
        number of crowders. Its associated keyword is 'nc'.
    mcrowd: float, default np.nan
        Mass of a crowder
    dwall: float, default 1
        Wall-forming particles diameter
    eps_hm: float, default 29
        The strength of attractive LJ interaction between hns poles and
        monomers. Its associated keyword is 'epshm'.
    eps_hc: float, default np.nan
        The strength of attractive LJ interaction between hns cores and
        crowders. Its associated keyword is 'epshc'.
    bend_mm: float, default np.nan
        Bending rigidity of DNA monomers.
    eps_others: float, default 1
        Unit of the LJ interaction strength
    dcyl: float, np.nan
        Size (or diameter) of the `cylindrical` or `slit` confinement, inferred
        from either 'r' keyword (the radius of a cylindrical confinement
        with open ends) or 'D' keyword (size of that confinement. Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl/2] inclusive, is the
        domain over which x and y cartesian coordinates are defined in the
        'cylindrical' geometry; however, `dcyl`, ranged from '[-dcyl/2,dcyl.2]
        inclusive, is the domain over which z cartesian coordinate is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the cylinder size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
    lcyl: float, np.nan
        Length of the cylindrical confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        cylindrical confinement along z axis.
    dt: float, default 0.005
        Simulation timestep. Its associated keyword is 'dt'.
    ndump: int, default 5000
        Frequency by which 'bug' configurations are dumped in a 'nucleoid'
        trajectory file. Its associated keyword is 'ndump'.
    adump: int, default 10000
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
    rho_hns_bulk: float, default np.nan
        Bulk number density fraction of hns proteins
    phi_hns_bulk: float, default np.nan
        Bulk volume fraction of hns proteins
    rho_c_bulk: float, default np.nan
        Bulk number density fraction of crowders
    phi_c_bulk: float, default np.nan
        Bulk volume fraction of crowders

    Class Attributes
    ----------------
    _groups: list of str
        Possible groups of the `HnsCyl` project.
    _lineage_attributes: dict of dict
        a dictionary of `lineage` names. For each `lineage`, a dictionary
        maps the keywords of physical attributes in that lineage to their
        corresponding attributes in this class.
    _physical_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created. This
        attributes are either used in the simulation/experiment/run but not use
        the `name`, or are created within this class.
    """
    _groups = ["nucleoid", "all"]
    _lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            #"eps_hc": "epshc",
            "nhns": "nh",
            "dcyl": "r",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            #"eps_hc": "epshc",
            "nhns": "nh",
            "dcyl": "r",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            #"eps_hc": "epshc",
            "nhns": "nh",
            "dcyl": "r",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcyl": "lz",
        },
        "ensemble": {
            "nmon": "N",
            "nhns": "nh",
            #"eps_hc": "epshc",
            "dcyl": "D",
            "dcrowd": "ac",
            "ncrowd": "nc",
        },
        "space": {
            "nmon": "N",
            "nhns": "nh",
            #"eps_hc": "epshc",
            "dcyl": "D",
            "dcrowd": "ac"
        },
    }
    _physical_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "eps_hm", "dt", "ndump", "adump"
            ],
        "whole": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "eps_hm", "dt", "ndump", "adump"
            ],
        "ensemble_long": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "eps_hm", "dt", "ndump", "adump"
            ],
        "ensemble": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others", "eps_hm",
            "dt", "ndump", "adump"
            ],
        "space": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others", "eps_hm",
            "dt", "ndump", "adump"
            ]
    }
    _geometry_error = "'HnsCyl' is used for the 'cylindrical' geometry."

    def __init__(
        self,
        name: str,
        lineage: str,
        geometry: str,
        group: str,
        topology: str,
        ispath: bool = True,
    ) -> None:
        invalid_keyword(geometry, ["cylindrical"], self._geometry_error)
        invalid_keyword(group, self._groups)
        super().__init__(name, lineage, geometry, group, topology, ispath)
        self._initiate_attributes()
        self._parse_lineage_name()
        self._set_parents()
        self.attributes = (
            list(self._lineage_attributes[self.lineage].keys())
            + self._physical_attributes[self.lineage]
        )
        self.genealogy: List[str] = super()._genealogy[self.lineage]
        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

    def _initiate_attributes(self) -> None:
        """
        defines and initiates the class attributes based on the physical
        attributes defined for the project.

        All negative initial values are unphysical.
        """
        # group attributes
        self.dmon: float = 1
        self.nmon: int = -1
        self.mmon: float = self.dmon**3
        self.phi_m_bulk: float = -1
        self.rho_m_bulk: float = -1
        self.dhns: float = 1
        self.dhns_patch: float = 0.178
        self.nhns: int = -1
        self.mhns: float = self.dhns**3
        self.phi_hns_bulk: float = -1
        self.rho_hns_bulk: float = -1
        self.dcrowd: float = -1
        self.ncrowd: int = -1
        self.mcrowd: float = -1
        self.phi_c_bulk: float = -1
        self.rho_c_bulk: float = -1
        # system attributes
        self.ensemble_id: int = -1
        self.segment_id: int = -1
        self.dt: float = 0.005
        self.ndump: int = 5000
        self.adump: int = 10000
        # cylindrical attributes
        self.bend_mm: float = -1
        self.eps_hm: float = 29
        self.eps_hc: float = -1
        self.eps_others: float = 1
        self.dcyl: float = -1
        self.lcyl: float = -1
        self.dwall: float = 1

    def _parse_lineage_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.lineage_name)
        attributes_float = ["bend_mm", "dcyl", "lcyl", "dcrowd"]
        for attr_n, attr_kw in self._lineage_attributes[self.lineage].items():
            try:
                attr_value = words[words.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                if attr_kw == "lz":
                    attr_value = 2 * attr_value
                if attr_kw == "r":
                    attr_value = 2 * attr_value - self.dwall
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
        self.mcrowd = self.dcrowd**3

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
            "N"
            + str(self.nmon)
            + "D"
            + str(self.dcyl)
            + "nh"
            + str(self.nhns)
            + "ac"
            + str(self.dcrowd)
            + "epshc"
            + str(self.eps_hc)
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
        vol_cell = np.pi * self.dcyl**2 * self.lcyl / 4
        vol_mon = np.pi * self.dmon**3 / 6
        self.rho_m_bulk = self.nmon / vol_cell
        self.phi_m_bulk = vol_mon * self.nmon / vol_cell
        vol_hns = np.pi * self.dhns**3 / 6
        self.rho_hns_bulk = self.nhns / vol_cell
        self.phi_hns_bulk = vol_hns * self.nhns / vol_cell
        vol_crowd = np.pi * self.dcrowd**3 / 6
        self.rho_c_bulk = self.ncrowd / vol_cell
        self.phi_c_bulk = self.rho_c_bulk * vol_crowd


class FloryChain(object):
    """
    extract the attributes of a Flory-type polymer chain from the `filename`
    of its size data and then clean the data itself.

    The chain size data are obtained solving a Flory equation for a linear
    chain with or without the three-body term. The excluded volume and
    three-body interaction coefficient, each, are model in two different
    ways.

    Parameters
    ----------
    name: str
        Name that is parsed for extracting information.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    dcyl: float, np.nan
        Size (or diameter) of the `cylindrical` or `slit` confinement, inferred
        from either 'r' keyword (radius of the cylindrical confinement with
        open ends) or 'D' keyword (size of that confinement). Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl.2] inclusive is the
        domain over which x and y cartesian coordinates are defined in the
        'cylindrical' geometry; however, `dcyl` ranged from '[-dcyl/2,dcyl.2]
        inclusive is the domain over which z cartesian coordinate is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
        Hence, `dcyl` is defined differently in different parser classes in
        this module.

    Class Attributes
    ----------------

    """

    _physical_attributes = {
        "deGennes": {
            "dimension": "dim",
            "w3body": "w",
            "nmon": "n",
            "dcyl": "d",
            "dcrowd": "ac",
            },
        "twoBodyOnly": {
            "dimension": "dim",
            "w3body": "w",
            "nmon": "n",
            "dcyl": "d",
            "dcrowd": "ac",
            },
        "shendruk": {
            "dimension": "dim",
            "nmon": "n",
            "dcyl": "d",
            "dcrowd": "ac"
            },
    }
    _column_names = {
        "deGennes": ["phi_c", "vexc", "swollen"],
        "twoBodyOnly": ["phi_c", "vexc", "swollen"],
        "shendruk": ["phi_c", "vexc", "w3body", "swollen"],
    }

    def __init__(self, name: str, limit: bool = False, ispath: bool = True
                 ) -> None:
        self.limit = limit
        self.ispath = ispath
        if self.ispath:
            self.filepath = name
            self.filename, _ = os.path.splitext(self.filepath)
            self.filename = self.filename.split("/")[-1]
        else:
            self.filepath = "N/A"
            self.filename = name
        self._initiate_attributes()
        self._parse_attributes()
        self._athermal_virial_coeffs()
        self._initial_chain_size()
        self._clean()
        self._expand()
        self._scale()

    def __str__(self) -> str:
        observation = f"""
        A Flory chain:
            Name: '{self.filename}',
        """
        return observation

    def __repr__(self) -> str:
        return f"A Flory chain('{self.filename}')"

    def _initiate_attributes(self):
        """
        defines and initiates the class attribute based on the physical
        attributes defined for the project.
        """
        self.w3body_model = "N/A"
        self.vexc_model = "N/A"
        self.dmon = 1
        self.nmon = np.nan
        self.dcyl = np.nan
        self.dcrowd = np.nan
        self.dimension = np.nan
        self.w3body = np.nan

    def _parse_attributes(self):
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        words = self.filename.split("-")
        self.w3body_model = words[1]
        self.vexc_model = words[2]
        attrs_pattern = re.compile(r"([a-zA-Z\-]+)")
        attrs = attrs_pattern.split(words[-1])
        attributes_float = ["dmon", "dcyl", "dcrowd", "three_body"]
        for attr_n, attr_kw in self._physical_attributes[
            self.w3body_model
        ].items():
            try:
                attr_value = attrs[attrs.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                setattr(self, attr_n, attr_value)
            except ValueError:
                print(
                    f"'{attr_kw}'"
                    " attribute keyword is not in "
                    f"'{self.filename}'"
                    " lineage name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not."
                )

    def _athermal_virial_coeffs(self):
        """
        compute the excluded volume (second virial coefficient) and
        three-body coefficient (third virial coefficient) for the
        equivalent athermal system.
        """
        if self.vexc_model == "AOExcVol":
            # The hard-sphere potential's virial coefficient
            self.vexc_athr = (4 / 3) * np.pi * self.dmon**3
            self.w3body_athr = (5 / 8) * self.vexc_athr**2
        elif self.vexc_model == "HaExcVol":
            # The fully-repulsive Weeks–Chandler–Anderson (WCA) potential
            # with r_cut = 2**(1/6) * sigma where sigma=self.dmon=1
            self.vexc_athr = 4.40944631
            self.w3body_athr = (5 / 8) * self.vexc_athr**2
        else:
            vexc_models = ["AOExcVol", "HaExcVol"]
            raise ValueError(
                f"'{self.vexc_model}' "
                "is not a valid model of the excluded volume of monomers."
                f"Please select one of these '{vexc_models}' models."
            )

    def _initial_chain_size(self):
        """
        compute the equilibrium chain size in the given dimension in
        the absence of crowders.
        """
        self.flory_exponent = 3 / (self.dimension + 2)
        if self.dimension == 3:
            self.r_flory = 1.12 * self.dmon * self.nmon**self.flory_exponent
        elif self.dimension == 2:
            print("In 2-dimensional space, the pre-factor is set to 1.")
            self.r_flory = (
                1
                * self.dmon
                * (
                    self.nmon**self.flory_exponent
                    * (self.dmon / self.dcyl) ** (1 / 4)
                )
            )
        elif self.dimension == 1:
            self.r_flory = (
                1.37
                * self.dmon
                * (
                    self.nmon**self.flory_exponent
                    * (self.dmon / self.dcyl) ** (2 / 3)
                )
            )
        else:
            raise ValueError(
                f"'{self.dimension}' "
                "is not a valid dimension for the system."
                f"Please select one of these '{list(range(1,4,1))}' models."
            )
        self.r_ideal = self.dmon * self.nmon**0.5

    def _clean(self) -> None:
        """
        clean a dataset of a chain's size based on its attributes.
        """
        if self.ispath:
            df = pd.read_csv(
                self.filepath,
                index_col=False,
                names=self._column_names[self.w3body_model],
            )
        else:
            raise ValueError(f"'{self.filepath}' " "is not a valid filepath.")
        df.dropna(inplace=True)
        df.sort_values(by=["phi_c"], inplace=True)
        df.reset_index(inplace=True, drop=True)
        if self.w3body_model == "deGennes":
            if self.nmon == 80:
                df = df.groupby("phi_c").min()
                df.reset_index(inplace=True)
        elif self.w3body_model == "shendruk":
            df = df.groupby("phi_c").min()
            df.reset_index(inplace=True)
        else:
            df = df.sort_values(by=["swollen"])
        self.chain_size = df

    def _expand(self):
        """
        expand `self.chain_size` by adding the physical attributes to it as
         new columns.
        """
        self.chain_size["w3body_model"] = self.w3body_model
        self.chain_size["vexc_model"] = self.vexc_model
        for attr_n in self._physical_attributes[self.w3body_model].keys():
            self.chain_size[attr_n] = getattr(self, attr_n)

    def _scale(self, limit=False):
        """
        scale a chain's attributes.
        """
        self.chain_size["phi_c_scaled"] = (
            self.dmon * self.chain_size["phi_c"] / self.dcrowd
        )
        # phi_c*a/a_c for a_c << a for maximum free energy gain:
        self.chain_size["vexc_scaled"] = self.chain_size["vexc"] / \
            self.vexc_athr
        self.chain_size["r"] = self.chain_size["swollen"] * self.r_ideal
        self.r_max = self.chain_size["r"].max()
        self.chain_size["r_scaled_max"] = self.chain_size["r"] / self.r_max
        self.chain_size["r_scaled_flory"] = self.chain_size["r"] / self.r_flory
        self.chain_size["w3body_scaled"] = self.chain_size["w3body"] / \
            self.w3body_athr
        if self.w3body_model == "shendruk":
            self.w3body = (
                self.chain_size["w3body"].min(),
                self.chain_size["w3body"].max(),
            )
        if limit:
            # Limit the data only to -1*exc_vol_athr <= exc_vol <= exc_vol_athr
            limit_range = (self.chain_size["vexc_scaled"] <= 1.5) & (
                self.chain_size["vexc_scaled"] >= -1
            )
            self.chain_size = self.chain_size[limit_range]


class LammpsDataTemplate(object):
    name: str
    ispath: bool = True
    """Parse the filename of a LAMMPS template data file based on one of the
    following file name templates:

    1. 'data_template-geometry-group-N#D#L#'

        where 'data_template' is the fixed header name of the LAMMPS template
        files, 'geometry' is the geometry of the space, 'group' is the name of
        file group, 'N" is the number of monomers, 'D' is the diameter (size)
        of the cylinder, and 'L' is the length of cylinder. This template is
        for 'cylindrical' and 'slit' geometries. It is important to note that
        the system is periodic along z direction with range [-L/2,L/2]
        inclusive in the 'cylindrical' geometry while it is periodic along x
        and y directions both with the same range [-D/2,D/2] in the 'slit'
        geometry. Here, x, y, and z are cartesian coordinates.

     2. 'data_template-geometry-group-N#L#'

        where 'data_template' is the fixed header name of the LAMMPS template
        files, 'geometry' is the geometry of the space, 'group' is the name of
        file group, 'N" is the number of monomers, and 'L' is the length
        (size) of the cubic box.


    To-do List
    ----------
    1. It would be great if the user set the input file name template for its
    project.

    Parameters
    ----------
    name: str
        The filename or filepath of the LAMMPS template data file.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    geometry : {'cylindrical', 'slit', 'cubic'}, default 'cylindrical'
        Shape of the simulation box
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    nmon: int, np.nan
        Number of monomers. Its associated keyword is 'N'.
    dcyl: float, np.nan
        Size (or diameter) of the cylindrical (cylindrical) confinement,
        inferred from 'D' keyword keyword. Caution: The size of cylindrical
        confinement is defined based on the LAMMPS' tango, so it is different
        from the size defined in other parsers defined in this module.

    Class Attributes
    ----------------
    _groups: list of str
        Possible groups of the `SumRule` project.
    _geometries_attributes: dict of dict
        a dictionary of `geometry` names. For each `geometry`, a dictionary
        maps the keywords of physical attributes in that geometry to their
        corresponding attributes in `DataTemplate` class.
    """
    _groups = ["bug", "all"]
    _geometries_attributes = {
        "cylindrical": {"nmon": "N", "dcyl": "D", "lcyl": "L"},
        "slit": {"nmon": "N", "dcyl": "D", "lcyl": "L"},
        "cubic": {
            "nmon": "N",
            "lcyl": "L",
        },
    }

    def __init__(self, name: str, ispath: bool = True):
        if ispath:
            self.filepath = name
            self.filename, _ = os.path.splitext(self.filepath)
            self.filename = self.filename.split("/")[-1]
        else:
            self.filepath = "N/A"
            self.filename = name
        self._initiate_attributes()
        self._parse_filename()

    def __str__(self) -> str:
        observation = f"""
        Lammps template data file:
            Name: '{self.filename}',
            Geometry: '{self.geometry},
        """
        return observation

    def __repr__(self) -> str:
        return f"LAMMPS template data('{self.filename}' in geometry \
            '{self.geometry}' from group '{self.group}'."

    def _initiate_attributes(self):
        """
        defines and initiates the class attribute based on the physical
        attributes defined for the project.
        """
        self.geometry = "N/A"
        self.group = "N/A"
        self.nmon = np.nan
        self.dcyl = np.nan
        self.lcyl = np.nan

    def _parse_filename(self):
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        words = self.filename.split("-")
        if len(words) == 4:
            raise ValueError(
                f"The number of words in the filename '{self.filename}'"
                " is not 4 when it is split by '-'. Please check whether"
                " the file name match the template data file name:"
                " 'data_template-geometry-group-N#D#' or not."
            )
        geometry = words[1]
        if geometry in self._geometries_attributes.keys():
            self.geometry = geometry
        else:
            geometries_string = (
                "'" + "', '".join(self._geometries_attributes.keys()) + "'"
            )
            raise ValueError(
                f"'{geometry}' "
                "is not a valid geometry. Please select one of "
                f"{geometries_string} geometries."
            )
        group = words[2]
        if group in self._groups:
            self.group = group
        else:
            groups_string = "'" + "', '".join(self._groups) + "'"
            raise ValueError(
                f"'{group}' "
                "is not a valid particle group. Please select one of "
                f"{groups_string} groups."
            )
        str_attrs = re.compile(r"([a-zA-Z\-]+)")
        attributes_float = ["dcyl", "lcyl"]
        attrs = str_attrs.split(words[3])
        for attr_n, attr_kw in self._geometries_attributes[
            self.geometry
        ].items():
            try:
                attr_value = attrs[attrs.index(attr_kw) + 1]
                if attr_n in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                setattr(self, attr_n, attr_value)
            except ValueError:
                print(
                    f"'{attr_kw}'"
                    " attribute keyword is not in "
                    f"'{self.filename}'"
                    " template data name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not."
                )


@dataclass
class ExcludedVolume(object):
    name: str
    ispath: bool = True
    """Parse a filepath/filename of a 'csv' file that contains the
    excluded-volume data of a monomer in crowded media and retrieve
    information about the exclude-volume data in that 'csv' file.

    The excluded-volume data are currently available for three models:

        AO: The Asakura-Oosawa depletion interaction
            The excluded volume of a large hard-spheres in a solution of
            smaller by using the Asakura-Oosawa depletion potential as the
            effective pair potential between a pair of monomers.

        Edwards: The Edwards-type interaction
            The excluded volume of a monomer in a solution of Edwards
            (self-avoiding) polymers and hard sphere (larger than monomers) is
            found by solving a coarse-grained Hamiltonian of the system.

        LJ: The Lannard-Jones interaction
            The excluded volume of a monomer is given by a function that is
            fitted to the numerically-measured excluded volume of a monomer
            in a bead-on-spring chain. Monomers and crowders are interacting
            by the Lannard-Jones interaction with various cut-off and
            interaction parameters.


    Parameters
    ----------
    name: str
        The name of the filepath or the filename.
    ispath: bool
        Whether the `name` is a filepath or not.

    Attributes
    ----------
    self.filepath: str
        Path to the data file.
    self.filename: str
        Name of the data file.
    self.dmon: float
        Size (diameter) of a monomer
    self.vexc_model: str
        Model by which the excluded-volume data in a crowded media is
        calculated.
    self.dcrowd: float
        Size of a crowder
    self.vexc_athr: float
        The excluded volume of a monomer in the absence of crowders in the high
        temperature (athermal) limit.
    self.vexc_df: pd.DataFrame
        A dataframe that contains the excluded-volume data.

    Class attributes
    ----------------
    self._vexc_models: list of str
        List if the defined excluded-volume models.
    """
    _vexc_models = ["AO", "LJ", "Edwards"]

    def __init__(self, name: str, ispath: bool = True):
        self.vexc_df: pd.DataFrame = pd.DataFrame()
        if ispath:
            self.filepath = name
            self.filename, _ = os.path.splitext(self.filepath)
            self.filename = self.filename.split("/")[-1]
        else:
            self.filepath = "N/A"
            self.filename = name
        self._initiate_attributes()
        self._parse_name()
        self._set_vexc_athermal()

    def _initiate_attributes(self) -> None:
        """defines and initiates the class attributes based on the physical
        attributes defined for the project.
        """
        self.vexc_model = "N/A"
        self.dmon = 1
        self.dcrowd = np.nan

    def _parse_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_attrs = re.compile(r"([a-zA-Z\-]+)")
        words = self.filename.split("-")
        if words[1] in self._vexc_models:
            self.vexc_model = words[1]  # model name
        else:
            vexc_models_string = "'" + "', '".join(self._vexc_models) + "'"
            raise ValueError(
                f"'{words[1]}' "
                "is not a valid excluded-volume model. Please check whether "
                "the data belongs to one of "
                f"{vexc_models_string} model or not."
            )
        new_words = str_attrs.split(words[2])  # find dcrowd
        self.dcrowd = float(new_words[new_words.index("ac") + 1])
        if (self.dcrowd == 1) & (self.vexc_model == "WCA"):
            print(
                f"The excluded data for 'a_c={self.dcrowd}'"
                " is computed via the fitting function for a_c<a in the"
                f"excluded-volume model '{self.vexc_model}'."
            )

    def _set_vexc_athermal(self) -> None:
        """set the athermal excluded-volume of a monomer in the absence of the
        crowders based on a given model.
        """
        _vexc_athrs = {
            "AO": (4 / 3) * np.pi * (self.dmon**3),
            "Edwards": (4 / 3) * np.pi * (self.dmon**3),
            "LJ": 4.40945,
        }
        self.vexc_athr = _vexc_athrs[self.vexc_model]

    def read_data(self) -> None:
        """Read the excluded-volume data"""
        if self.ispath is True:
            self.vexc_df = pd.read_csv(self.filepath, names=[
                                       "phi_c_bulk", "vexc"])
        else:
            raise ValueError(
                "The excluded volume data not found:"
                f"'{self.filename}' is not a valid filepath"
            )

    def scale(self, limit: bool = True) -> None:
        """Rescale the bulk volume fraction of crowders 'phi_c_bulk' and
        excluded-volume 'vexc' and add them as two new columns to
        `self.vexc_df` attribute.

        Parameters
        ----------
        limit: bool, default True
            Whether limit the excluded-volume data to [-1*vexc_athr,
            vexc_athr] or not.
        """
        # phi_c is rescaled as phi_c*a/a_c when a_c << a to gain the maximum
        # depletion free energy:
        self.vexc_df["phi_c_bulk_scaled"] = (
            self.dmon * self.vexc_df["phi_c_bulk"] / self.dcrowd
        )
        self.vexc_df["vexc_scaled"] = self.vexc_df["vexc"] / self.vexc_athr
        # Limit the vexc data to [-1*vexc_athr, vexc_athr]:
        if limit is True:
            self.vexc_df = self.vexc_df[
                (self.vexc_df["vexc_scaled"] <= 1)
                & (self.vexc_df["vexc_scaled"] >= (-1 * 1))
            ]

    def add_model_info(self):
        """Add `self.vexc_model` and `self.dcrowd` data as two new columns to
        the `self.vexc_df` attribute.
        """
        self.vexc_df["vexc_model"] = self.vexc_model
        self.vexc_df["dcrowd"] = self.dcrowd


@dataclass
class FreeEnergyVirial(object):
    name: str
    ispath: bool = True
    """Parse a filepath/filename of a 'csv' file that contains the
    free energy approach data of a polymer in crowded cylindrical confinement
    and retrieve information about the chain size, inner and outer density of
    crowders data in that 'csv' file.

    The free energy approach data are currently available for four models that
    are created by combining two different models for depletion volumes of
    two monomers and two different models for tail interactions between
    monomers.

    Models for the maximum depletion volume:
        de Vries:
        Ha:
    Models for the tail interaction:
        Virial three-body term:
        LJ power-6-law term:

    Issues
    ------
    Currently the class works only for the cylindrical system; what about bulk/
    cubic and slit geometries?

    define "_set_other_attributes" for other parsing classes.

    Parameters
    ----------
    name: str
        The name of the filepath or the filename.
    ispath: bool
        Whether the `name` is a filepath or not.

    Attributes
    ----------
    ?

    Class attributes
    ----------------
    ?
    """
    _tail_models = ["ThreeBody", "LJ"]
    _vdep_models = ["deVries", "Ha"]

    def __init__(self, name: str, ispath: bool = True):
        self.r_chain_df = pd.DataFrame()
        if ispath:
            self.filepath = name
            self.filename, _ = os.path.splitext(self.filepath)
            self.filename = self.filename.split("/")[-1]
        else:
            self.filepath = "N/A"
            self.filename = name
        self._initiate_attributes()
        self._parse_name()
        self._set_other_attributes()

    def _initiate_attributes(self) -> None:
        """defines and initiates the class attributes based on the physical
        attributes defined for the project.
        """
        self.tail_model = "N/A"
        self.vdep_model = "N/A"
        self.dmon = 1
        self.vmon = np.pi * self.dmon**3 / 6
        self.nmon = np.nan
        self.dcyl = np.nan
        self.dcrowd = np.nan
        self.vcrowd = np.nan

    def _parse_name(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_attrs = re.compile(r"([a-zA-Z\-]+)")
        words = self.filename.split("-")
        if words[0] in self._tail_models:
            self.tail_model = words[0]  # model name
        else:
            tail_models_string = "'" + "', '".join(self._tail_models) + "'"
            raise ValueError(
                f"'{words[0]}' "
                "is not a valid tail model. Please check whether "
                "the data belongs to one of "
                f"{tail_models_string} model or not."
            )
        if words[1] in self._vdep_models:
            self.vdep_model = words[1]  # model name
        else:
            vdep_models_string = "'" + "', '".join(self._vdep_models) + "'"
            raise ValueError(
                f"'{words[1]}' "
                "is not a valid tail model. Please check whether "
                "the data belongs to one of "
                f"{vdep_models_string} model or not."
            )
        new_words = str_attrs.split(words[2])  # find dcrowd
        self.nmon = float(new_words[new_words.index("N") + 1])
        self.dcyl = float(new_words[new_words.index("D") + 1])
        self.dcrowd = float(new_words[new_words.index("ac") + 1])

    def _set_other_attributes(self) -> None:
        """set any other possible attributes that should be defined based on
        parsed attributes.
        """
        self.vcrowd = np.pi * self.dcrowd**3 / 6

    def read_data(self) -> None:
        """Read the free energy approach data."""
        if self.ispath is True:
            self.r_chain_df = pd.read_csv(
                self.filepath, names=["rho_c_out", "rho_c_in", "r_chain"]
            )
        else:
            raise ValueError(
                "The free energy approach data not found:"
                f"'{self.filename}' is not a valid filepath"
            )

    def scale(self, phi_c_out_cap: float = 0.45) -> None:
        """Rescale the bulk number density of crowders inside and outside the
        chain-occupying region ('rho_c_out' and 'rho_c_in'), and the chain size
        'r_chain' to create new columns.

        `scale` adds the volume fraction of crowders inside and outside the
        chain-occupying region ('phi_c_out' and 'phi_c_in') and normalized
        chain size 'r_scaled' to `self.r_chain_df` dataset.

        Parameters
        ----------
        phi_c_out_cap: float, default 0.45
            The upper limit over the volume fraction of crowders outside the
            chain-occupying region.
        """
        # phi_c is rescaled as phi_c*a/a_c when a_c is much smaller than a to
        # gain the maximum depletion free energy:
        if (phi_c_out_cap <= 0) or (phi_c_out_cap > 1):
            raise ValueError(
                f"'{phi_c_out_cap}' "
                "should be larger than 0 and smaller or equal to 1."
            )
        self.r_chain_df["phi_c_out"] = \
            self.r_chain_df["rho_c_out"] * self.vcrowd
        self.r_chain_df["phi_c_in"] = \
            self.r_chain_df["rho_c_in"] * self.vcrowd
        r_chain_max = self.r_chain_df["r_chain"].max()
        self.r_chain_df["r_scaled"] = self.r_chain_df["r_chain"] / r_chain_max
        cond = self.r_chain_df["phi_c_out"] <= phi_c_out_cap
        self.r_chain_df = self.r_chain_df.loc[cond, :]
        # Limit the vexc data to [-1*vexc_athr, vexc_athr]:

    def add_model_info(self) -> None:
        """Add the parsed attributed as the new columns to
        the `self.r_chain_df` dataset.
        """
        self.r_chain_df["tail_model"] = self.tail_model
        self.r_chain_df["vdep_model"] = self.vdep_model
        self.r_chain_df["nmon"] = self.nmon
        self.r_chain_df["dcyl"] = self.dcyl
        self.r_chain_df["dcrowd"] = self.dcrowd


class Snapshot:
    file: InputT
    """
    Read a single snapshot from a dump `file` object.
    """

    def __init__(self, file: InputT) -> None:
        self.file = file
        # time, n_atoms, and n_props are default to invalid numbers.
        self.time = -1  # time stamp
        self.n_atoms = 0  # number of atoms
        self.boxstr = ''  # box information
        self.is_scaled = -1  # unknown status for coordinates
        self.triclinic: bool = False  # default: orthogonal box
        self.xlo = 0.0  # box information
        self.xhi = 0.0
        self.xy = 0.0
        self.ylo = 0.0
        self.yhi = 0.0
        self.xz = 0.0
        self.zlo = 0.0
        self.zhi = 0.0
        self.yz = 0.0
        # names and col indices of per-atom properties:
        self.props: Dict[str, int] = {}
        self.n_props = 0  # number of per atom properties in the dump file
        self.atoms = np.empty((self.n_atoms, self.n_props), dtype=np.float64)
        self._read()

    def _read(self) -> None:
        """
        Read a single snapshot and assigns column names (file must be
        self-describing). Moreover, it changes `self.is_scaled` from -1
        (unknown) to 0 (unscaled) or 1 (scaled), and converts 'xs' and 'xu' to
        'x' in `self.props`.
        """
        _ = self.file.readline()
        # just grab 1st field
        self.time = int(self.file.readline().split()[0])
        _ = self.file.readline()
        self.n_atoms = int(self.file.readline())
        self.boxstr = self.file.readline()
        words = self.boxstr.split()
        self.triclinic = words == 9  # ITEM: BOX BOUNDS
        if self.triclinic:
            self.xlo, self.xhi, self.xy = \
                map(float, self.file.readline().split())
            self.ylo, self.yhi, self.xz = \
                map(float, self.file.readline().split())
            self.zlo, self.zhi, self.yz = \
                map(float, self.file.readline().split())
        else:
            self.xlo, self.xhi = map(float, self.file.readline().split())
            self.ylo, self.yhi = map(float, self.file.readline().split())
            self.zlo, self.zhi = map(float, self.file.readline().split())

        x_flag = y_flag = z_flag = -1
        self.props = {
            c: i for i, c in enumerate(str(self.file.readline().split()[2:]))
            }  # dump per atom props
        props = list(self.props.keys())
        if ("x" in props) or ("xu" in props):
            x_flag = 0
        elif ("xs" in props) or ("xsu" in props):
            x_flag = 1
        elif ("y" in props) or ("yu" in props):
            y_flag = 0
        elif ("ys" in props) or ("ysu" in props):
            y_flag = 1
        elif ("z" in props) or ("zu" in props):
            z_flag = 0
        elif ("zs" in props) or ("zsu" in props):
            z_flag = 1
        if x_flag == 0 and y_flag == 0 and z_flag == 0:
            self.is_scaled = 0
        if x_flag == 1 and y_flag == 1 and z_flag == 1:
            self.is_scaled = 1
        words = self.file.readline().split()
        self.n_props = len(words)
        self.atoms = np.zeros((self.n_atoms, self.n_props), dtype=np.float64)
        self.atoms[0, :] = np.array(list(map(float, words)))  # read first atom
        if self.n_props != len(self.props):
            raise ValueError(
                f"Number of per atom columns '{self.n_props}' is not equal "
                f"to the number of column names '{len(self.props)}'."
                )
        for i in range(1, self.n_atoms):
            self.atoms[i, :] = \
                np.array(list(map(float, self.file.readline().split())))


class Dump:
    """
    Read, write, and manipulate dump files and particle attributes.

    Examples
    --------
    d = dump("dump.one")              read in a single dump.
    d = dump(["dump.1", "dump.2"])	  can be more than one dump.
    d = dump(["dump.1", "dump.2.gz"]) can be gzipped

    `Dump` automatically delete incomplete and duplicate snapshots and
    unscaled coordinates if they are stored in files as scaled.
    """
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.names: Dict[str, int] = {}
        self.is_scaled = -1  # -1/0/1 mean unknown/unscale/scale coordinates.
        self.snaps: List[Snapshot] = []
        self.n_snaps: int = 0  # total number of snapshots
        self.increment = 1
        self.eof = 0

    def read_all(self) -> None:
        """
        Read all snapshots from each file; test for gzipped files.
        """
        with openany_context(self.filepath) as dfile:
            snap = Snapshot(dfile)
            if self.names == {}:
                self.names = snap.props
                print("Assigned columns: " + self.names2str())
                self.is_scaled = snap.is_scaled
            while snap:
                self.snaps.append(snap)
                snap = Snapshot(dfile)
        # sort entries by timestep, cull duplicates
        self.snaps.sort(key=self.compare_time)
        self.cull()
        self.n_snaps = len(self.snaps)
        print(f"{self.n_snaps} snapshots were read.")
        # select all timesteps and atoms
        # if snapshots are scaled, unscale them
        if ("x" not in self.names.keys()) or \
            ("y" not in self.names.keys()) or \
                ("z" not in self.names.keys()):
            print("Dump scaling status is unknown.")
        elif self.n_snaps > 0:
            if self.is_scaled == 1:
                self.unscale()
            elif self.is_scaled == 0:
                print("dump is already unscaled")
        else:
            print("Dump scaling status is unknown")

    def names2str(self) -> str:
        """
        Convert column names to a string, by column order.

        Returns
        -------
        col_str: str
            The string of column names. `col_str` is repeated in the header
            section of a snapshot, determining per atom properties in that
            snapshot.
        """
        names_res: Dict[int, str] = \
            dict(((v, k) for k, v in self.names.items()))
        names_res = OrderedDict(sorted(names_res.items()))
        col_str = " ".join(names_res.values())
        return col_str

    def cull(self) -> None:
        """
        Delete successive snapshots with duplicate time stamp.
        """
        i = 1
        while i < len(self.snaps):
            if self.snaps[i].time == self.snaps[i-1].time:
                del self.snaps[i]
            else:
                i += 1

    def find_snapshot(self, ts: int) -> int:
        """
        Find a snapshot with timestep `ts`.

        Parameters
        ----------
        ts: int
            The snapshot to be found.

        Returns
        -------
        i: int
            The index of the snapshot in `self.n_snaps` for which timestep is
            equal to `ts`.

        Raises
        ------
        ValueError
            raise error  if `ts` is not found.
        """
        for i in range(self.n_snaps):
            if self.snaps[i].time == ts:
                return i
        raise ValueError(f"No step '{ts}' exists.")

    @staticmethod
    def compare_time(snap: Snapshot) -> int:
        """
        Sort snapshots on time stamp.
        """
        return snap.time

    @staticmethod
    def unit_cell(snap: Snapshot) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute elements of the h-matrix for a tricilinic or orthogonal unit
        cell and its inverse matrix
        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.

        Return
        ------
        h_mat : numpy array
            An array of length 6 containing the elements of h_matrix.
        h_inv : numpy array
            An array of length 6 containing the elements of the inverse of
            the h_matrix.
        """
        xlo_bound = snap.xlo
        xhi_bound = snap.xhi
        ylo_bound = snap.ylo
        yhi_bound = snap.yhi
        zlo_bound = snap.zlo
        zhi_bound = snap.zhi
        xy = snap.xy
        xz = snap.xz
        yz = snap.yz
        xlo = xlo_bound - min((0.0, xy, xz, xy+xz))
        xhi = xhi_bound - max((0.0, xy, xz, xy+xz))
        ylo = ylo_bound - min((0.0, yz))
        yhi = yhi_bound - max((0.0, yz))
        zlo = zlo_bound
        zhi = zhi_bound
        h_mat = np.zeros(6, dtype=np.float64)
        h_mat[0] = xhi - xlo
        h_mat[1] = yhi - ylo
        h_mat[2] = zhi - zlo
        h_mat[3] = yz
        h_mat[4] = xz
        h_mat[5] = xy
        h_inv = np.zeros(6, dtype=np.float64)
        h_inv[0] = 1.0 / h_mat[0]
        h_inv[1] = 1.0 / h_mat[1]
        h_inv[2] = 1.0 / h_mat[2]
        h_inv[3] = yz / (h_mat[1]*h_mat[2])
        h_inv[4] = (h_mat[3]*h_mat[5] - h_mat[1]*h_mat[4]) / \
            (h_mat[0]*h_mat[1]*h_mat[2])
        h_inv[5] = xy / (h_mat[0]*h_mat[1])
        return h_mat, h_inv

    def scale(self, timestep: int = -1) -> None:
        """
        Scale coordinates from box size to [0,1] for all snapshots or for
        snapshot `timestep`, using the h-matrix to treat orthogonal or
        triclinic boxes.

        Parameters
        ----------
        timestep : int, optional
            The timestep of the snapshot being unscaled.
        """
        if timestep == -1:
            print("Scaling dump ...")
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            for snap in self.snaps:
                self.scale_one(snap, x, y, z)
        else:
            i = self.find_snapshot(timestep)
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            self.scale_one(self.snaps[i], x, y, z)

    def scale_one(self, snap: Snapshot, x: int, y: int, z: int) -> None:
        """
        Scale a single snapshot.

        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.
        x : int
            Column id of the x coordinate.
        y : int
            Column id of the y coordinate.
        z : int
            Column id of the z coordinate.
        """
        # scale coordinates in a single `snap`, using
        if not snap.triclinic:
            xprdinv = 1.0 / (snap.xhi - snap.xlo)
            yprdinv = 1.0 / (snap.yhi - snap.ylo)
            zprdinv = 1.0 / (snap.zhi - snap.zlo)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = (atoms[:, x] - snap.xlo) * xprdinv
                atoms[:, y] = (atoms[:, y] - snap.ylo) * yprdinv
                atoms[:, z] = (atoms[:, z] - snap.zlo) * zprdinv
        else:
            # Creating h-matrix:
            _, h_inv = self.unit_cell(snap)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = (atoms[:, x] - snap.xlo)*h_inv[0] + \
                    (atoms[:, y] - snap.ylo)*h_inv[5] + \
                    (atoms[:, z] - snap.zlo)*h_inv[4]
                atoms[:, y] = (atoms[:, y] - snap.ylo)*h_inv[1] + \
                    (atoms[:, z] - snap.zlo)*h_inv[3]
                atoms[:, z] = (atoms[:, z] - snap.zlo)*h_inv[2]

    def unscale(self, timestep: int = -1) -> None:
        """
        Unscale coordinates from [0,1] to box size for all snapshots or
        for snapshot `timestep`, using the h-matrix to treat orthogonal or
        triclinic boxes.

        Parameters
        ----------
        timestep : int, optional
            The timestep of the snapshot being unscaled.
        """
        if timestep == -1:
            print("Unscaling dump ...")
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            for snap in self.snaps:
                self.unscale_one(snap, x, y, z)
        else:
            i = self.find_snapshot(timestep)
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            self.unscale_one(self.snaps[i], x, y, z)

    def unscale_one(self, snap: Snapshot, x: int, y: int, z: int) -> None:
        """
        Unscale a single snapshot.

        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.
        x : int
            Column id of the x coordinate.
        y : int
            Column id of the y coordinate.
        z : int
            Column id of the z coordinate.
        """
        if not snap.triclinic:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x]*xprd
                atoms[:, y] = snap.ylo + atoms[:, y]*yprd
                atoms[:, z] = snap.zlo + atoms[:, z]*zprd
        else:
            # Creating h-matrix:
            h_mat, _ = self.unit_cell(snap)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x] * \
                    h_mat[0] + atoms[:, y]*h_mat[5] + atoms[:, z]*h_mat[4]
                atoms[:, y] = snap.ylo + atoms[:, y]*h_mat[1] + \
                    atoms[:, z]*h_mat[3]
                atoms[:, z] = snap.zlo + atoms[:, z]*h_mat[2]

    def wrap(self) -> None:
        """
        Wraps coordinates from outside box to inside.
        """
        print("Wrapping dump ...")
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        ix = self.names["ix"]
        iy = self.names["iy"]
        iz = self.names["iz"]
        for snap in self.snaps:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            atoms[:, x] -= atoms[:, ix]*xprd
            atoms[:, y] -= atoms[:, iy]*yprd
            atoms[:, z] -= atoms[:, iz]*zprd

    def unwrap(self) -> None:
        """
        Unwraps coordinates from inside box to outside.
        """
        print("Unwrapping dump ...")
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        ix = self.names["ix"]
        iy = self.names["iy"]
        iz = self.names["iz"]
        for snap in self.snaps:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            atoms[:, x] += atoms[:, ix]*xprd
            atoms[:, y] += atoms[:, iy]*yprd
            atoms[:, z] += atoms[:, iz]*zprd

    @staticmethod
    def write_header(snap: Snapshot, file: IO, cols: str) -> None:
        """
        Write the 'header' section of a dump file.

        Parameters
        ----------
        snap : Snapshot
            The snapshot for which header is written
        file : TextIO
            Pointer to the file to which the header is written
        cols : str
            The space-separated column names
        """
        file.write("ITEM: TIMESTEP\n")
        file.write(str(snap.time) + "\n")
        file.write("ITEM: NUMBER OF ATOMS\n")
        file.write(str(snap.n_atoms) + "\n")
        if snap.boxstr:
            file.write(snap.boxstr)
        else:
            file.write("ITEM: BOX BOUNDS\n")
        if snap.triclinic:
            file.write(f"{snap.xlo} {snap.xhi} {snap.xy}\n")
            file.write(f"{snap.ylo} {snap.yhi} {snap.xz}\n")
            file.write(f"{snap.zlo} {snap.zhi} {snap.yz}\n")
        else:
            file.write(f"{snap.xlo} {snap.xhi}\n")
            file.write(f"{snap.ylo} {snap.yhi}\n")
            file.write(f"{snap.zlo} {snap.zhi}\n")
        file.write("ITEM: ATOMS " + cols + "\n")

    def write(self, filename: str, mode: str = "w"
              ) -> None:
        """
        write a single dump file from current selection.

        Parameters
        ----------
        filename: str
            _description_
        mode : str, optional
            Whether write to a new file or append to an existing file.
        """
        col_str = self.names2str()
        if "id" in self.names:
            atom_id_col = self.names["id"]
        else:
            raise ValueError("atom id column not found!")
        if "type" in self.names:
            atom_type_col = self.names["type"]
        else:
            raise ValueError("atom type column not found!")
        with open(filename, mode) as snapshot:
            for snap in self.snaps:
                self.write_header(snap, snapshot, col_str)
                for atom in snap.atoms:
                    line = ""
                    for j in range(snap.n_props):
                        if j == atom_id_col or j == atom_type_col:
                            line += str(int(atom[j])) + " "
                        else:
                            line += str(atom[j]) + " "
                    snapshot.write(line)
                    snapshot.write("\n")
        print(f"{self.n_snaps} snapshots are written to a single dump file.")

    def scatter(self, prefix: str) -> None:
        """
        Write one dump file per snapshot from the current selection, using
        `prefix` as the prefix for filenames.

        Parameters
        ----------
        prefix: str
            The prefix used in the filenames.
        """
        col_str = self.names2str()
        for snap in self.snaps:
            print(snap.time, flush=True)
            filename = prefix + "." + str(snap.time) + '.lammpstrj'
            with open(filename, "w") as snapshot:
                self.write_header(snap, snapshot, col_str)
                for atom in snap.atoms:
                    line = ""
                    for j in range(snap.n_props):
                        if j < 2:
                            line += str(int(atom[j])) + " "
                        else:
                            line += str(atom[j]) + " "
                    snapshot.write(line)
        print(f"{self.n_snaps} snapshot(s) are written to "
              f"{self.n_snaps} file(s).")
