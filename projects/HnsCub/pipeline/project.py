import os
import re
import warnings
from typing import Literal, TypeVar, IO, Tuple, Dict, List, Union
from abc import ABC, abstractmethod
from dataclasses import dataclass
from collections import OrderedDict
import numpy as np
import pandas as pd

from polyphys.manage.utilizer import invalid_keyword, openany_context, InputT
from polyphys.manage.parser import ParserBase


class HnsCub(ParserBase):
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
        space: N#kbmm#nh#ac#epshc#ring
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

    def __init__(
        self,
        name: str,
        lineage: str,
        lineage_attributes: Dict[str,  Dict[str, str]],
        physical_attributes: List[str],
        geometry: str,
        group: str,
        topology: str
    ) -> None:
        super().__init__(name, lineage, lineage_attributes)
        self._geometry = geometry
        self._group = group
        self._topology = topology
        self._extract_attribute_from_filename()
        self._set_parents()
        self._physical_attributes = physical_attributes
        self.attributes = (
            self._lineage_attributes[self.lineage]
            + self._physical_attributes[self.lineage]
        )

        if self.lineage in ["segment", "whole", "ensemble_long"]:
            self._bulk_attributes()

    def _extract_attribute_from_filename(self) -> None:
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        str_lineages = re.compile(r"([a-zA-Z\-]+)")
        words = str_lineages.split(self.filename)
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
        for nama in self.GENEALOGY[self.lineage]:
            
        
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