import os
import re
from typing import TypeVar
from dataclasses import dataclass
import numpy as np
import pandas as pd

TExcludedVolume = TypeVar("TExcludedVolume", bound="ExcludedVolume")


@dataclass
class ExcludedVolume:
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
        self.ispath = ispath
        if self.ispath:
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
