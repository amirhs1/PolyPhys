import os
import re
from typing import TypeVar
from dataclasses import dataclass
import numpy as np
import pandas as pd

TFreeEnergyVirial = TypeVar("TFreeEnergyVirial", bound="FreeEnergyVirial")

@dataclass
class FreeEnergyVirial:
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
