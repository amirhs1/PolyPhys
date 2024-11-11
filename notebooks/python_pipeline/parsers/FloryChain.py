import os
import re
from typing import TypeVar
import numpy as np
import pandas as pd

TExcludedVolume = TypeVar("TExcludedVolume", bound="ExcludedVolume")
TFreeEnergyVirial = TypeVar("TFreeEnergyVirial", bound="FreeEnergyVirial")


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
            "dcrowd": "ac",
        },
    }
    _column_names = {
        "deGennes": ["phi_c", "vexc", "swollen"],
        "twoBodyOnly": ["phi_c", "vexc", "swollen"],
        "shendruk": ["phi_c", "vexc", "w3body", "swollen"],
    }

    def __init__(
        self, name: str, limit: bool = False, ispath: bool = True
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
        self.chain_size["vexc_scaled"] = (
            self.chain_size["vexc"] / self.vexc_athr
        )
        self.chain_size["r"] = self.chain_size["swollen"] * self.r_ideal
        self.r_max = self.chain_size["r"].max()
        self.chain_size["r_scaled_max"] = self.chain_size["r"] / self.r_max
        self.chain_size["r_scaled_flory"] = self.chain_size["r"] / self.r_flory
        self.chain_size["w3body_scaled"] = (
            self.chain_size["w3body"] / self.w3body_athr
        )
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
