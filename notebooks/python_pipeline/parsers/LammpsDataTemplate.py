import os
import re
import numpy as np


class LammpsDataTemplate:
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
        if geometry in self._geometries_attributes:
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

