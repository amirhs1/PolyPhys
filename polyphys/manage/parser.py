import numpy as np
import re


class SumRule(object):
    name: str
    geometry: str = 'biaxial'
    group: str = 'bug'
    lineage: str = 'segment'
    ispath: bool = True
    """
    parses a `lineage_name` to extract information about the 'lineage' of \
    that 'lineag_name', based on the following 'lineage' patterns:

    segment: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
        One of multiple chunks of a complete simulation or measurement.
    whole: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
        A complete simulation or measurement; a collection of 'segments'.
    ensemble: N#D#ac#nc#
        A collection of 'wholes' (complete simulations) that differs only \
        in their initial conditions (e.g., random number seed).
    ensemble_long: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
        Long name of an ensemble.
    space: N#D#ac#
        A collection of ensembles.

    In the above lineages, the keywords are attributes where their values \
    (shown by "#" sign) are float or integer number. If a lineage does not \
    have an attribute, then the value of that attribute is set to numpy.nan.
    These are the attributes with "numpy.nan" values for different lineages:
    whole: 'j'
    ensemble_long: 'ens', and 'j'
    ensemble: 'lz', 'dt', 'bdump', 'adump', 'ens', and 'j'
    space:  'lz', 'nc', 'dt', 'bdump', 'adump', 'ens' , and 'j'

    There are some diffeence between the keywords of physical attributes and \
    their associated attributes in the `SumeRule` class. Below, the these two \
    types of attributes are explained.

    To-do List
    ----------
    1. This class can be split to 4 classes in the following order of \
    inheritance: parent->child: space -> ensemble -> whole -> segment.
    2. Using @attribute for attributes to get, set, or delete them.

    Parameters
    ----------
    name: str
        Name that is parsed for extracting information.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group. 'bug' is used for a single polymer. \
        'all' is used for all the particles/atoms in the system.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble', \
        'space'}, default 'segment'
        Type of the lineage of the name.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    group: {'bug', 'all'}
        Type of the particle group.  'bug' is used for a single polymer. \
        'all' is used for all the particles/atoms in the system.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble', \
        'space'}, default whole
        Type of the lineage of the name.
    pathname: str, default "N/A"
        Equal to `name` if `name` is a filepath, otherwise "N/A".
    filename: str
        Name of a the file reffered to by `name` if `name` is a filepath, \
        otherwise the `name` itself.
    lineage_name: str,
        The unique name of type extracted from self.fullname
    dmon: float, default 1.0
        Size (diameter) of a monomer
    nmon: int, np.nan
        Number of monomers. Its associated keyword is 'N'.
    dcyl: float, np.nan
        Size (or diameter) of the biaxial confinement, inferred from either \
        'r' keyword (radius of the biaxial confinement, a cylinder with open \
        ends) or 'D' keyword (size of that confinement
    lcyl: float, np.nan
        Length of the biaxial confinement along z axis (the periodic, \
        direction), inferred from 'lz' keyword (half of the length of the \
        biaxial confinement along z axis.
    epsilon: float, np.nan
        Wall-particle LJ interaction strength. Its associated keyword is \
        'epsilon' keyword.
    dcrowd: float, np.nan
        Size (diamter) of a crowder. Its associated keyword is 'sig' or 'ac'.
    ncrowd: int, np.nan
        number of crowders
    ensemble_id: int, np.nan
        The ensemble number of a 'whole' simulation in an ensemble.
    segment_id: int, np.nan
        The 'segment_id' keyword starts with 'j', ends with a 'padded' \
        number such as '05' or '14', showing the succesion of segments \
        in a whole file.
    dt: float, np.nan
        Simulation timestep. Its associated keyword is 'dt'.
    bdump: int
        Frequency by which 'bug' configurations are dumped in a 'bug' \
        trajectory file. Its associated keyword is 'bdump'.
    adump: int, default np.nan
        Frequency by which 'all' configurations are dumped in a 'segment' \
        trajectory file. Its associated keyword is 'adump'.
    mmon: float, default 1.0
        Mass of a monomer
    eps_others: float, default 1.0
        Other LJ interaction strengths
    mcrowd: float, default 1.0
        Mass of a crowder
    dwall: float, default 1.0
        Wall-forming particles diameter
    space: str
        A space's name
    ensemble: str or "N/A"
        An ensemble's name if applicable, otherwise "N/A"
    ensemble_long: str or "N/A"
        The name of ensemble derived from 'whole' name if applicable, \
        otherwise "N/A"
    whole: str or "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment: str or "N/A"
        A segment's name if applicable, otherwise "N/A"
    self.rho_m_bulk: float, default np.nan
        Bulk number desnity fraction of monomers
    self.phi_m_bulk: float, default np.nan
        Bulk volume fraction of monomers
    self.rho_c_bulk: float, default np.nan
        Bulk number density fraction of crowders
    self.phi_c_bulk: float, default np.nan
        Bulk volume fraction of crowders

    Class Attributes
    ----------------
    _geomteries: list of str
        Possible geometries of a simulation box
    _groups: list of str
        Possible groups of the `SumRule` project.
    _lineage_attributes: dict of dict
        a dictionary of `lineage` names. For each `lineage`, a dictionary \
        maps the keywords of physical attributes in that lineage to their \
        corresponding attributes in `SumRule` class.
    _lineage_private_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of \
        class attributes that are NOT "N/A" or np.nan is created that can be\
        used in the simulation/experiment/run report files (called \
        "*-properties.csv")
    """
    _geometries = ['biaxial', 'slit', 'box']
    _groups = ['bug', 'all']
    _lineage_attributes = {
        'segment': {
            'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
            'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'ensemble_id': 'ens', 'segment_id': 'j'
            },
        'whole': {
            'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
            'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump', 'ensemble_id': 'ens'
            },
        'ensemble_long': {
            'nmon': 'N', 'epsilon': 'epsilon', 'dcyl': 'r', 'lcyl': 'lz',
            'dcrowd': 'sig', 'ncrowd': 'nc', 'dt': 'dt', 'bdump': 'bdump',
            'adump': 'adump'
            },
        'ensemble': {
            'nmon': 'N', 'dcyl': 'D', 'dcrd': 'ac', 'ncrowd': 'nc'
            },
        'space': {
            'nmon': 'N', 'dcyl': 'D', 'dcrowd': 'ac'
            }
        }
    _physical_attributes = {
        'segment': [
            'dmon', 'phi_m_bulk', 'rho_m_bulk', 'phi_c_bulk', 'rho_c_bulk'
            ],
        'whole': [
            'dmon', 'phi_m_bulk', 'rho_m_bulk', 'phi_c_bulk', 'rho_c_bulk'
            ],
        'ensemble_long': [
            'dmon', 'phi_m_bulk', 'rho_m_bulk', 'phi_c_bulk', 'rho_c_bulk'
            ],
        'ensemble': [
            'dmon'
            ],
        'space': [
            'dmon'
            ]
        }
    _genealogy = {
        'segment': [
            'lineage_name', 'segment', 'whole', 'ensemble_long',
            'ensemble', 'space'
            ],
        'whole': [
            'lineage_name', 'whole', 'ensemble_long', 'ensemble',
            'space'
            ],
        'ensemble_long': [
            'lineage_name', 'ensemble_long', 'ensemble', 'space',
            ],
        'ensemble': [
            'lineage_name', 'ensemble', 'space'
            ],
        'space': [
            'lineage_name', 'space'
            ]
        }

    def __init__(
        self,
        name: str,
        geometry: str = 'biaxial',
        group: str = 'bug',
        lineage: str = 'segment',
        ispath: bool = True
    ):
        if geometry in self._geometries:
            self.geometry = geometry
        else:
            geometries_string = "'" + "', '".join(
                self._geometries) + "'"
            raise ValueError(
                f"'{geometry}' "
                "is not a valid geometry. Please select one of "
                f"{geometries_string} geometries.")
        if group in self._groups:
            self.group = group
        else:
            groups_string = "'" + "', '".join(
                self._groups) + "'"
            raise ValueError(
                f"'{group}' "
                "is not a valid particle group. Please select one of "
                f"{groups_string} groups.")
        if lineage in self._lineage_attributes.keys():
            self.lineage = lineage
        else:
            types_string = "'" + "', '".join(
                self._lineage_attributes.keys()) + "'"
            raise ValueError(
                f"'{type}' "
                "is not a valid name type. Please select one of "
                f"{types_string} types.")
        if ispath:
            self.filepath = name
            self.filename = name.split("/")[-1]
        else:
            self.filepath = "N/A"
            self.filename = name
        self._find_lineage_name()
        self._initiate_attributes()
        self._parser_lineage_name()
        self._set_parents()
        self._bulk_attributes()
        self.attributes = list(self._lineage_attributes[self.lineage].keys())\
            + self._physical_attributes[self.lineage]
        self.genealogy = self._genealogy[self.lineage]

    def __str__(self) -> str:
        observation = f"""
        Observation:
            Name: '{self.filename}',
            Geometry: '{self.geometry},
            Group: '{self.group}',
            Lineage: '{self.lineage}')
        """
        return observation

    def __repr__(self) -> str:
        return f"Observation('{self.filename}' in geometry '{self.geometry} \
            from group '{self.group}' with lineage '{self.lineage}'"

    def _find_lineage_name(self):
        """
        parses the unique lineage_name (the first substring of filename \
        and/or the segment keyword middle substring) of a filename.
        """
        if self.lineage in ['segment', 'whole']:
            # a 'segment' lineage only used in 'probe' phase
            # a 'whole' lineage used in 'probe' or 'analyze' phases\
            # so its lineage_name is either ended by 'group' keyword or "-".
            # these two combined below:
            self.lineage_name = \
                self.filename.split("." + self.group)[0].split("-")[0]
        else:  # 'ensemble' or 'space' lineages
            self.lineage_name = self.filename.split('-')[0]

    def _initiate_attributes(self):
        """
        defines and initiates the class attribute based on the physical \
        attributes defined for the project.
        """
        self.dmon = 1.0
        self.nmon = np.nan
        self.dcyl = np.nan
        self.lcyl = np.nan
        self.epsilon = np.nan
        self.dcrowd = np.nan
        self.ncrowd = np.nan
        self.ensemble_id = np.nan
        self.segment_id = np.nan
        self.dt = np.nan
        self.bdump = np.nan
        self.adump = np.nan
        self.mmon = 1.0
        self.eps_others = 1.0
        self.mcrowd = 1.0
        self.dwall = 1.0
        self.phi_m_bulk = np.nan
        self.rho_m_bulk = np.nan
        self.phi_c_bulk = np.nan
        self.rho_c_bulk = np.nan

    def _parser_lineage_name(self):
        """
        parses a lineage_name based on a list of keywords of physical \
        attributes.
        """
        str_lineages = re.compile(r'([a-zA-Z\-]+)')
        words = str_lineages.split(self.lineage_name)
        attributes_float = ['dmon', 'dcyl', 'lcyl', 'epsilon', 'dcrowd', 'dt']
        for attr_name, attr_keyword in \
                self._lineage_attributes[self.lineage].items():
            try:
                attr_value = words[words.index(attr_keyword)+1]
                if attr_name in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                if attr_keyword == 'lz':
                    attr_value = 2 * attr_value
                if attr_keyword == 'r':
                    attr_value = 2 * attr_value - 1
                setattr(self, attr_name, attr_value)
            except ValueError:
                print(
                    f"'{attr_keyword}'"
                    " attribute keyword is not in "
                    f"'{self.lineage_name}'"
                    " lineage name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not.")

    def _set_parents(self):
        """
        set to parent names for a lineage_name based on it lineage.


        The following map is used for seeting relationships:
        'segment' lineage is a child of 'whole' lineage.
        'whole' lineage is a child of 'ensemble' lineage.
        'ensemble' lineage is a child of 'space' lineage.
        'space' is the root of other lineages.
        """
        if self.lineage == 'space':
            self.space = self.lineage_name
            self.ensemble = "N/A"
            self.ensemble_long = "N/A"
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == 'ensemble':
            self.space = self.lineage_name.split('nc')[0]
            self.ensemble = self.lineage_name
            self.ensemble_long = "N/A"
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == 'ensemble_long':
            self.space = 'N' + str(self.nmon) + 'D' + str(self.dcyl) \
                + 'ac' + str(self.dcrowd)
            self.ensemble = self.space + 'nc' + str(self.ncrowd)
            self.ensemble_long = self.lineage_name
            self.whole = "N/A"
            self.segment = "N/A"
        elif self.lineage == 'whole':
            self.space = 'N' + str(self.nmon) + 'D' + str(self.dcyl) \
                + 'ac' + str(self.dcrowd)
            self.ensemble = self.space + 'nc' + str(self.ncrowd)
            self.ensemble_long = self.lineage_name.split('ens')[0]
            self.whole = self.lineage_name
            self.segment = "N/A"
        else:
            self.space = 'N' + str(self.nmon) + 'D' + str(self.dcyl) \
                + 'ac' + str(self.dcrowd)
            self.ensemble = self.space + 'nc' + str(self.ncrowd)
            self.ensemble_long = self.lineage_name.split('ens')[0]
            self.whole = self.lineage_name.split(".j")[0]
            self.segment = self.lineage_name

    def _bulk_attributes(self):
        """
        computes some physicla attributes of a lineage based on its \
        primary attributes.
        """
        if self.lineage in ['segment', 'whole']:
            vol_cell = np.pi * self.dcyl**2 * self.lcyl / 4.0
            vol_mon = np.pi * self.dmon**3 / 6
            self.rho_m_bulk = self.nmon / vol_cell
            self.phi_m_bulk = self.rho_m_bulk * vol_mon
            vol_crowd = np.pi * self.dcrowd**3 / 6
            self.rho_c_bulk = self.ncrowd / vol_cell
            self.phi_c_bulk = self.rho_c_bulk * vol_crowd
