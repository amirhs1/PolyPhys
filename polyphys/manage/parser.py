import numpy as np
import pandas as pd
import os
import re


class SumRule(object):
    name: str
    geometry: str = 'biaxial'
    group: str = 'bug'
    lineage: str = 'segment'
    ispath: bool = True
    """
    parses a `lineage_name` to extract information about the 'lineage' oftrj
    that 'lineage_name', based on the following 'lineage' patterns:

    segment: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#.j#
        One of multiple chunks of a complete simulation or measurement.
    whole: N#epsilon#r#lz#sig#nc#dt#bdump#adump#ens#
        A complete simulation or measurement; a collection of 'segments'.
    ensemble: N#D#ac#nc#
        A collection of 'wholes' (complete simulations) that differs only
        in their initial conditions (e.g., random number seed).
    ensemble_long: N#epsilon#r#lz#sig#nc#dt#bdump#adump#
        Long name of an ensemble.
    space: N#D#ac#
        A collection of ensembles.

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

    To-do List
    ----------
    1. This class can be split to 4 classes in the following order of
    inheritance: parent->child: space -> ensemble -> whole -> segment.
    2. Using @attribute for attributes to get, set, or delete them.

    Parameters
    ----------
    name: str
        Name that is parsed for extracting information.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group. 'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}, default 'segment'
        Type of the lineage of the name.
    ispath: bool, default True
        Whether the name is a filepath or a simple name.

    Attributes
    ----------
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    group: {'bug', 'all'}
        Type of the particle group.  'bug' is used for a single polymer.
        'all' is used for all the particles/atoms in the system.
    lineage: {'segment', 'whole', 'ensemble_long', 'ensemble',
        'space'}, default whole
        Type of the lineage of the name.
    pathname: str, default "N/A"
        Equal to `name` if `name` is a filepath, otherwise "N/A".
    filename: str
        Name of a the file referred to by `name` if `name` is a filepath,
        otherwise the `name` itself.
    lineage_name: str,
        The unique name of type extracted from self.fullname
    dmon: float, default 1.0
        Size (diameter) of a monomer
    nmon: int, np.nan
        Number of monomers. Its associated keyword is 'N'.
    dcyl: float, np.nan
        Size (or diameter) of the `biaxial` or `slit` confinement, inferred
        from either 'r' keyword (radius of the biaxial confinement; a cylinder
        with open ends) or 'D' keyword (size of that confinement. Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl.2] inculsive is the
        domain over which x and y cartesian coordiantes are defined in the
        'biaxial' geometry; however, `dcyl` ranged from '[-dcyl/2,dcyl.2]
        inculsive is the domain over which z cartesian coordiante is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
        Hence, `dcyl` is defined differenty in different parser classes in
        this module.
    lcyl: float, np.nan
        Length of the biaxial confinement along z axis (the periodic,
        direction), inferred from 'lz' keyword (half of the length of the
        biaxial confinement along z axis.
    epsilon: float, np.nan
        Wall-particle LJ interaction strength. Its associated keyword is
        'epsilon' keyword.
    dcrowd: float, np.nan
        Size (diameter) of a crowder. Its associated keyword is 'sig' or 'ac'.
    ncrowd: int, np.nan
        number of crowders
    ensemble_id: int, np.nan
        The ensemble number of a 'whole' simulation in an ensemble.
    segment_id: int, np.nan
        The 'segment_id' keyword starts with 'j', ends with a 'padded'
        number such as '05' or '14', showing the succession of segments
        in a whole file.
    dt: float, np.nan
        Simulation timestep. Its associated keyword is 'dt'.
    bdump: int
        Frequency by which 'bug' configurations are dumped in a 'bug'
        trajectory file. Its associated keyword is 'bdump'.
    adump: int, default np.nan
        Frequency by which 'all' configurations are dumped in a 'segment'
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
        The name of ensemble derived from 'whole' name if applicable,
        otherwise "N/A"
    whole: str or "N/A"
        A whole's name if applicable, otherwise "N/A"
    segment: str or "N/A"
        A segment's name if applicable, otherwise "N/A"
    self.rho_m_bulk: float, default np.nan
        Bulk number density fraction of monomers
    self.phi_m_bulk: float, default np.nan
        Bulk volume fraction of monomers
    self.rho_c_bulk: float, default np.nan
        Bulk number density fraction of crowders
    self.phi_c_bulk: float, default np.nan
        Bulk volume fraction of crowders

    Class Attributes
    ----------------
    _geometries: list of str
        Possible geometries of a simulation box
    _groups: list of str
        Possible groups of the `SumRule` project.
    _lineage_attributes: dict of dict
        a dictionary of `lineage` names. For each `lineage`, a dictionary
        maps the keywords of physical attributes in that lineage to their
        corresponding attributes in `SumRule` class.
    _lineage_private_attributes: dict of lists
        a dictionary of `lineage` names. For each `lineage`, a list of
        class attributes that are NOT "N/A" or np.nan is created that can be
        used in the simulation/experiment/run report files (called
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
        self._parse_lineage_name()
        self._set_parents()
        if self.lineage in ['segment', 'whole', 'ensemble_long']:
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
            Lineage: '{self.lineage}'
        """
        return observation

    def __repr__(self) -> str:
        return f"Observation('{self.filename}' in geometry '{self.geometry}' \
            from group '{self.group}' with lineage '{self.lineage}'"

    def _find_lineage_name(self):
        """
        parses the unique lineage_name (the first substring of filename
        and/or the segment keyword middle substring) of a filename.
        """
        if self.lineage in ['segment', 'whole']:
            # a 'segment' lineage only used in 'probe' phase
            # a 'whole' lineage used in 'probe' or 'analyze' phases
            # so its lineage_name is either ended by 'group' keyword or "-".
            # these two combined below:
            self.lineage_name = \
                self.filename.split("." + self.group)[0].split("-")[0]
        else:  # 'ensemble' or 'space' lineages
            self.lineage_name = self.filename.split('-')[0]

    def _initiate_attributes(self):
        """
        defines and initiates the class attribute based on the physical
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

    def _parse_lineage_name(self):
        """
        parses a lineage_name based on a list of keywords of physical
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
                    attr_value = 2 * attr_value - self.dmon
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


        The following map is used for setting relationships:
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
        Size (or diameter) of the `biaxial` or `slit` confinement, inferred
        from either 'r' keyword (radius of the biaxial confinement; a cylinder
        with open ends) or 'D' keyword (size of that confinement. Following
        LAMMPS' tango, `dcyl` ranged from '[-dcyl/2,dcyl.2] inculsive is the
        domain over which x and y cartesian coordiantes are defined in the
        'biaxial' geometry; however, `dcyl` ranged from '[-dcyl/2,dcyl.2]
        inculsive is the domain over which z cartesian coordiante is defined
        in the 'slit' geometry. It is important to note that `dcyl` is
        different from the size defined in LAMMPS input file if the
        wall-forming particle are defined; in this case:
            `self.dcyl` = LAMMPS.dcyl - `self.dwall`
        Hence, `dcyl` is defined differenty in different parser classes in
        this module.

    Class Attributes
    ----------------

     """
    _physical_attributes = {
        'deGennes': {
            'dimension': 'dim',
            'w3body': 'w',
            'nmon': 'n',
            'dcyl': 'd',
            'dcrowd': 'ac'
        },
        'twoBodyOnly': {
            'dimension': 'dim',
            'w3body': 'w',
            'nmon': 'n',
            'dcyl': 'd',
            'dcrowd': 'ac'
        },
        'shendruk': {
            'dimension': 'dim',
            'nmon': 'n',
            'dcyl': 'd',
            'dcrowd': 'ac'
        }
    }
    _column_names = {
        'deGennes': ['phi_c', 'vexc', 'swollen'],
        'twoBodyOnly': ['phi_c', 'vexc', 'swollen'],
        'shendruk': ['phi_c', 'vexc', 'w3body', 'swollen']
    }

    def __init__(
        self,
        name: str,
        limit: bool = False,
        ispath: bool = True
    ):
        self.limit = limit
        self.ispath = ispath
        if self.ispath:
            self.filepath = name
            self.extention = os.path.splitext(name.split("/")[-1])[1]
            self.filename = os.path.splitext(name.split("/")[-1])[0]
        else:
            self.filepath = "N/A"
            self.extention = "N/A"
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
        self.dmon = 1.0
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
        words = self.filename.split('-')
        self.w3body_model = words[1]
        self.vexc_model = words[2]
        attrs_pattern = re.compile(r'([a-zA-Z\-]+)')
        attrs = attrs_pattern.split(words[-1])
        attributes_float = ['dmon', 'dcyl', 'dcrowd', 'three_body']
        for attr_name, attr_keyword in \
                self._physical_attributes[self.w3body_model].items():
            try:
                attr_value = attrs[attrs.index(attr_keyword)+1]
                if attr_name in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                setattr(self, attr_name, attr_value)
            except ValueError:
                print(
                    f"'{attr_keyword}'"
                    " attribute keyword is not in "
                    f"'{self.filename}'"
                    " lineage name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not.")

    def _athermal_virial_coeffs(self):
        """
        compute the excluded volume (second viriral coefficient) and
        three-body coefficient (thrid virial coefficient) for the
        equivalent athermal system.
        """
        if self.vexc_model == 'AOExcVol':
            # The hard-sphere potenial's virial coefficient
            self.vexc_athr = (4/3) * np.pi * self.dmon**3
            self.w3body_athr = (5/8) * self.vexc_athr**2
        elif self.vexc_model == 'HaExcVol':
            # The fully-repulsive Weeks–Chandler–Anderson (WCA) potential
            # with r_cut = 2**(1/6) * sigma where sigma=self.dmon=1.0
            self.vexc_athr = 4.40944631
            self.w3body_athr = (5/8) * self.vexc_athr**2
        else:
            vexc_models = ['AOExcVol', 'HaExcVol']
            raise ValueError(
                f"'{self.vexc_model}' "
                "is not a valid model of the excluded volume of monomers."
                f"Please select one of these '{vexc_models}' models."
            )

    def _initial_chain_size(self):
        """
        compute the equlibirium chain size in the given dimension in
        the absent of crowders.
        """
        self.flory_exponent = 3 / (self.dimension + 2)
        if self.dimension == 3:
            self.r_flory = 1.12 * self.dmon * self.nmon ** self.flory_exponent
        elif self.dimension == 2:
            print("In 2-dimensional space, the pre-factor is set to 1.0.")
            self.r_flory = 1.0 * self.dmon * (
                self.nmon ** self.flory_exponent *
                (self.dmon/self.dcyl) ** (1/4)
                )
        elif self.dimension == 1:
            self.r_flory = 1.37 * self.dmon * (
                self.nmon ** self.flory_exponent *
                (self.dmon/self.dcyl) ** (2/3)
                )
        else:
            raise ValueError(
                f"'{self.dimension}' "
                "is not a valid dimension for the system."
                f"Please select one of these '{list(range(1,4,1))}' models."
            )
        self.r_ideal = self.dmon * self.nmon ** 0.5

    def _clean(self):
        """
        clean a dataset of a chain's size based on its attributes.
        """
        if self.ispath:
            df = pd.read_csv(
                self.filepath,
                index_col=False,
                names=self._column_names[self.w3body_model]
            )
        else:
            raise ValueError(
                f"'{self.filepath}' "
                "is not a valid filepath."
            )
        df.dropna(inplace=True)
        df.sort_values(by=['phi_c'], inplace=True)
        df.reset_index(inplace=True, drop=True)
        if self.w3body_model == 'deGennes':
            if self.nmon == 80:
                df = df.groupby('phi_c').min()
                df.reset_index(inplace=True)
        elif self.w3body_model == 'shendruk':
            df = df.groupby('phi_c').min()
            df.reset_index(inplace=True)
        else:
            df = df.sort_values(by=['swollen'])
        self.chain_size = df

    def _expand(self):
        """
        expand `self.chain_size` by adding the physical attributes to it as
         new columns.
        """
        self.chain_size['w3body_model'] = self.w3body_model
        self.chain_size['vexc_model'] = self.vexc_model
        for attr_name in self._physical_attributes[self.w3body_model].keys():
            self.chain_size[attr_name] = getattr(self, attr_name)

    def _scale(self, limit=False):
        """
        scale a chain's attributes.
        """
        self.chain_size['phi_c_scaled'] = \
            self.dmon * self.chain_size['phi_c'] / self.dcrowd
        # phi_c*a/a_c for a_c << a for maximum free energy gain:
        self.chain_size['vexc_scaled'] = \
            self.chain_size['vexc'] / self.vexc_athr
        self.chain_size['r'] = self.chain_size['swollen'] * self.r_ideal
        self.r_max = self.chain_size['r'].max()
        self.chain_size['r_scaled_max'] = self.chain_size['r'] / self.r_max
        self.chain_size['r_scaled_flory'] = self.chain_size['r'] / self.r_flory
        self.chain_size['w3body_scaled'] = \
            self.chain_size['w3body'] / self.w3body_athr
        if self.w3body_model == 'shendruk':
            self.w3body = (
                self.chain_size['w3body'].min(),
                self.chain_size['w3body'].max()
                )
        if limit:
            # Limit the data only to -1*exc_vol_athr <= exc_vol <= exc_vol_athr
            limit_range = \
                (self.chain_size['vexc_scaled'] <= 1.5) & (
                 self.chain_size['vexc_scaled'] >= -1.0)
            self.chain_size = self.chain_size[limit_range]


class DataTemplate(object):
    name: str
    ispath: bool = True
    """Parse the filename of a LAMMPS template data file based on one of the
    following file name templates:

    1. 'data_template-geometry-group-N#D#L#'

        where 'data_template' is the fixed header name of the LAMMPS template
        files, 'geometry' is the geoemtery of the space, 'group' is the name of
        file group, 'N" is the number of monomers, 'D' is the diameter (size)
        of the cylinder, and 'L' is the length of cylinder. This template is
        for 'biaxial' and 'slit' geomteries. It is important to note that the
        system is peroidic along z direction with range [-L/2,L/2] inclusive
        in the 'biaxial' geoemtry while it is periodic along x and y directions
        both with the same range [-D/2,D/2] in the 'slit' geometry. Here, x, y,
        and z are cartesian coordinates.

     2. 'data_template-geometry-group-N#L#'

        where 'data_template' is the fixed header name of the LAMMPS template
        files, 'geometry' is the geoemtery of the space, 'group' is the name of
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
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    nmon: int, np.nan
        Number of monomers. Its associated keyword is 'N'.
    dcyl: float, np.nan
        Size (or diameter) of the biaxial (cylindrical) confinement, inferred
        from 'D' keyword keyword. Cuation: The size of cylindrical confinement
        is defined based on the LAMMPS' tango, so it is different from the size
        defined in other parsers defined in this module.

    Class Attributes
    ----------------
    _groups: list of str
        Possible groups of the `SumRule` project.
    _geometries_attributes: dict of dict
        a dictionary of `geometry` names. For each `geometry`, a dictionary
        maps the keywords of physical attributes in that geometry to their
        corresponding attributes in `DataTemplate` class.
    """
    _groups = ['bug', 'all']
    _geometries_attributes = {
        'biaxial': {
            'nmon': 'N', 'dcyl': 'D', 'lcyl': 'L'
        },
        'slit': {
            'nmon': 'N', 'dcyl': 'D', 'lcyl': 'L'
        },
        'box': {
            'nmon': 'N',  'lcyl': 'L',
        }
    }

    def __init__(
        self,
        name: str,
        ispath: bool = True
    ):
        if ispath:
            self.filepath = name
            self.filename = name.split("/")[-1]
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
        self.geometry = 'N/A'
        self.group = 'N/A'
        self.nmon = np.nan
        self.dcyl = np.nan
        self.lcyl = np.nan

    def _parse_filename(self):
        """
        parses a lineage_name based on a list of keywords of physical
        attributes.
        """
        words = self.filename.split('-')
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
            geometries_string = "'" + "', '".join(
                self._geometries_attributes.keys()) + "'"
            raise ValueError(
                f"'{geometry}' "
                "is not a valid geometry. Please select one of "
                f"{geometries_string} geometries.")
        group = words[2]
        if group in self._groups:
            self.group = group
        else:
            groups_string = "'" + "', '".join(
                self._groups) + "'"
            raise ValueError(
                f"'{group}' "
                "is not a valid particle group. Please select one of "
                f"{groups_string} groups.")
        str_attrs = re.compile(r'([a-zA-Z\-]+)')
        attributes_float = ['dcyl', 'lcyl']
        attrs = str_attrs.split(words[3])
        for attr_name, attr_keyword in \
                self._geometries_attributes[self.geometry].items():
            try:
                attr_value = attrs[attrs.index(attr_keyword)+1]
                if attr_name in attributes_float:
                    attr_value = float(attr_value)
                else:
                    attr_value = int(float(attr_value))
                setattr(self, attr_name, attr_value)
            except ValueError:
                print(
                    f"'{attr_keyword}'"
                    " attribute keyword is not in "
                    f"'{self.filename}'"
                    " template data name. Please check whether "
                    f"'{self.filename}'"
                    " is valid name or not.")
