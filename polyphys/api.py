from glob import glob
from typing import (
    List,
    Union,
    Callable,
    Optional
)
import warnings
import numpy as np
import pandas as pd
import itertools
from polyphys.manage.parser import (SumRuleCyl, TransFociCub, TransFociCyl,
                                    HnsCub, HnsCyl)
from polyphys.analyze import analyzer
from polyphys.analyze import measurer
from polyphys.manage.utilizer import round_up_nearest

PROJECTS_DETAILS = {
    'SumRuleCyl': {
        'group': 'bug',
        'geometry': 'cylindrical',
        'topology': 'linear',
        'parser': SumRuleCyl,
        'divisor': 0.025,
        'space_pat': 'N*D*ac*',
        'hierarchy': 'N*',
        'directions': ['r', 'z'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon', 'dcyl',
                       'dcrowd', 'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'rfloryTMon', 'shapeTMon', 'transSizeTMon'
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['space', 'ensemble_long', 'ensemble', 'nmon',
                             'dcyl', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round'
                             ],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'fsdMon-mean', 'fsdMon-var',
                             'fsdMon-sem', 'gyrMon-mean', 'gyrMon-var',
                             'gyrMon-sem', 'rfloryMon-mean', 'rfloryMon-var',
                             'rfloryMon-sem', 'shapeMon-mean', 'shapeMon-var',
                             'shapeMon-sem', 'transSizeMon-mean',
                             'transSizeMon-var', 'transSizeMon-sem'
                             ],
        'rhosPhisNormalizedScaled': [('Mon', 'dmon'), ('Crd', 'dcrowd')]
    },
    'TransFociCyl': {
        'group': 'bug',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'parser': TransFociCyl,
        'divisor': 0.025,
        'space_pat': 'ns*nl*al*D*ac*',
        'hierarchy': 'eps*',
        'directions': ['r', 'z'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'ns*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon_small',
                       'nmon_large', 'dmon_large', 'dcyl', 'dcrowd',
                       'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'shapeTMon'
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'dcyl',
                             'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'
                             ],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'fsdMon-mean',
                             'fsdMon-var', 'fsdMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem'
                             ],
        'rhosPhisNormalizedScaled': [('Mon', 'dmon_small'), ('Crd', 'dcrowd'),
                                     ('Foci', 'dmon_large')
                                     ]
    },
    'TransFociCub': {
        'group': 'bug',
        'geometry': 'cubic',
        'topology': 'ring',
        'parser': TransFociCub,
        'divisor': 0.025,
        'space_pat': 'ns*nl*al*ac*',
        'hierarchy': 'al*',
        'directions': ['r'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'ns*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon_small',
                       'nmon_large', 'dmon_large', 'dcrowd', 'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'gyrTMon', 'shapeTMon'],
        'equil_measures': [np.mean],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space',
                             'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'
                             ],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem'
                             ],
        'rhosPhisNormalizedScaled': [('Mon', 'dmon_small'), ('Crd', 'dcrowd'),
                                     ('Foci', 'dmon_large')
                                     ]
    },
    'HnsCub': {
        'group': 'nucleoid',
        'geometry': 'cubic',
        'topology': 'ring',
        'parser': HnsCub,
        'divisor': 0.01,
        'space_pat': 'N*epshm*nh*ac*',
        'hierarchy': 'N*',
        'directions': ['r'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'eps_hm',
                       'nmon', 'nhns', 'dcrowd', 'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'gyrTMon', 'shapeTMon'],
        'equil_measures': [np.mean],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space',
                             'eps_hm', 'nmon', 'nhns', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round'
                             ],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem'
                             ],
        'rhosPhisNormalizedScaled': [('Mon', 'dmon'), ('Crd', 'dcrowd'),
                                     ('Hns', 'dhns')]
    },
    'HnsCyl': {
        'group': 'nucleoid',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'parser': HnsCyl,
        'divisor': 0.01,
        'space_pat': 'N*D*nh*ac*',
        'hierarchy': 'N*',
        'directions': ['r', 'z'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'dcyl',
                       'nmon', 'nhns', 'dcrowd', 'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'shapeTMon', 'nBoundHnsPatchT',
                               'nFreeHnsPatchT', 'nEngagedHnsPatchT',
                               'nFreeHnsCoreT', 'nBridgeHnsCoreT',
                               'nDangleHnsCoreT', 'nBoundMonT',
                               'nCisMonT', 'nTransMonT', 'bondLengthVecMon'
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'dcyl',
                             'nhns', 'nmon', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round'],
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'fsdMon-mean', 'fsdMon-var', 'fsdMon-sem', 'gyrMon-mean',
            'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean', 'shapeMon-var',
            'shapeMon-sem', 'nBoundHnsPatch-mean', 'nBoundHnsPatch-var',
            'nBoundHnsPatch-sem', 'nFreeHnsPatch-mean', 'nFreeHnsPatch-var',
            'nFreeHnsPatch-sem', 'nEngagedHnsPatch-mean',
            'nEngagedHnsPatch-var', 'nEngagedHnsPatch-sem',
            'nFreeHnsCore-mean', 'nFreeHnsCore-var', 'nFreeHnsCore-sem',
            'nBridgeHnsCore-mean', 'nBridgeHnsCore-var', 'nBridgeHnsCore-sem',
            'nDangleHnsCore-mean', 'nDangleHnsCore-var', 'nDangleHnsCore-sem',
            'nBoundMon-mean', 'nBoundMon-var', 'nBoundMon-sem', 'nCisMon-mean',
            'nCisMon-var', 'nCisMon-sem', 'nTransMon-mean', 'nTransMon-var',
            'nTransMon-sem', 'bondLengthMon-mean', 'bondLengthMon-var',
            'bondLengthMon-sem'],
        'rhosPhisNormalizedScaled': [('Mon', 'dmon'), ('Crd', 'dcrowd'),
                                     ('Hns', 'dhns')]
    }
}


def all_in_one_equil_tseries(
    project: str,
    analysis_db: str,
    group: str,
    spaces: List[str],
    properties: List[str],
    measures: List[Callable],
    kind: str = 'dataframe',
    topology: Optional[str] = None,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    save_space: Optional[bool] = False,
    save_to: Optional[str] = None,
) -> pd.DataFrame:
    """
    Performs a group of `measures` on a collection of physical `properties`
    of a `species` in a `group` in all the `spaces` in the "ens" phase of the
    "analysis" stage of a `project` and merges the resulting dataframes as an
    "allInOne" dataframe.

    This function adds "phi_c_bulk_round" as anew column to the final dataset.

    Each statistical measure is applied to each "whole" *times series* (a
    column in an "ensemble" data frame) in each "ensemble" of a given physical
    property in a space.

    Parameters
    ----------
    project: str
        The name of a project.
    analysis_db: str
        The  path to the "analysis" directory of the project.
    group: str
        The name of the group to which an `species` belong.
    spaces: str
        The names of all the spaces in a project.
    properties: list of str
        The names of physical properties.
    measures: list of Callable
        The list of applying measures/functions.
    kind: {'dataframe', 'array'}, default 'dataframe'
        Type of 'properties' file.
    topology : str, optional
        The topology of the polymer, if any.
    divisor: float, default 0.025
        The step by which the values of "phi_c_bulk" attribute are rounded.
    round_to: int, default 3
        The number of significant decimal digits in the values of "phi_c_bulk"
        attribute.
    save_space: bool, default False
        Whether save the "space" dataframes or not.
    save_to : str, default None
        Absolute or relative path to which the output is written.
    """
    if save_space:
        save_to_space = save_to
    else:
        save_to_space = None
    all_in_one_equil_props = []
    for space in spaces:
        space_db = analysis_db + "-".join([space, group, "ens"])
        whole_stamps = space_db + "/*stamps*.csv"
        whole_stamps = glob(whole_stamps)
        if len(whole_stamps) > 1:
            raise ValueError(
                "More than one 'whole' stamps dataset found. Which of the"
                f" following is the correct one? '{whole_stamps}'")
        whole_stamps = pd.read_csv(whole_stamps[0], header=0)
        space_equil_props = analyzer.equilibrium_wholes(
            space,
            space_db,
            properties,
            measures,
            whole_stamps,
            output_type=kind,
            topology=topology,
            output_path=save_to_space
        )
        all_in_one_equil_props.append(space_equil_props)
    all_in_one_equil_props = pd.concat(all_in_one_equil_props)
    # add rounded phi_crds to the dataset
    all_in_one_equil_props['phi_c_bulk_round'] = \
        all_in_one_equil_props['phi_c_bulk'].apply(
            round_up_nearest,
            args=[divisor, round_to]
        )
    if save_to is not None:
        output = "-".join(["allInOne", project, group, "equilProps-whole.csv"])
        all_in_one_equil_props.to_csv(save_to + output, index=False)
    return all_in_one_equil_props


def all_in_one_equil_tseries_ens_avg(
    project: str,
    project_db: Union[pd.DataFrame, str],
    group: str,
    properties: List[str],
    attributes: List[str],
    save_to: Optional[str] = None
):
    """Perform ensemble-averaging and then normalization on the equilibrium
    properties in a `project`.

    Parameters
    ----------
    project: str
        The name of a project.
    project_db: str or pd.Dataframe
        The project "whole" "allInOne" dataframe or the path to it.
    group: str
        The name of the group to which an `species` belong.
    properties: list of str
        The names of physical properties.
    attributes: list of str
        The name of physical attributes.
    save_to : str, default None
        Absolute or relative path to which the output is written.

    Return
    ------
    ens_avg: pd.DataFrame
        The dataframe of ensemble-averaged properties.

    Requirements
    ------------
    Numpy, Pandas.
    """
    if isinstance(project_db, str):
        project_path = "-".join(
            ["allInOne", project, group, "equilProps-whole.csv"]
        )
        project_db = pd.read_csv(project_db + project_path, header=0)
    cols_to_drop = list(
        set(project_db.columns).difference(set(attributes + properties))
    )
    project_db.drop(columns=cols_to_drop, inplace=True)
    # Ensemble-averaging all the measures of all the properties:
    ens_avg = project_db.groupby(attributes).agg(np.mean)
    ens_avg.reset_index(inplace=True)
    # Normalizing the mean values of each property in each ensemble in a
    # space by the value of the property ensemble with phi_c_bulk_round=0 in
    # that space:
    # Here, the normalization is only performed for the "mean" measure not, the
    # "std" or "sem" measures.
    spaces = ens_avg['space'].unique()
    norm_props = [
        prop.split('-')[0] for prop in properties if prop.endswith('mean')
    ]
    for space, prop in itertools.product(spaces, norm_props):
        space_con = ens_avg['space'] == space
        phi_c_con = ens_avg['phi_c_bulk_round'] == 0
        prop_0 = ens_avg.loc[space_con & phi_c_con, prop + "-mean"].values[0]
        if prop_0 != 0:
            ens_avg.loc[space_con, prop + "-norm"] = \
                ens_avg.loc[space_con, prop + "-mean"] / prop_0
        else:
            warnings.warn(
                f"The '{prop}' value in the absence of crowders "
                f"(phi_c=0) for space '{space}, the values of '{prop}-norm' "
                "at all the values of phi_c are set to 'np.nan'.",
                UserWarning
                )
            ens_avg.loc[space_con, prop + "-norm"] = np.nan
    if save_to is not None:
        output = "-".join(
            ["allInOne", project, group, "equilProps-ensAvg.csv"]
        )
        ens_avg.to_csv(save_to + output, index=False)
    return ens_avg
