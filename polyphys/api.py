from glob import glob
from typing import (
    List,
    Union,
    Callable,
    Optional
)
import warnings
import numpy as np
import itertools
import pandas as pd
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
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'fsdMon-mean', 'fsdMon-var', 'fsdMon-sem',
            'gyrMon-mean', 'gyrMon-var', 'gyrMon-sem',
            'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem',
            'rfloryMon-var', 'rfloryMon-sem', 'shapeMon-mean',
            'transSizeMon-mean', 'transSizeMon-var', 'transSizeMon-sem'
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
        'attributes': ['space', 'ensemble_long', 'ensemble',
                       'nmon_small', 'nmon_large', 'dmon_large',
                       'dcyl', 'dcrowd', 'phi_c_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'shapeTMon'
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space',
                             'dcyl', 'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'
                             ],
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'fsdMon-mean', 'fsdMon-var', 'fsdMon-sem',
            'gyrMon-mean', 'gyrMon-var', 'gyrMon-sem',
            'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem',
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
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space',
                             'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'
                             ],
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'gyrMon-mean', 'gyrMon-var', 'gyrMon-sem',
            'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem'
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
        'space_pat': 'N*kbmm*nh*ac*epshc*',
        'hierarchy': 'N*',
        'directions': ['r'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble',
                       'bend_mm', 'eps_hc', 'nmon', 'nhns', 'dcrowd',
                       'phi_c_bulk', 'rho_hns_bulk',
                       ],
        'time_varying_props': ['asphericityTMon', 'gyrTMon',
                               'shapeTMon', 'bondLengthVecMon',
                               'nBoundTHnsPatch', 'nFreeTHnsPatch',
                               'nEngagedTHnsPatch', 'nFreeTHnsCore',
                               'nBridgeTHnsCore', 'nDangleTHnsCore',
                               'nCisTHnsCore', 'nTransTHnsCore',
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'eps_hc',
                             'bend_mm', 'nmon', 'nhns', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round', 'rho_hns_bulk'],
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'gyrMon-mean', 'gyrMon-var', 'gyrMon-sem',
            'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem',
            'bondLengthMon-mean', 'bondLengthMon-var', 'bondLengthMon-sem',
            'nBoundHnsPatch-mean', 'nBoundHnsPatch-var', 'nBoundHnsPatch-sem',
            'nFreeHnsPatch-mean', 'nFreeHnsPatch-var', 'nFreeHnsPatch-sem',
            'nEngagedHnsPatch-mean', 'nEngagedHnsPatch-var',
            'nEngagedHnsPatch-sem',
            'nFreeHnsCore-mean', 'nFreeHnsCore-var', 'nFreeHnsCore-sem',
            'nBridgeHnsCore-mean', 'nBridgeHnsCore-var', 'nBridgeHnsCore-sem',
            'nDangleHnsCore-mean', 'nDangleHnsCore-var', 'nDangleHnsCore-sem',
            'nCisHnsCore-mean', 'nCisHnsCore-var', 'nCisHnsCore-sem',
            'nTransHnsCore-mean', 'nTransHnsCore-var', 'nTransHnsCore-sem'
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
        'space_pat': 'N*D*nh*ac*epshc*',
        'hierarchy': 'N*',
        'directions': ['r', 'z'],
        'props': ['Rho', 'Phi'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'dcyl',
                       'nmon', 'nhns', 'dcrowd', 'phi_c_bulk', 'rho_hns_bulk'
                       ],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'shapeTMon',
                               'gyrTMon', 'transSizeTMon',
                               'bondLengthVecMon',
                               'nBoundTHnsPatch', 'nFreeTHnsPatch',
                               'nEngagedTHnsPatch', 'nFreeTHnsCore',
                               'nBridgeTHnsCore', 'nDangleTHnsCore',
                               'nCisTHnsCore', 'nTransTHnsCore',
                               ],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'dcyl',
                             'nhns', 'nmon', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round', 'rho_hns_bulk'],
        'equil_properties': [
            'asphericityMon-mean', 'asphericityMon-var', 'asphericityMon-sem',
            'fsdMon-mean', 'fsdMon-var', 'fsdMon-sem',
            'gyrMon-mean', 'gyrMon-var', 'gyrMon-sem',
            'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem',
            'transSizeMon-mean', 'transSizeMon-var', 'transSizeMon-sem',
            'bondLengthMon-mean', 'bondLengthMon-var', 'bondLengthMon-sem',
            'nBoundHnsPatch-mean', 'nBoundHnsPatch-var', 'nBoundHnsPatch-sem',
            'nFreeHnsPatch-mean', 'nFreeHnsPatch-var', 'nFreeHnsPatch-sem',
            'nEngagedHnsPatch-mean', 'nEngagedHnsPatch-var',
            'nEngagedHnsPatch-sem',
            'nFreeHnsCore-mean', 'nFreeHnsCore-var', 'nFreeHnsCore-sem',
            'nBridgeHnsCore-mean', 'nBridgeHnsCore-var', 'nBridgeHnsCore-sem',
            'nDangleHnsCore-mean', 'nDangleHnsCore-var', 'nDangleHnsCore-sem',
            'nCisHnsCore-mean', 'nCisHnsCore-var', 'nCisHnsCore-sem',
            'nTransHnsCore-mean', 'nTransHnsCore-var', 'nTransHnsCore-sem'
        ],
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


def load_project_db(
    project: str,
    project_db: Union[pd.DataFrame, str]
) -> pd.DataFrame:
    """Load the project database.

    Parameters
    ----------
    project : str
        The name of the project.
    project_db : str or pd.DataFrame
        The project "whole" "allInOne" dataframe or the path to it.

    Returns
    -------
    pd.DataFrame
        The project database.
    """
    if isinstance(project_db, str):
        project_path = "-".join(["allInOne", project, "equilProps-whole.csv"])
        project_db = pd.read_csv(project_db + project_path, header=0)

    return project_db


def drop_unnecessary_columns(
    df: pd.DataFrame,
    properties: List[str],
    attributes: List[str]
) -> pd.DataFrame:
    """Drop columns from the dataframe that are not properties or attributes.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to process.
    properties : list of str
        The names of physical properties.
    attributes : list of str
        The name of physical attributes.

    Returns
    -------
    pd.DataFrame
        The processed dataframe.
    """
    cols_to_drop = list(
        set(df.columns).difference(set(attributes + properties)))
    df.drop(columns=cols_to_drop, inplace=True)

    return df


def normalize_data(
    df: pd.DataFrame,
    project: str,
    properties: List[str]
) -> pd.DataFrame:
    """Normalize the data in the dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to process.
    project : str
        The name of the project.
    properties : list of str
        The names of physical properties.

    Returns
    -------
    pd.DataFrame
        The processed dataframe.
    """ 
    df_copy = df.copy()
    norm_props = [
        prop.split('-')[0] for prop in properties if prop.endswith('mean')]
    spaces = df_copy['space'].unique()
    if project in ['HnsCyl', 'HnsCub']:
        for space, prop in itertools.product(spaces, norm_props):
            space_con = df_copy['space'] == space
            phi_c_con = df_copy['phi_c_bulk_round'] == 0
            prop_0 = \
                df_copy.loc[
                    space_con & phi_c_con, prop + "-mean"
                    ].to_numpy()[0]  # type: ignore
            if prop_0 != 0:
                df_copy.loc[space_con, prop + "-norm"] = \
                    df_copy.loc[space_con, prop + "-mean"] / prop_0
            else:
                warnings.warn(
                    f"The '{prop}' value in the absence of crowders "
                    f"(phi_c=0) for space '{space}, the values of '{prop}-norm' "
                    "at all the values of phi_c are set to 'np.nan'.",
                    UserWarning
                )
                df_copy.loc[space_con, prop + "-norm"] = np.nan
        n_patch_per_cor = 2
        hpatch_cols = ['nBoundHnsPatch', 'nFreeHnsPatch', 'nEngagedHnsPatch']
        hcore_cols = ['nFreeHnsCore', 'nBridgeHnsCore', 'nDangleHnsCore',
                      'nCisHnsCore', 'nTransHnsCore']
        for col in hpatch_cols:
            df_copy[col + '-norm'] = \
                df_copy[col + '-mean'] / (df_copy['nhns'] * n_patch_per_cor)
        for col in hcore_cols:
            df_copy[col + '-norm'] = df_copy[col + '-mean'] / df_copy['nhns']
    else:
        for space, prop in itertools.product(spaces, norm_props):
            space_con = df_copy['space'] == space
            phi_c_con = df_copy['phi_c_bulk_round'] == 0
            prop_0 = \
                df_copy.loc[
                    space_con & phi_c_con, prop + "-mean"
                    ].to_numpy()[0]  # type: ignore
            if prop_0 != 0:
                df_copy.loc[space_con, prop + "-norm"] = \
                    df_copy.loc[space_con, prop + "-mean"] / prop_0
            else:
                warnings.warn(
                    f"The '{prop}' value in the absence of crowders "
                    f"(phi_c=0) for space '{space}, the values of '{prop}-norm' "
                    "at all the values of phi_c are set to 'np.nan'.",
                    UserWarning
                )
                df_copy.loc[space_con, prop + "-norm"] = np.nan
    return df_copy


def all_in_one_equil_tseries_ens_avg(
    project: str,
    project_db: Union[pd.DataFrame, str],
    group: str,
    properties: List[str],
    attributes: List[str],
    save_to: Optional[str] = None
) -> pd.DataFrame:
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

    Returns
    -------
    pd.DataFrame
        The dataframe of ensemble-averaged properties.

    Raises
    ------
    ValueError
        If `properties` or `attributes` is empty.
    KeyError
        If required columns are not found in the `project_db` dataframe.

    Requirements
    ------------
    Numpy, Pandas.
    """
    if not properties:
        raise ValueError("The `properties` list is empty.")
    if not attributes:
        raise ValueError("The `attributes` list is empty.")

    project_db = load_project_db(project, project_db)

    required_cols = attributes + properties
    missing_cols = set(required_cols) - set(project_db.columns)
    if missing_cols:
        raise KeyError(
            f"Missing required columns in the `project_db` dataframe: "
            f"{missing_cols}"
            )

    ens_avg = project_db.copy()
    ens_avg = drop_unnecessary_columns(ens_avg, properties, attributes)
    # Ensemble-averaging all the measures of all the properties:
    ens_avg = ens_avg.groupby(attributes).agg(np.mean).reset_index()
    # Normalizing the mean values of each property in each ensemble in a
    # space by the value of the property ensemble with phi_c_bulk_round=0 in
    # that space:
    # Here, the normalization is only performed for the "mean" measure not, the
    # "std" or "sem" measures.
    ens_avg = normalize_data(ens_avg, project, properties)

    if save_to is not None:
        output = "-".join(
            ["allInOne", project, group, "equilProps-ensAvg.csv"]
            )
        ens_avg.to_csv(save_to + output, index=False)

    return ens_avg
