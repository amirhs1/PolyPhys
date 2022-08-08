from glob import glob
from typing import (
    List,
    Union,
    Callable,
    Optional
)
import numpy as np
import pandas as pd
import itertools
from polyphys.analyze import analyzer
from polyphys.manage.utilizer import round_up_nearest


def allInOne_equil_tseries(
    project: str,
    analysis_db: str,
    group: str,
    spaces: List[str],
    properties: List[str],
    measures: List[Callable],
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    save_space: Optional[bool] = False,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """Performs a group of `measures` on a collection of physical `properties`
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
        The names of all the spaces in a prorject.
    properties: list of str
        The names of physical properties.
    measures: list of Callable
        The list of applying measures/functions.
    divisor: float, default 0.025
        The step by which the values of "phi_c_bulk" attribute are rounded.
    round_to: int, default 3
        The number of significant decimal digits in the values of "phi_c_bulk"
        attribute.
    save_space: bool, default False
        Whether save the "space" dataframes or not.
    save_to : str, default None
        Absolute or relative path to which the output is wrriten.

    Requirements
    ------------
    polyphys, Pandas.
    """
    if save_space:
        save_to_space = save_to
    else:
        save_to_space = False
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
        spac_equil_props = analyzer.equilibrium_tseries_wholes(
            space,
            space_db + "/*.csv",
            properties,
            measures,
            whole_stamps,
            save_to=save_to_space
        )
        all_in_one_equil_props.append(spac_equil_props)
    all_in_one_equil_props = pd.concat(all_in_one_equil_props)
    # add rounded phi_crds to the dataset
    all_in_one_equil_props['phi_c_bulk_round'] = \
        all_in_one_equil_props['phi_c_bulk'].apply(
            round_up_nearest,
            args=[divisor, round_to]
        )
    output = "-".join(
        ["allInOne", project, group, "equilProps-whole.csv"]
    )
    all_in_one_equil_props.to_csv(save_to + output, index=False)
    return all_in_one_equil_props


def allInOne_equil_tseries_ensAvg(
    project: str,
    project_db: Union[pd.DataFrame, str],
    group: str,
    properties: List[str],
    attributes: List[str],
    save_to: Optional[str] = None
):
    """Perform ensemble-avergaing and then normalization on the equilibrium
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
        Absolute or relative path to which the output is wrriten.

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
    for prop in norm_props:
        ens_avg[prop + "-norm"] = 0
    for space, prop in itertools.product(spaces, norm_props):
        space_con = ens_avg['space'] == space
        phi_c_con = ens_avg['phi_c_bulk_round'] == 0
        prop_0 = ens_avg.loc[space_con & phi_c_con, prop + "-mean"].values[0]
        ens_avg.loc[space_con, prop + "-norm"] = \
            ens_avg.loc[space_con, prop + "-mean"] / prop_0
    if save_to is not None:
        output = "-".join(
            ["allInOne", project, group, "equilProps-ensAvg.csv"]
        )
        ens_avg.to_csv(save_to + output, index=False)
    return ens_avg
