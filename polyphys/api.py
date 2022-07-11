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
from polyphys.manage import organizer
from polyphys.analyze import analyzer
from polyphys.manage.parser import SumRule
from polyphys.manage.utilizer import round_up_nearest


def allInOne_equil_tseries(
    project: str,
    analysis_db: str,
    species: str,
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
    species: str
        The name of species.
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
        ["allInOne", project, group, species, "equilProps-whole.csv"]
    )
    all_in_one_equil_props.to_csv(save_to + output, index=False)
    return all_in_one_equil_props


def allInOne_equil_tseries_ensAvg(
    project: str,
    project_db: Union[pd.DataFrame, str],
    species: str,
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
    species: str
        The name of species.
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
            ["allInOne", project, group, species, "equilProps-whole.csv"]
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
            ["allInOne", project, group, species, "equilProps-ensAvg.csv"]
        )
        ens_avg.to_csv(save_to + output, index=False)
    return ens_avg


def all_in_one_distributions(
                dist_tuples, properties, geometry, direction, save_to=None,
                round_to=4):
    """takes ensemble-averaged distributions and performs three operations
    on them. First, it concatenates the ensemble-averaged distributions of
    all the species in an ensebmle into one dataframe.

    In this dataframe, bin centers are indices and the number of columns
    are equal to the number of species. Next, it adds the properties of
    the ensemble to the merged distributions. Finally, it combines all
    the ensebles from all the gropus into one dataframe.

    Cautions:
    For each species, a distribution is a dataframe with bin centers as
    indices and a column of frequencies.
    A simulation group usually results in a graph or curve for the project
    and refers to a collection of simulations that all have the same values
    for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations
    that differs only in their random number seeds, initial conditions, or
    boundary conditions but have the same input parameters. In standard
    statitatical mechanical approach, an ensmeble is equivalent a simulation,
    but here we used it to reffer to all the thermodynamically-equivalent
    simulations.
    An ensemble-averaged group is an average over all the simulations in an
    ensemble and usually gives a data point.
    If there are N ensembles, each with M simulations, then there are N
    ensemble-average groups and N*M simulations in the simulation group.
    By default, the indexes of an ensemble (dataframe) are timesteps.
    For some properties such as histograms, however, the indexes are bin
    centers. As a result, simulations are tuples where each tuple can have
    one or more members. For histogram-like properties, each tuple has two
    members where the second member is the bin edges.

    Issue:
    1. Current implementation only work for a pair of species
    (e.g., a monomer type and a crowder type).
    2. This function only works when we have all the three distributions for
    a species: local histograms, densities, and volume fractions.

    Parameters:
    dist_tuples (list of tuples): A list of tuples in which each tuple has
    all the A*B ensemble-averaged distributions in one ensemble;
    here, A is the number of distribution types in this direction
    (histogram of particle numbers, number desnity, and volume fraction) ,
    and M is the number of different types of particles in that ensemble.
    properties (Pandas dataframe): The database of ensemble-average properties.
    geomtery (str): the geometry of the simulation box.
    direction (str): Direction along which the distributions are computed.
    save_to (bool): whether save to file or not
    round_to (int): rounding the whole dataframe to the given round-off number.

    Return;
    A pandad databse in which all the dsitbrutions of all the simulation
    groups are merged.
    """
    group_species = []
    properties[['ens_name', 'garbage']] = properties.filename.str.split(
        pat='-', expand=True)  # find the ensemble names.
    selected_cols = [
        'filename', 'nmon', 'dcyl', 'lcyl', 'phi_m_bulk', 'rho_m_bulk',
        'dcrowd', 'phi_c_bulk', 'rho_c_bulk', 'phi_c_bulk_normalized',
        'phi_c_bulk_eff', 'phi_c_bulk_eff_normalized', 'ens_name']
    for hist_crd, hist_mon, phi_crd, phi_mon, rho_crd, rho_mon in dist_tuples:
        hist_mon_df = pd.read_csv(hist_mon, index_col=0)
        ens_name = list(hist_mon_df.columns)[0].split('-')[0]
        hist_crd_df = pd.read_csv(
            hist_crd, names=['hist_crd_' + direction], skiprows=1, index_col=0)
        hist_mon_df = pd.read_csv(
            hist_mon, names=['hist_mon_' + direction], skiprows=1, index_col=0)
        rho_crd_df = pd.read_csv(
            rho_crd, names=['rho_crd_' + direction], skiprows=1, index_col=0)
        rho_mon_df = pd.read_csv(
            rho_mon, names=['rho_mon_' + direction], skiprows=1, index_col=0)
        phi_crd_df = pd.read_csv(
            phi_crd, names=['phi_crd_' + direction], skiprows=1, index_col=0)
        phi_mon_df = pd.read_csv(
            phi_mon, names=['phi_mon_' + direction], skiprows=1, index_col=0)
        ens_species = pd.concat([
            hist_mon_df, hist_crd_df, rho_mon_df, rho_crd_df, phi_mon_df,
            phi_crd_df], axis=1)
        ens_species.reset_index(inplace=True)
        ens_species.rename(columns={'index': direction}, inplace=True)
        # add the properties of the ensemble to merged distributions
        for col in selected_cols:
            cond = properties['ens_name'] == ens_name
            ens_species[col] = properties[cond][col].values[0]
        cell_attrs = SumRule(ens_name, geometry, warning=False)
        # normalizing the bin centers (or directions) in three different ways:
        bin_center_norm_box = {
            'r': cell_attrs.dcyl,
            'z': cell_attrs.lcyl
            }
        ens_species[direction+'_norm'] = \
            2 * ens_species[direction] / bin_center_norm_box[direction]
        ens_species[direction+'_norm_mon'] = ens_species[direction]  \
            # This should be divided by 'dmon' but this operation is not\
        # done since dmon = 1.0.
        ens_species[direction+'_norm_crd'] = \
            ens_species[direction] / cell_attrs.dcrowd
        # normalizing local volume fractions in different directions
        if direction == 'r':
            phi_mon_origin = ens_species.loc[0, 'phi_mon_r']
            phi_crd_infinity = ens_species.loc[0, 'phi_crd_r']  \
                # it is actually origin but for the sake of similiarity with\
            # z direction.
            rho_mon_origin = ens_species.loc[0, 'rho_mon_r']
            rho_crd_infinity = ens_species.loc[0, 'rho_crd_r']  \
                # it is actually origin but for the sake of similiarity with\
            # z direction
        if direction == 'z':
            around_max = 20
            non_zero_phi_mon_z_max_idx = ens_species[
                ens_species['phi_mon_z'] != 0]['phi_mon_z'].idxmax()
            phi_mon_origin = ens_species.loc[
                non_zero_phi_mon_z_max_idx - around_max:
                non_zero_phi_mon_z_max_idx + around_max, 'phi_mon_z'].mean()
            phi_crd_infinity = ens_species[
                ens_species['phi_mon_z'] == 0.0]['phi_crd_z'].mean()
            non_zero_rho_mon_z_max_idx = ens_species[
                ens_species['rho_mon_z'] != 0]['rho_mon_z'].idxmax()
            rho_mon_origin = ens_species.loc[
                non_zero_rho_mon_z_max_idx - around_max:
                non_zero_rho_mon_z_max_idx + around_max, 'rho_mon_z'].mean()
            rho_crd_infinity = ens_species[
                ens_species['rho_mon_z'] == 0.0]['rho_crd_z'].mean()
        ens_species['phi_mon_' + direction + '_uniform'] = phi_mon_origin
        ens_species['phi_mon_' + direction + '_uniform_size'] = \
            phi_mon_origin  \
            # This should be divided by 'dmon' but this operation is not\
        # done since dmon = 1.0
        ens_species['phi_mon_' + direction + '_uniform_size_norm'] = 1.0
        ens_species['phi_mon_' + direction + '_norm'] = \
            ens_species['phi_mon_' + direction] / phi_mon_origin
        ens_species['phi_mon_' + direction + '_norm_size'] = \
            ens_species['phi_mon_' + direction + '_norm']  \
            # This should be divided by 'dmon' but this operation is not\
        # done since dmon = 1.0
        ens_species['phi_mon_' + direction + '_size'] = \
            ens_species['phi_mon_' + direction]  \
            # This should be divided by 'dmon' but this operation is not\
        # done since dmon = 1.0
        ens_species['phi_crd_' + direction + '_uniform'] = phi_crd_infinity
        ens_species['phi_crd_' + direction + '_uniform_size'] = \
            phi_crd_infinity / cell_attrs.dcrowd
        ens_species['phi_crd_' + direction + '_uniform_size_norm'] = 1.0
        ens_species['phi_crd_' + direction + '_norm'] = \
            ens_species['phi_crd_' + direction] / phi_crd_infinity
        ens_species['phi_crd_' + direction + '_norm_size'] = \
            ens_species['phi_crd_' + direction + '_norm'] / cell_attrs.dcrowd
        ens_species['phi_crd_' + direction + '_size'] = \
            ens_species['phi_crd_' + direction] / cell_attrs.dcrowd
        ens_species['phi_' + direction + '_uniform_sum'] = \
            ens_species['phi_mon_' + direction + '_uniform_size']
        + ens_species['phi_crd_' + direction + '_uniform_size']
        ens_species['phi_' + direction + '_uniform_sum_norm'] = 1.0
        ens_species['rho_mon_' + direction + '_uniform'] = rho_mon_origin
        ens_species['rho_mon_' + direction + '_uniform_size'] = \
            rho_mon_origin  \
            # This should be mutiplied by 'dmon'**2 but this operation is\
        # not done since dmon = 1.0
        ens_species['rho_mon_' + direction + '_norm'] = \
            ens_species['rho_mon_' + direction] / rho_mon_origin
        ens_species['rho_mon_' + direction + '_norm_size'] = \
            ens_species['rho_mon_' + direction + '_norm'] \
            # This should be mutiplied by 'dmon'**2 but this operation is\
        # not done since dmon = 1.0
        ens_species['rho_mon_' + direction + '_size'] = \
            ens_species['rho_mon_' + direction] \
            # This should be mutiplied by 'dmon'**2 but this operation is\
        # not done since dmon = 1.0
        ens_species['rho_crd_' + direction + '_uniform'] = rho_crd_infinity
        ens_species['rho_crd_' + direction + '_uniform_size'] = \
            (rho_crd_infinity * cell_attrs.dcrowd**2)
        ens_species['rho_crd_' + direction + '_norm'] = \
            (ens_species['rho_crd_' + direction] / rho_crd_infinity)
        ens_species['rho_crd_' + direction + '_norm_size'] = \
            (ens_species['rho_crd_' + direction + '_norm']
                * cell_attrs.dcrowd**2)
        ens_species['rho_crd_' + direction + '_size'] = \
            (ens_species['rho_crd_' + direction] * cell_attrs.dcrowd**2)
        ens_species['rho_' + direction + '_uniform_sum'] = \
            (ens_species['rho_mon_' + direction + '_uniform_size']
                + ens_species['rho_crd_' + direction + '_uniform_size'])
        # sum rule: The monomer and crowder volume fractions normalized\
        # by their particle size and then added to create the sum rule:
        ens_species['phi_sumrule_' + direction] = \
            ens_species['phi_mon_' + direction + '_size']
        + ens_species['phi_crd_' + direction + '_size']
        ens_species['phi_sumrule_' + direction + '_norm'] = \
            (ens_species['phi_sumrule_' + direction]
                / (phi_mon_origin
                    + (phi_crd_infinity / cell_attrs.dcrowd)))  \
            # The first term in the denominator should be divided by 'dmon' \
        # but this operation isnot done since dmon = 1.0
        ens_species['rho_sumrule_' + direction] = \
            (ens_species['rho_mon_' + direction + '_size']
                + ens_species['rho_crd_' + direction + '_size'])
        ens_species['rho_sumrule_' + direction + '_norm'] = \
            (ens_species['rho_sumrule_' + direction]
                / (rho_mon_origin
                    + (rho_crd_infinity * cell_attrs.dcrowd**2))) \
            # The first term in the denominator should be mutiplied by\
        # 'dmon'**2 but this operation isnot done since dmon = 1.0
        # sume rule phi: Amir's way
        if direction == 'z':
            ens_species['phi_sumrule_' + direction] = \
                (ens_species['phi_mon_' + direction + '_size']
                    + ens_species['phi_crd_' + direction + '_size'])
            ens_species['phi_sumrule_' + direction + '_norm_crowd_only'] = \
                (ens_species['phi_sumrule_' + direction]
                    / (phi_crd_infinity / cell_attrs.dcrowd))
        if ens_species['phi_c_bulk'].any() == 0:
            ens_species['phi_crd_' + direction + '_norm'] = 0.0
            ens_species['phi_crd_' + direction + '_norm_size'] = 0.0
            ens_species['rho_crd_' + direction + '_norm'] = 0.0
            ens_species['rho_crd_' + direction + '_norm_size'] = 0.0
            if direction == 'z':
                ens_species['phi_sumrule_'
                            + direction + '_norm_crowd_only'] = 0.0
        # Defining concise name for ensembles anf groups
        ens_species['ens_name'] = f"N{cell_attrs.nmon}D{cell_attrs.dcyl}\
            ac{cell_attrs.dcrowd}nc{cell_attrs.ncrowd}"
        ens_species['group_name'] = f"N{cell_attrs.nmon}D{cell_attrs.dcyl}\
            ac{cell_attrs.dcrowd}"
        ens_species.drop(['filename'], axis=1, inplace=True)
        group_species.append(ens_species)
    species_database = pd.concat(group_species)
    species_database = species_database.round(round_to)
    species_database.reset_index(inplace=True, drop=True)
    if save_to is not None:
        species_database.to_csv(save_to)
    return species_database


def generator_dist_all_in_one(
                dist_db, properties_db, simulation_type, geometry, direction):
    """reads distribution files, groups them based on their types, orders \
    them based on their file names, and pass it to the \
    'distributions_all_in_one' function.

    Parameters:
    dist_db (str): Path to distribution files.
    properties_db (str): Path to ensemble-averaged properties.
    simulation_type (str): the type of particles for which the analysis is \
    one. For now, it is either 'all' or 'bug'. 'all' means all the species \
    in the system. 'bug' means the whole polymeric chain.
    geometry (str): Geomtery of the simulation box.
    direction (str): Direction along which the distributions are computed.

    Return:
    a pandas dataframe in which all the dsitributions from all the \
    simulation group are merged.
    """
    directions_shortened = {
        'radial': 'r',
        'longitudinal': 'z',
        'azimuth': 'theta'
        }
    dist_names = [
        'HistsCrd', 'HistsMon', 'PhisCrd', 'PhisMon', 'RhosCrd', 'RhosMon']
    dist_extentions = [
        directions_shortened[direction] + dist_name + '-ens_avg.csv' for
        dist_name in dist_names]
    dist_csvs = glob(dist_db)
    dist_files = organizer.sort_filenames(dist_csvs, fmts=dist_extentions)
    properties = pd.read_csv(properties_db, index_col=0)
    dist_save_to = properties_db.split(
        'all_in_one')[0] + "all_in_one-" + simulation_type + "-distributions-"
    + direction + "-ens_avg.csv"
    return all_in_one_distributions(
        dist_files, properties, geometry, directions_shortened[direction],
        save_to=dist_save_to)
