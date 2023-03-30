from glob import glob
from typing import (
    Callable,
    Dict,
    Tuple,
    Optional
)
import inspect
import numpy as np
import pandas as pd
from scipy import optimize
import statsmodels.tsa.stattools as tsas
import MDAnalysis as mda

from ..manage.utilizer import invalid_keyword
from ..manage.organizer import (
    save_parent,
    sort_filenames
)
from ..manage.typer import ParserT


def acf_of_wholes(
    ensemble: pd.DataFrame,
    nlags: int,
    alpha: float
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculates the `nlags` auto correlation function (ACF) of 'wholes'
    in the `ensemble` and their lower and upper confidence limits with
    accuracy 1-`alpha`.

    Parameters
    ----------
    ensemble: DataFrame
        A ensemble's dataframe in which columns are wholes of that ensemble.
    nlags :int
        Maximum lag in the ACF.
    alpha: float
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett's formula.

    Return
    ------
    acfs: DataFrame
        A dataframe in which the columns are wholes' ACFs.
    lower_cls: DataFrame
        A dataframe in which the columns are the lower limits of the
        confidence intervals of the acfs of wholes.
    upper_cls: DataFrame
        A dataframe in which the columns are the upper limits of the
        confidence intervals of the acfs of wholes.
    """
    acfs = pd.DataFrame(columns=ensemble.columns)
    lower_cls = pd.DataFrame(columns=ensemble.columns)
    upper_cls = pd.DataFrame(columns=ensemble.columns)
    for whole in ensemble.columns:
        acf_t = tsas.acf(
            ensemble.loc[:, whole].values,
            nlags=nlags,
            alpha=alpha
        )
        acfs[whole] = acf_t[0]
        lower_cls[whole] = acf_t[1][:, 0] - acf_t[0]
        upper_cls[whole] = acf_t[1][:, 1] - acf_t[0]
    return acfs, lower_cls, upper_cls


def acf_generator(
    property_: str,
    ensembles: Dict[str, pd.DataFrame],
    nlags: int,
    alpha: float,
    group: str,
    save_to: Optional[str] = None
) -> Tuple[Dict[str, pd.DataFrame],
           Dict[str, pd.DataFrame],
           Dict[str, pd.DataFrame]]:
    """Generates the auto-correlation function (ACF) of `nlags` lags of the
    'wholes' in each of the `ensembles` of the given attribute or property
    (property_) with their associate confidence intervals (CIs) (alpha).

    Parameters
    ----------
    property_ : str
        The physical property.
    ensembles: dict of DataFrame
        a dictionary of ensembles in which each key is an ensemble name and
        its associated value is the dataframe of its wholes. The names of
        the columns in each dataframe are wholes' names.
    nlags :int
        Maximum lag in the ACF.
    alpha: float
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett's formula.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    acfs: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are
        ensembles' acfs (dataframes). In each ensemble dataframe, columns
        are wholes names.
    lower_cls: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are the
        lower limits of the confidence intervals of the acfs of ensembles.
        In each ensemble dataframe, columns are wholes' names.
    upper_cls: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are the
        upper limits of the confidence intervals of the acfs of ensembles.
        In each ensemble dataframe, columns are wholes' names.

    These three dictionaries are also written to memory as csv files if
    save_to is not None.
    """
    invalid_keyword(group, ['bug', 'nucleoid', 'all'])
    acfs = {}
    lower_cls = {}
    upper_cls = {}
    for ens_name, ens in ensembles.items():
        acf, lower_cl, upper_cl = acf_of_wholes(ens, nlags, alpha)
        acfs[ens_name] = acf
        lower_cls[ens_name] = lower_cl
        upper_cls[ens_name] = upper_cl
        if save_to is not None:
            _ = dict(
                map(
                    lambda ens: (ens[0],
                                 save_parent(
                                    ens[0],
                                    ens[1],
                                    property_ + '-acf',
                                    save_to,
                                    group=group,
                                    ext='csv'
                                    )
                                 ),
                    acfs.items()
                    )
                )
            _ = dict(
                map(
                    lambda ens: (ens[0],
                                 save_parent(
                                    ens[0],
                                    ens[1],
                                    property_ + '-acfLowerCi',
                                    save_to,
                                    group=group,
                                    ext='csv'
                                    )
                                 ),
                    lower_cls.items()
                    )
                )
            _ = dict(
                map(
                    lambda ens: (ens[0],
                                 save_parent(
                                    ens[0],
                                    ens[1],
                                    property_ + '-acfUpperCi',
                                    save_to,
                                    group=group,
                                    ext='csv'
                                    )
                                 ),
                    upper_cls.items()
                    )
                )
    return acfs, lower_cls, upper_cls


def mono_unit_exp(x: np.ndarray, omega: float, alpha: float) -> np.ndarray:
    """Calculates the mono-exponential decay with the unit amplitude, decay
    coefficient `omega` and decay exponent `alpha` for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    omega: float
        Decay coefficient
    alpha: float
        Decay exponent

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return np.exp(-1 * omega * x ** alpha)


def mono_unit_exp_tau(x: np.ndarray, tau: float, alpha: float) -> np.ndarray:
    """Calculates the mono-exponential decay with the unit amplitude, decay
    time `tau` and decay exponent `alpha` for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    tau: float
        Decay coefficient
    alpha: float
        Decay exponent

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return np.exp(-1 * (x / tau) ** alpha)


def mono_exp_tau_res(
    x: np.ndarray, tau: float, alpha: float, amp: float, residue: float
) -> np.ndarray:
    """Calculates the mono-exponential decay with the amplitude `amp`, decay
    time `tau`, decay exponent `alpha` and `residue` for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    tau: float
        Decay coefficient
    alpha: float
        Decay exponent
    amp: float
        The amplitude
    residue: float
        The residue value.

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return amp * np.exp(-1 * (x/tau) ** alpha) + residue


def mono_exp_res(
    x: np.ndarray, omega: float, alpha: float, amp: float, residue: float
) -> np.ndarray:
    """Calculates the mono-exponential decay with the amplitude `amp`,
    decay coefficient `tau`, decay exponent `alpha` and `residue`
    for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    omega: float
        The inverse of exponential decay time.
    alpha: float
        Decay exponent
    amp: float
        The amplitude
    residue: float
        The residue value.

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return amp * np.exp(-1 * omega * x ** alpha) + residue


def mono_exp(
    x: np.ndarray, omega: float, alpha: float, amp: float
) -> np.ndarray:
    """Calculates the mono-exponential decay with the amplitude `amp`,
    decay coefficient `omega`, and decay exponent `alpha` for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    omega: float
        The inverse of exponential decay time.
    alpha: float
        Decay exponent
    amp: float
        The amplitude

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return amp * np.exp(-1 * omega * x ** alpha)


def mono_exp_tau(
    x: np.ndarray, tau: float, alpha: float, amp: float
) -> np.ndarray:
    """Calculates the mono-exponential decay with the amplitude `amp`,
    decay time `tau`, and decay exponent `alpha` for array `x`.

    Parameters
    ----------
    x: array-like
        Input data
    tau: float
        Decay coefficient
    alpha: float
        Decay exponent
    amp: float
        Amplitude of decay.

    Return
    ------
    array-like:
        Measured exponential values for the input data `x`.
    """
    return amp * np.exp(-1 * (x/tau) ** alpha)


def fit_wholes(
    space: str,
    property_path: str,
    property_: str,
    func_name: str,
    property_pattern: str,
    parser: Callable,
    group: str,
    geometry: str,
    topology: str,
    x_type: str = 'index',
    scale: Optional[str] = None,
    length: int = 50000,
    save_to: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """Takes the `property_path` to the directory in which the ensemble-average
    timeseries of a given physical `property_` of a given `group` in a given
    `geometry`, and performs the following operations in the `orient` of
    interest: First, it concatenates the timeseries into one dataframe along
    the 0 or 'row' or 'index' in pandas lingo, and then adds the physical
    `attributes` of interest as the name columns to the concatenated
    timeseries.

    In each 'ensemble-averaged' dataframe, there are 3 columns with
    this name patter:
    column name = 'long_ensemble-group-property_[-measure]-stat'
    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'.

    Parameters
    ----------
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    property_path: str
        Path to the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    func_name: str
        Function fit to the data
    property_pattern: str
        The pattern by which the filenames of timeseries are started with.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    group: str in {'bug', 'all'}
        The type of the particle group.
    geometry: str in {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    topology: str in {'ring', 'linear'}
        The topology of the polymer.
    x_type: {'index', 'time'}, default 'index'
        Either use the 'index' of the data set as x variable or use the real
        'time' of simulation as x variable.
    scale: {'zscore', 'minmax}, default None,
        Whether scaled the data before fitting or not. If data is scaled, one
        of these two options used:

        'zscore':
        Scaling data by subtracting the mean of data and then dividing by the
        standard deviation; this is, using the z-score of data instead of the
        original data in fitting process.

        'minmax':
        Scaling data to [0,1] range by the min-max normalization; that is,
        reducing the min value form data and diving it by the pick-to-pick
        value (= data_max - data_min).

    length: int = 50000,
        The length of data to which `fit_func` is fit.
    save_to : str, default None
        An absolute or a relative path of a directory to which outputs are
        saved.
    **kwargs :
        Keyword arguments passed to `scipy.optimize.fit` method.

    Return
    ------
    all_in_one: pandas.DataFrame
        a dataframe in which all the timeseries are concatenated along `orient`
        of interest, and "properties and attributes" of interest are added to
        it as the new columns.

    Requirements:
    Scipy, Numpy, PolyPhys, Pandas
    """
    invalid_keyword(scale, ['zscore', 'minmax', None])
    fit_funcs = {
        'mono_unit_exp': {
            'func': mono_unit_exp,
            'params': ['omega', 'alpha']
        },
        'mono_unit_exp_tau': {
            'func': mono_unit_exp_tau,
            'params': ['tau', 'alpha']
        },
        'mono_exp_tau_res': {
            'func': mono_exp_tau_res,
            'params': ['tau', 'alpha', 'amp', 'residue']
        },
        'mono_exp_res': {
            'func': mono_exp_res,
            'params': ['omega', 'alpha', 'amp', 'residue']
        },
        'mono_exp': {
            'func': mono_exp,
            'params': ['omega', 'alpha', 'amp']
        },
        'mono_exp_tau': {
            'func': mono_exp_tau,
            'params': ['tau', 'alpha', 'amp']
        }
    }
    invalid_keyword(func_name, fit_funcs.keys())
    fit_func = fit_funcs[func_name]['func']
    fit_params = fit_funcs[func_name]['params']
    property_ext = '-' + property_ + '.csv'
    property_csvs = glob(property_path + property_pattern + property_ext)
    property_csvs = sort_filenames(property_csvs, fmts=[property_ext])
    params_std = [param + '-std' for param in fit_params]
    cols = ['whole', 'convergence'] + fit_params + params_std
    fit_data = []
    for property_csv in property_csvs:
        property_df = pd.read_csv(property_csv[0], header=0)
        # the first column of property_df is used to extract
        # the information about the property and the space it
        # belongs to.
        for col in property_df.columns:
            whole_name = col.split('-')[0]
            whole_info = parser(
                whole_name,
                'whole',
                geometry,
                group,
                topology,
                ispath=False,
            )
            whole_data = [whole_name]
            y = property_df.loc[:length, col].values
            if x_type == 'time':
                x = (np.arange(len(y)) + 1.0) * whole_info.dt
            else:
                x = (np.arange(len(y)) + 1.0)  # lags or index
            if scale == 'zscore':
                y_mean = y.mean()
                y_std = y.std()
                y = (y - y_mean) / y_std
            elif scale == 'minmax':
                y_min = y.min()
                y_max = y.max()
                y = (y - y_min) / (y_max - y_min)
            try:
                params, cov_mat = optimize.curve_fit(
                    fit_func,
                    x,
                    y,
                    **kwargs
                )
                whole_data.extend([True])
                whole_data.extend(params)
                whole_data.extend(np.diag(cov_mat))
                fit_data.append(whole_data)
            except RuntimeError:
                print("could not fit " + whole_name)
                # convergence parameter is the second element
                # in whole_data list:
                whole_data.extend([False])
                whole_data.extend(np.zeros(4))
                whole_data.extend(np.zeros(4))
                fit_data.append(whole_data)
                continue
    fit_df = pd.DataFrame(data=fit_data, columns=cols)
    if save_to is not None:
        output = '-'.join(
            ['fitReport', func_name, property_, space, 'length' + str(length)]
         )
        fit_df.to_csv(save_to + output + ".csv", index=False)
    return fit_df


def set_x_property(
    x_property: str,
    group: str,
    y_length: int,
    whole_info: ParserT
) -> np.ndarray:
    """Defines variable `x` based on `x_property` and the time difference
    between two consecutive values of a "whole" time series of a given
    simulation for which the details are given by `whole_info`.
    """
    if x_property == 'time' and group == 'bug':
        x = (np.arange(y_length) + 1.0) * whole_info.dt * whole_info.bdump
    elif x_property == 'time' and group == 'all':
        x = (np.arange(y_length) + 1.0) * whole_info.dt * whole_info.adump
    elif x_property == 'index':
        x = (np.arange(y_length) + 1.0)  # lags or index
    else:
        raise ValueError(
            f"'{x_property}' is an invalid 'x_property' or "
            f"'{group}' is an invalid 'group'."
        )
    return x


def fit_exp_wholes(
    space: str,
    space_path: str,
    property_: str,
    fit_func: Callable,
    parser: Callable,
    geometry: str,
    group: str,
    topology: str,
    alpha: float,
    omega_lag: int,
    amp_lag: int,
    res_lag: int,
    x_property: str = 'index',
    property_ext: str = 'csv',
    save_to: str = None,
    **kwargs
) -> pd.DataFrame:
    """Takes the `property_path` to the directory in which the ensemble-average
    timeseries of a given physical `property_` of a given `group` in a given
    `geometry`, and performs the following operations in the `orient` of
    interest: First, it concatenates the timeseries into one dataframe along
    the 0 or 'row' or 'index' in pandas lingo, and then adds the physical
    `attributes` of interest as the name columns to the concatenated
    timeseries.

    In each 'ensemble-averaged' dataframe, there are 3 columns with
    this name patter:
    column name = 'long_ensemble-group-property_[-measure]-stat'
    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'.

    Issues
    ------
    Currently, the `p0` argument of `scipy.optimize.fit_curve` is based on the
    assumption that the second and third arguments of the `fit_func` are the
    amplitude of

    Requirements:
    Scipy, Numpy, PolyPhys, Pandas
    """
    space_attributes = [
        'space', 'ensemble_long', 'ensemble', 'whole', 'nmon', 'dcyl',
        'dcrowd', 'phi_c_bulk']
    fit_name = ''
    if inspect.isfunction(fit_func):
        fit_name = fit_func.__name__
        # Assuming: 1. All arguments are positional. 2. The first argument is
        # the independent one. 3. Other positional arguments should be found by
        # fitting -- See scipy.optimize.fit_curve
        fit_params = inspect.getfullargspec(fit_func).args[1:]
    property_pat = '-' + property_ + '.' + property_ext
    property_dbs = glob(space_path)
    property_dbs = sort_filenames(property_dbs, fmts=[property_pat])
    params_std = ["-".join([property_, param, 'std']) for param in fit_params]
    params_name = ["-".join([property_, param]) for param in fit_params]
    stats_names = [
        "-".join([property_, measure]) for measure in ['mean', 'var', 'sem']
        ]
    cols = space_attributes + stats_names + ['convergence']
    cols = cols + params_name + params_std
    fit_data = []  # the fitting data for all wholes.
    for property_db in property_dbs:
        property_df = pd.read_csv(property_db[0], header=0)
        # the first column of property_df is used to extract
        # the information about the property and the space it
        # belongs to.
        for col in property_df.columns:
            whole_name = col.split('-')[0]
            whole_info = parser(
                whole_name,
                'whole',
                geometry,
                group,
                topology,
                ispath=False
            )
            whole_data = []  # contains fitting info of a whole
            for attr in space_attributes:
                attr_value = getattr(whole_info, attr)
                whole_data.append(attr_value)
            y = property_df.loc[:, col].to_numpy()
            y_mean = np.mean(y)
            y_var = np.var(y, ddof=1)
            y_sem = np.std(y, ddof=1) / len(y) ** 0.5
            whole_data.extend([y_mean, y_var, y_sem])
            x = set_x_property(x_property, group, len(y), whole_info)
            # Issue: p0_guess only works for
            # polyphys.analyze.correlations.mono_exp_res
            # omega sets the x interval over which fit_func decays:
            omega = 1 / x[omega_lag]
            amp = np.mean(y[:amp_lag])
            res = np.mean(y[-1 * res_lag:])
            p0_guess = [omega, alpha, amp, res]
            try:
                params_values, cov_mat = optimize.curve_fit(
                    fit_func,
                    x,
                    y,
                    p0=p0_guess,
                    **kwargs
                )
                whole_data.extend([True])
                whole_data.extend(params_values)
                whole_data.extend(np.diag(cov_mat))
                fit_data.append(whole_data)
            except RuntimeError:
                print("could not fit " + whole_name)
                # convergence parameter is the second element
                # in whole_data list:
                whole_data.extend([False])
                whole_data.extend(np.zeros(4))
                whole_data.extend(np.zeros(4))
                fit_data.append(whole_data)
                continue
    fit_df = pd.DataFrame(data=fit_data, columns=cols)
    if save_to is not None:
        output = '-'.join(
            ['fitReport', space, group, fit_name, property_]
         )
        fit_df.to_csv(save_to + output + ".csv", index=False)
    return fit_df


def bond_info(
    ag: mda.AtomGroup,
    topology: str,
    wrap: Optional[bool] = False,
    unwrap: Optional[bool] = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate all the bond lengths and the cosines of the angles between all
    the pairs of bond vectors based on the `positions` of the monomers in a
    given polymer with a given `topology`.

    `n_bonds=n_atoms-1` for 'linear' topology while `n_bonds=n_atoms` for
    'ring' topology.

    Parameters
    ----------
    ag: mda.Atomgroup
        The atomgroup for which the bond information is measured.
    topology: str
        The topology of the polymer.
    wrap : bool, optional
        If ``True``,  move all atoms to the primary unit cell before
        calculation.
    unwrap : bool, optional
        If ``True``, compounds will be unwrapped before computing their

    Return
    ------
    bond_lengths: np.ndarray
        An array of size `n_bonds` containing the bond lengths.
    cosine_corrs: np.ndarray
        An array of size `n_bonds` containing the cosines of the angles between
        all the pairs of vectors. The cosines_ij are accumulated based on their
        lags; for example, there are `nbonds` cosines with lga `i-j=1`. Hence,
        `cosine_corr[k]` is the sum of `k` cosines with lag `i-j=k`.
    """
    chain_topos: Dict[str, str] = {'linear': 'raise', 'ring': 'wrap'}
    atoms = ag.atoms
    n_bonds = atoms.n_atoms
    com = atoms.center_of_mass(wrap=wrap, unwrap=unwrap)
    pos = ag.positions
    pos = pos - com
    bonds = pos - np.take(
        pos, range(1, n_bonds+1), axis=0, mode=chain_topos[topology])
    cosine_corrs = np.zeros(n_bonds, dtype=np.float64)
    bond_lengths = np.linalg.norm(bonds, axis=1).reshape(n_bonds, 1)
    bonds = np.divide(bonds, bond_lengths)
    cosines = np.dot(bonds, bonds.T)
    for i in range(n_bonds):
        cosine_corrs[:n_bonds-i] += cosines[i, i:]
    return bond_lengths, cosine_corrs
