from glob import glob
from typing import (
    Dict,
    Tuple,
)

import numpy as np
import pandas as pd
from scipy import optimize
import statsmodels.tsa.stattools as tsas
from polyphys.manage.organizer import (
    invalid_keyword,
    save_parent,
    sort_filenames
)
from polyphys.manage.parser import SumRule


def acf_of_wholes(
    ensemble: Dict[str, pd.DataFrame],
    nlags: int,
    alpha: float
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    calculates the `nlags` auto correlation function (ACF) of 'wholes' \
    in the `ensemble` and thier lower and upper confidence limits with \
    accuracy 1-`alpha`.

    Parameters
    ----------
    ensemble: DataFrame
        A ensemble's dataframe in which columns are wholes of that ensemble.
    nlags :int
        Maximum lag in the ACF.
    alpha: float
        If a number is given, the confidence intervals for the given level \
        are returned. For instance if alpha=.05, 95 % confidence intervals \
        are returned where the standard deviation is computed according to \
        Bartlett”s formula.

    Return
    ------
    acfs: DataFrame
        A dataframe in which the columns are wholes' ACFs.
    lower_cls: DataFrame
        A dataframe in which the columns are the lower limits of the \
        confidence intervals of the acfs of wholes.
    upper_cls: DataFrame
        A dataframe in which the columns are the upper limits of the \
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
    group: str = 'bug',
    save_to: str = None
) -> Tuple[Dict[str, pd.DataFrame],
           Dict[str, pd.DataFrame],
           Dict[str, pd.DataFrame]]:
    """
    generates the autocorrelation function (ACF) of `nlags` lags of the \
    'wholes' in each of the `ensembles` of the given attribute or property \
    (property_) with their associate confidence intervals (CIs) (alpha).

    Parameters
    ----------
    property_ : str
        The physical property.
    whole: dict of DataFrame
        a dictionary of ensembles in which each key is a ensemble' name and \
        its associated value is the dataframe of its wholes. The names of \
        the columns in each dataframe are wholes' names.
    nlags :int
        Maximum lag in the ACF.
    alpha: float
        If a number is given, the confidence intervals for the given level \
        are returned. For instance if alpha=.05, 95 % confidence intervals \
        are returned where the standard deviation is computed according to \
        Bartlett”s formula.
     group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    acfs: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are \
        ensembles' acfs (dataframes). In each ensemble dataframe, columns \
        are wholes names.
    lower_cls: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are the \
        lower limits of the confidence intervals of the acfs of ensembles. \
        In each ensemble dataframe, columns are wholes's names.
    upper_cls: dict of DataFrame
        Dict of ensembles where keys are ensembles' names and values are the \
        upper limits of the confidence intervals of the acfs of ensembles. \
        In each ensemble dataframe, columns are wholes's names.

    These three dictionaries are also written to memory as csv files if \
    save_to is not None.
    """
    invalid_keyword(group, ['bug', 'all'])
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


def mono_unit_exp(x, omega, alpha):
    """calculate the mono-exponential decay with the unit amplitude, decay
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
        Measured exponentail values for the input data `x`.
    """
    return np.exp(-1 * omega * x ** alpha)


def mono_unit_exp_tau(x, tau, alpha):
    """calculate the mono-exponential decay with the unit amplitude, decay
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
        Measured exponentail values for the input data `x`.
    """
    return np.exp(-1 * (x/tau) ** alpha)


def mono_exp_tau_res(x, tau, alpha, amp, residue):
    """calculate the mono-exponential decay with the amplitude `amp`, decay
    time `tau`, decay exponent `alpha` and 'residue` for array `x`.

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
        Measured exponentail values for the input data `x`.
    """
    return amp * np.exp(-1 * (x/tau) ** alpha) + residue


def mono_exp_res(x, omega, alpha, amp, residue):
    """calculate the mono-exponential decay with the amplitude `amp`,
    decay coefficient `tau`, decay exponent `alpha` and 'residue`
    for array `x`.

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
        Measured exponentail values for the input data `x`.
    """
    return amp * np.exp(-1 * omega * x ** alpha) + residue


def mono_exp(x, omega, alpha, amp):
    """calculate the mono-exponential decay with the amplitude `amp`,
    decay coefficient `omega`, and decay exponent `alpha` for array `x`.

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
        Measured exponentail values for the input data `x`.
    """
    return amp * np.exp(-1 * omega * x ** alpha)


def mono_exp_tau(x, tau, alpha, amp):
    """calculate the mono-exponential decay with the amplitude `amp`,
    decay time `tau`, and decay exponent `alpha` for array `x`.

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
        Measured exponentail values for the input data `x`.
    """
    return amp * np.exp(-1 * (x/tau) ** alpha)


def fit_wholes(
    space: str,
    property_path: str,
    property_: str,
    func_name: str,
    x_type: str = 'index',
    scale: str = None,
    length: int = 50000,
    property_pattern: str = 'N*',
    save_to: str = None,
    **kwargs
) -> pd.DataFrame:
    """take the `property_path` to the directory in which the ansemble-average
    timeseries of a given physical `property_` of a given `group` in a given
    `geometry`, and performs the following operations in the `orient` of
    interest: First, it concatenates the timeseries into one dataframe along
    the 0 or 'row' or 'index' in pandas's lingo, and thenadds the physical
    `attributes` of interest as the name columns to the concatenated
    timeseries.

    In each 'ensemble-averaged' dataframe, there are 3 columns with
    this name patter:
    column name = 'long_ensemble-group-porperty_[-measure]-stat'
    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'.

    Parameters
    ----------
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    property_path: str
        Path to the the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    func_name: str
        Function fit to the data
    x_type: {'index', 'time'}, default 'index'
        Whether use the 'index' of the data set as x variable or use the real
        'time' of simulation as x vaiable.
    scale: {'zscore', 'minmax}, default None,
        Whether scaled the data before fitting or not. If data is scaled, one
        of these two optionsis used:

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
    property_pattern: str, default 'N*'
        The pattern by which the filenames of timeseries are started with.
    save_to : str, default None
        An/a absolute/relative path of a directory to which outputs are saved.
    **kwargs :
        Keyword argumnts passed to `scipy.optimize.fit` method.

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
        # the first column of porperty_df is used to extract
        # the information about the property and the space it
        # belongs to.
        for col in property_df.columns:
            whole_name = col.split('-')[0]
            whole_info = SumRule(
                whole_name,
                geometry='biaxial',
                group='bug',
                lineage='whole',
                ispath=False
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
                # the'convergance' parameter is the second element
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
