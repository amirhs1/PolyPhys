from typing import Dict, Tuple
import pandas as pd
import statsmodels.tsa.stattools as tsas
from polyphys.manage.organizer import (
    invalid_keyword,
    save_parent
)


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
        acf_t = tsas.acf(ensemble.loc[:, whole], nlags=nlags, alpha=alpha)
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
    for ens_name, ens in ensembles.items():  \
            # bunch of simulations in a group/ensemble/ensemble-averaged group
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
                                    property_ + '-acfLowerCl',
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
                                    property_ + '-acfUpperCl',
                                    save_to,
                                    group=group,
                                    ext='csv'
                                    )
                                 ),
                    upper_cls.items()
                    )
                )
    return acfs, lower_cls, upper_cls
