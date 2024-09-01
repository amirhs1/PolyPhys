from typing import (
    Callable,
    List,
    Tuple,
    Optional,
    Union
)
import os
import glob
import warnings
import numpy as np
import pandas as pd
from ..manage.utilizer import invalid_keyword
from ..manage.organizer import (
    sort_filenames,
    database_path,
    whole_from_file,
    whole_from_segment,
    whole_from_dist_mat_t,
    ensemble,
    ensemble_avg,
    children_stamps,
    parents_stamps
)
from ..manage.parser import (
    SumRuleCyl, TransFociCyl, TransFociCub, HnsCub, HnsCyl,
    SumRuleCubHeteroRing, SumRuleCubHeteroLinear
)
from ..manage.typer import (
    ParserT,
    TransFociT,
    TimeSeriesT,
    HistogramT,
    NonScalarHistT,
    NonScalarMatT,
    EdgeT
)
from .distributions import distributions_generator
from .correlations import acf_generator


ANALYSIS_DETAILS_NUCLEOID = {
    'SumRuleCylSegment': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'group': 'bug',
        'geometry': 'cylindrical',
        'topology': 'linear',
        'is_segment': True,
        'has_stamp': True,
        'nonscalar_mat_t_properties': [
            # property_, species, group
            ('principalT', 'Mon', 'bug'),
        ],
        'acf_tseries_properties': [
            # property_, species, group
            ('fsdT', 'Mon', 'bug'),
            ('gyrT', 'Mon', 'bug'),
            ('transSizeT', 'Mon', 'bug'),
            ('rfloryT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'SumRuleCylWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'group': 'bug',
        'geometry': 'cylindrical',
        'topology': 'linear',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_mat_t_properties': [
            # property_, species, group
            ('principalT', 'Mon', 'bug'),
        ],
        'acf_tseries_properties': [
            # property_, species, group
            ('fsdT', 'Mon', 'bug'),
            ('gyrT', 'Mon', 'bug'),
            ('transSizeT', 'Mon', 'bug'),
            ('rfloryT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'TransFociCylWhole': {
        'hierarchy': 'eps*/eps*',  # dir/file
        'parser': TransFociCyl,
        'group': 'bug',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'Foci', 'bug', 0),
            ('clustersHistT', 'Foci', 'bug', 0)
        ],
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('distMatT', 'Foci', 'bug'),
            ('principalT', 'Mon', 'bug')
        ],
        'acf_tseries_properties': [
            # property_, species, group
            ('fsdT', 'Mon', 'bug'),
            ('gyrT', 'Mon', 'bug'),
            ('transSizeT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'TransFociCubWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'bug',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'Foci', 'bug', 0),
            ('clustersHistT', 'Foci', 'bug', 0)
        ],
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('distMatT', 'Foci', 'bug'),
            ('principalT', 'Mon', 'bug')
        ],
        'acf_tseries_properties': [
            ('gyrT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'SumRuleCubHeteroLinearWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': SumRuleCubHeteroLinear,
        'group': 'bug',
        'geometry': 'cubic',
        'topology': 'linear',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'Foci', 'bug', 0),
            ('clustersHistT', 'Foci', 'bug', 0)
        ],
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('distMatT', 'Foci', 'bug'),
            ('principalT', 'Mon', 'bug')
        ],
        'acf_tseries_properties': [
            ('gyrT', 'Mon', 'bug'),
            ('rfloryT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'SumRuleCubHeteroRingWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': SumRuleCubHeteroRing,
        'group': 'bug',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'Foci', 'bug', 0),
            ('clustersHistT', 'Foci', 'bug', 0)
        ],
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('distMatT', 'Foci', 'bug'),
            ('principalT', 'Mon', 'bug')
        ],
        'acf_tseries_properties': [
            ('gyrT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'TransFociCubSegment': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'bug',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': True,
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'Foci', 'bug', 0),
            ('clustersHistT', 'Foci', 'bug', 0)
        ],
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('distMatT', 'Foci', 'bug'),
            ('principalT', 'Mon', 'bug')
        ],
        'acf_tseries_properties': [
            ('gyrT', 'Mon', 'bug'),
            ('shapeT', 'Mon', 'bug'),
            ('asphericityT', 'Mon', 'bug')
        ]
    },
    'HnsCubWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': HnsCub,
        'group': 'nucleoid',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_mat_t_properties': [
            # property_, species, group
            ('principalT', 'Mon', 'nucleoid'),
        ],
        'nonscalar_hist_t_properties': [
            # property_, species, group, avg_axis
            ('bondsHistT', 'HnsCore', 'nucleoid', 0),
            ('bondsHistDangleT', 'HnsCore', 'nucleoid', 0),
            ('bondsHistBridgeT', 'HnsCore', 'nucleoid', 0),
            ('clustersHistT', 'HnsCore', 'nucleoid', 0),
            ('clustersHistDangleT', 'HnsCore', 'nucleoid', 0),
            ('clustersHistBridgeT', 'HnsCore', 'nucleoid', 0),
            ('bondsHistDirDepT', 'HnsCore', 'nucleoid', 0),
            ('clustersHistDirDepT', 'HnsCore', 'nucleoid', 0)
        ],
        'tseries_properties': [
            # treat these two as timeseries!
            ('loopLengthHist', 'Mon', 'nucleoid'),
            ('bondCosineCorrVec', 'Mon', 'nucleoid'),
            ('bondLengthVec', 'Mon', 'nucleoid'),
            ('nBoundT', 'HnsPatch', 'nucleoid'),
            ('nFreeT', 'HnsPatch', 'nucleoid'),
            ('nEngagedT', 'HnsPatch', 'nucleoid'),
            ('nFreeT', 'HnsCore', 'nucleoid'),
            ('nBridgeT', 'HnsCore', 'nucleoid'),
            ('nDangleT', 'HnsCore', 'nucleoid'),
            ('nCisT', 'HnsCore', 'nucleoid'),
            ('nTransT', 'HnsCore', 'nucleoid')
        ],
        'acf_tseries_properties': [
            # property_, species, group
            ('gyrT', 'Mon', 'nucleoid'),
            ('shapeT', 'Mon', 'nucleoid'),
            ('asphericityT', 'Mon', 'nucleoid')
        ]
    },
    'HnsCylWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': HnsCyl,
        'group': 'nucleoid',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'is_segment': False,
        'has_stamp': True,
        'nonscalar_mat_t_properties': [
            # property_, species, group, avg_axis
            ('principalT', 'Mon', 'nucleoid')
        ],
        'tseries_properties': [
            # treat these two as timeseries!
            ('loopLengthHist', 'Mon', 'nucleoid'),
            ('bondCosineCorrVec', 'Mon', 'nucleoid'),
            ('bondLengthVec', 'Mon', 'nucleoid'),
            ('nBoundT', 'HnsPatch', 'nucleoid'),
            ('nFreeT', 'HnsPatch', 'nucleoid'),
            ('nEngagedT', 'HnsPatch', 'nucleoid'),
            ('nFreeT', 'HnsCore', 'nucleoid'),
            ('nBridgeT', 'HnsCore', 'nucleoid'),
            ('nDangleT', 'HnsCore', 'nucleoid'),
            ('nCisT', 'HnsCore', 'nucleoid'),
            ('nTransT', 'HnsCore', 'nucleoid')
        ],
        'acf_tseries_properties': [
            # property_, species, group
            ('fsdT', 'Mon', 'nucleoid'),
            ('gyrT', 'Mon', 'nucleoid'),
            ('transSizeT', 'Mon', 'nucleoid'),
            ('shapeT', 'Mon', 'nucleoid'),
            ('asphericityT', 'Mon', 'nucleoid')
        ]
    }
}
ANALYSIS_DETAILS_ALL = all_details = {
    'SumRuleCylSegment': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'group': 'all',
        'geometry': 'cylindrical',
        'topology': 'linear',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('z', 'Crd', 'all'),
            ('z', 'Mon', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('theta', 'Crd', 'all'),
            ('theta', 'Mon', 'all'),
        ],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all'),
        ]
    },
    'SumRuleCylWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'group': 'all',
        'geometry': 'cylindrical',
        'topology': 'linear',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('z', 'Crd', 'all'),
            ('z', 'Mon', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('theta', 'Crd', 'all'),
            ('theta', 'Mon', 'all'),
        ],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all'),
        ]
    },
    'TransFociCylWhole': {
        'hierarchy': 'eps*/eps*',  # dir/file
        'parser': TransFociCyl,
        'group': 'all',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
            ('z', 'Crd', 'all'),
            ('z', 'Mon', 'all'),
            ('z', 'Foci', 'all')
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
            ('z', 'Dna', 'all'),
            ('theta', 'Crd', 'all'),
            ('theta', 'Mon', 'all'),
            ('theta', 'Dna', 'all'),
            ('theta', 'Foci', 'all')
        ],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Dna', 'all'),
            ('xy', 'Foci', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Dna', 'all'),
            ('xz', 'Foci', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Dna', 'all'),
            ('yz', 'Foci', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all'),
        ]
    },
    'TransFociCubWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'all',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all')],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Dna', 'all'),
            ('xy', 'Foci', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Dna', 'all'),
            ('xz', 'Foci', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Dna', 'all'),
            ('yz', 'Foci', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all')
        ]
    },
    'SumRuleCubHeteroLinearWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'all',
        'geometry': 'cubic',
        'topology': 'linear',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all')],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Dna', 'all'),
            ('xy', 'Foci', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Dna', 'all'),
            ('xz', 'Foci', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Dna', 'all'),
            ('yz', 'Foci', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all')
        ]
    },
    'SumRuleCubHeteroRingWhole': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'all',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all')],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Dna', 'all'),
            ('xy', 'Foci', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Dna', 'all'),
            ('xz', 'Foci', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Dna', 'all'),
            ('yz', 'Foci', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all')
        ]
    },
    'TransFociCubSegment': {
        'hierarchy': 'al*/al*',  # dir/file
        'parser': TransFociCub,
        'group': 'all',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all')],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Dna', 'all'),
            ('xy', 'Foci', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Dna', 'all'),
            ('xz', 'Foci', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Dna', 'all'),
            ('yz', 'Foci', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all')
        ]
    },
    'HnsCubWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': HnsCub,
        'group': 'all',
        'geometry': 'cubic',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Hns', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Mon', 'all'),
        ],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Hns', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Hns', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Hns', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all')
        ]
    },
    'HnsCylWhole': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': HnsCyl,
        'group': 'all',
        'geometry': 'cylindrical',
        'topology': 'ring',
        'is_segment': True,
        'has_stamp': False,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Hns', 'all'),
            ('z', 'Crd', 'all'),
            ('z', 'Mon', 'all'),
            ('z', 'Hns', 'all')
        ],
        'hist_properties': [
            # direction, species, group
            ('theta', 'Crd', 'all'),
            ('theta', 'Mon', 'all'),
            ('theta', 'Hns', 'all')
        ],
        'hist2d_properties': [
            # direction, species, group
            ('xy', 'Crd', 'all'),
            ('xy', 'Mon', 'all'),
            ('xy', 'Hns', 'all'),
            ('xz', 'Crd', 'all'),
            ('xz', 'Mon', 'all'),
            ('xz', 'Hns', 'all'),
            ('yz', 'Crd', 'all'),
            ('yz', 'Mon', 'all'),
            ('yz', 'Hns', 'all'),
        ],
        'hist2d_edges': [
            # direction, group
            ('x', 'all'),
            ('y', 'all'),
            ('z', 'all'),
        ]
    }
}


def time_series(
    observations: List[str],
    parser: ParserT,
    geometry: str,
    topology: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    tseries_properties: Optional[List[TimeSeriesT]] = None,
    acf_tseries_properties: Optional[List[TimeSeriesT]] = None,
    nlags: int = 7000,
    alpha: float = 0.05
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `TimeSeriesT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `is_segment` is `True`, `observations` are "segments" and
    are "vertically" merged to create "wholes".

    Issue
    -----
    If the order of the acf_generator and the first ensemble_avg is
    change, the column names in each ensemble dataframe in 'ensembles'
    changes.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-averages are saved.
    tseries_properties: list of HistogramT
        A list of tuples where each tuple has three members: the property
        name, species, and group of a 'time-series' property.
    acf_tseries_properties: list of HistogramT
        A list of tuples where each tuple has three members: the property
        name, species, and group of a 'time-series' property. For
        `cls_tseries_properties`, the auto correlation function (AFC) is
        also computed.
    nlags: int, default 7000
        Maximum lag in the auto correlation function (AFC).
    alpha: float, default 0.05
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett's formula.
    """
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    if tseries_properties is not None:
        for property_, species, group in tseries_properties:
            tseries = sort_filenames(
                observations,
                fmts=['-' + property_ + species + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology,
                    relation='tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    topology,
                    'vector',
                    save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
            )
    if acf_tseries_properties is not None:
        for property_, species, group in acf_tseries_properties:
            tseries = sort_filenames(
                observations,
                fmts=['-' + property_ + species + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology,
                    'tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    topology,
                    'vector',
                    save_to=save_to_ens
                )
            acfs, lower_cls, upper_cls = acf_generator(
                property_ + species,
                ensembles,
                nlags,
                alpha,
                group,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acf',
                    acfs,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acfLowerCi',
                    lower_cls,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acfUpperCi',
                    upper_cls,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
                )


def histograms(
    observations: List[str],
    parser: ParserT,
    geometry: str,
    topology: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    hist_properties: Optional[List[HistogramT]] = None,
    hist_properties_no_edge: Optional[List[HistogramT]] = None,
    rho_phi_hist_properties: Optional[List[HistogramT]] = None
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `HistogramT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `segment` is `True`, `observations` are "segments" and
    are "horizontally" merged to create "wholes".

    Issue
    -----
    HThis function only work for spatial distributions created by
    SpatialDistribution class.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool, default False
        Whether `observations` are 'segment' or 'whole'
    save_to: tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-averages are saved.
    rho_phi_hist_properties: list of HistogramT, default None
        A list of tuples where each tuple has four members: the direction,
        direction long name, species, and group of a 'histogram' property.
        These histogram properties are then used to calculate the local
        number density and volume fraction.
    hist_properties: list of HistogramT default None
        A list of tuples where each tuple has three members: the direction,
        species, and group of a 'histogram' property.
    hist_properties_no_edge: list of HistogramT default None
        A list of tuples where each tuple has three members: the direction,
        species, and group of a 'histogram' property. This type of histrograms
        does not have an accopanying edge.
    """
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    # Histograms:
    # Two types of histograms with and without rhos and phis:
    # rho: local number density
    # phi: local volume fraction
    if rho_phi_hist_properties is not None:
        for direction, species, group in rho_phi_hist_properties:
            hists = sort_filenames(
                observations,
                fmts=['-' + direction + 'Hist' + species + '.npy']
            )
            edges = sort_filenames(
                observations,
                fmts=['-' + direction + 'Edge' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    topology,
                    'histogram',
                    save_to=save_to_whole
                )
                edge_wholes = whole_from_segment(
                    direction + 'Edge' + species,
                    edges,
                    parser,
                    geometry,
                    group,
                    topology,
                    'bin_edge',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group,
                    topology
                )
                edge_wholes = whole_from_file(
                    edges,
                    parser,
                    geometry,
                    group,
                    topology
                )
            # 'whole' dataframes, each with a 'whole' columns
            rho_wholes, phi_wholes = distributions_generator(
                wholes,
                edge_wholes,
                group,
                species,
                geometry,
                topology,
                direction,
                parser,
                save_to=save_to_whole
            )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                topology,
                'vector',
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                geometry,
                group,
                'dataframe',
                save_to=save_to_ens_avg
            )
            ensembles = ensemble(
                direction + 'Rho' + species,
                rho_wholes,
                parser,
                geometry,
                group,
                topology,
                'vector',
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Rho' + species,
                ensembles,
                geometry,
                group,
                'dataframe',
                save_to=save_to_ens_avg
            )
            ensembles = ensemble(
                direction + 'Phi' + species,
                phi_wholes,
                parser,
                geometry,
                group,
                topology,
                'vector',
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Phi' + species,
                ensembles,
                geometry,
                group,
                'dataframe',
                save_to=save_to_ens_avg
            )
    if hist_properties is not None:
        for direction, species, group in hist_properties:
            hists = sort_filenames(
                observations,
                fmts=['-' + direction + 'Hist' + species + '.npy']
            )
            edges = sort_filenames(
                observations,
                fmts=['-' + direction + 'Edge' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    topology,
                    'histogram',
                    save_to=save_to_whole
                )
                edge_wholes = whole_from_segment(
                    direction + 'Edge' + species,
                    edges,
                    parser,
                    geometry,
                    group,
                    topology,
                    'bin_edge',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group,
                    topology
                )
                edge_wholes = whole_from_file(
                    edges,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                topology,
                'vector',
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                geometry,
                group,
                'dataframe',
                save_to=save_to_ens_avg
            )
    if hist_properties_no_edge is not None:
        for direction, species, group in hist_properties_no_edge:
            hists = sort_filenames(
                observations,
                fmts=[direction + 'Hist' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    topology,
                    'histogram',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                topology,
                'vector',
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                geometry,
                group,
                'dataframe',
                save_to=save_to_ens_avg
            )


def histograms_2d(
    observations: List[str],
    parser: ParserT,
    geometry: str,
    topology: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    hist2d_properties: Optional[List[HistogramT]] = None,
    hist2d_edges: Optional[List[EdgeT]] = None,
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `HistogramT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `segment` is `True`, `observations` are "segments" and
    are "horizontally" merged to create "wholes".

    Issue
    -----
    HThis function only work for spatial distributions created by
    SpatialDistribution class.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool, default False
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-averages are saved.
    hist2d_properties: list of HistogramT default None
        A list of tuples where each tuple has three members: the direction,
        species, and group of a 'histogram' property.
    hist2d_edges: list of HistogramT default None
        A list of tuples where each tuple has three members: the direction,
        species, and group of a 'edge' property.
    """
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    if hist2d_properties is not None:
        for direction, species, group in hist2d_properties:
            hists = sort_filenames(
                observations,
                fmts=[direction + 'Hist' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    topology,
                    'histogram',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                topology,
                'matrix',
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                geometry,
                group,
                'ndarray',
                save_to=save_to_ens_avg
            )
    if hist2d_edges is not None:
        for direction, group in hist2d_edges:
            edges = sort_filenames(
                observations,
                fmts=[direction + 'Edge' + '.npy']
            )
            if is_segment is True:
                edge_wholes = whole_from_segment(
                    direction + 'Edge',
                    edges,
                    parser,
                    geometry,
                    group,
                    topology,
                    'bin_edge',
                    save_to=save_to_whole
                )
            else:
                edge_wholes = whole_from_file(
                    edges,
                    parser,
                    geometry,
                    group,
                    topology
                )
            ensembles = ensemble(
                direction + 'Edge',
                edge_wholes,
                parser,
                geometry,
                group,
                topology,
                'bin_edge',
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Edge',
                ensembles,
                geometry,
                group,
                'bin_edge',
                save_to=save_to_ens_avg
            )


def nonscalar_time_series(
    observations: List[str],
    parser: ParserT,
    geometry: str,
    topology: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    nonscalar_hist_t_properties: Optional[List[NonScalarHistT]] = None,
    nonscalar_mat_t_properties: Optional[List[NonScalarMatT]] = None
) -> None:
    """Runs overs all 'segment' `observations` of each of
    "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
    a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
    and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
    associated analyses to file. `nonscalar_hist_properties` are time-varying
    histograms while `nonscalar_matrix_properties` are time-varying matrices.

    A non_scalar property is either a time-varying 1D array (a vector) or a
    time-varying 2D one (a matrix). See the "Notes" and "Issues" below.

    Notes
    -----
    Currently, the vector-like properties are histograms collected over time so
    there are passed by `nonscalar_hist_properties` argument. The "avg_axis"
    ,defined below , sets the axis showing the time evolution, thus allowing
    performing time-averaging over a vector and finding the averages of its
    components.

    If the `is_segment` is `True`, `observations` are "segments" and
    are "vertically" merged to create "wholes".

    Issues
    ------
    1. Based on the Notes above, different type of nonscalar properties should
    be passed to this function and new functions should be defined for
    matrix-like nonscalar properties, so it can be possible to define the
    ensemble and ensemble-averaged versions.

    2. `nonscalar_hist_properties` is currently list of time-varying histogram
    that do not need any more process as those passed to "histograms".

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-averages are saved.
    nonscalar_hist_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has four string members. The
        first string is the name of a physical property, the second one is
        the particle type, the third one is `group` type, and the last one
        is the axis over which these physical properties are all of
        nonscalar form.
    nonscalar_mat_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    """
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    if nonscalar_hist_t_properties is not None:
        for property_, species, group, avg_axis in nonscalar_hist_t_properties:
            tseries = sort_filenames(
                    observations,
                    fmts=['-' + property_ + species + '.npy']
                )
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology,
                    'tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group,
                    topology
                )
            # Cleaning process on "clustersHistTFoci":
            # See `polyphys.analyze.clusters.count_clusters` for why we have to
            # drop the first item and what it means.
            # See `polyphys.probe.prober.trans_fuci_bug` for how
            # `count_clusters` applied.
            cluster_list = \
                ['clustersHistT', 'clustersHistDangleT', 'clustersHistBridgeT',
                 'clustersHistDirDepT']
            if property_ in cluster_list:
                wholes = {
                    whole_name: whole_array[:, 1:]
                    for whole_name, whole_array in wholes.items()
                }
            # Time-averaging process:
            wholes = {
                whole_name: np.mean(whole_array, axis=avg_axis)
                for whole_name, whole_array in wholes.items()
            }
            # changing property_ name after averaging:
            property_old = property_
            property_ = ''.join(property_.split('T'))
            warnings.warn(
                f"property name '{property_old}' changed to"
                f" '{property_}' after averaging over time.",
                UserWarning
            )
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    topology,
                    'vector',
                    save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    geometry,
                    group,
                    'dataframe',
                    save_to=save_to_ens_avg
            )
    if nonscalar_mat_t_properties is not None:
        for property_, species, group in nonscalar_mat_t_properties:
            whole_paths = sort_filenames(
                observations,
                fmts=['-' + property_ + species + '.npy']
            )
            # changing property_ name after averaging:
            # Here not only "T" is dropped but also the "property name" is
            # fully changed.
            if property_ == 'distMatT' and \
                parser.__name__ in ['TransFociCub', 'TransFociCyl']:
                wholes_hists, wholes_rdfs, wholes_tseries = \
                    whole_from_dist_mat_t(
                        whole_paths,
                        parser,
                        geometry,
                        group,
                        topology
                        )
                # Foci hists:
                # "wholes_hists" are of 'dataframe' whole_type, we directly get
                # its ensemble-averaged properties.
                # changing property_ name after averaging:
                property_new = 'pairDistHist'
                warnings.warn(
                    f"property name '{property_}' changed to"
                    f" '{property_new}' after averaging over time.",
                    UserWarning
                )
                _ = ensemble(
                        property_new + species + '-ensAvg',
                        wholes_hists,
                        parser,
                        geometry,
                        group,
                        topology,
                        'dataframe',
                        save_to=save_to_ens_avg
                )
                # Foci rdfs:
                # "wholes_rdfs" are of 'dataframe' whole_type, we directly get
                # its ensemble-averaged properties.
                # changing property_ name after averaging:
                property_new = 'pairDistRdf'
                warnings.warn(
                    f"property name '{property_}' changed to"
                    f" '{property_new}' after averaging over time.",
                    UserWarning
                )
                _ = ensemble(
                        property_new + species + '-ensAvg',
                        wholes_rdfs,
                        parser,
                        geometry,
                        group,
                        topology,
                        'dataframe',
                        save_to=save_to_ens_avg
                )
                # Foci time series:
                # "wholes_tseries" are of 'dataframe' whole_type, we directly
                # get its ensemble-averaged properties.
                # changing property_ name after averaging:
                property_new = 'pairDistT'
                warnings.warn(
                    f"property name '{property_}' changed to"
                    f" '{property_new}' after averaging over time.",
                    UserWarning
                )
                # For 'dataframe' whole_type, we directly get the
                # ensemble-averaged properties.
                _ = ensemble(
                        property_new + species + '-ensAvg',
                        wholes_tseries,
                        parser,
                        geometry,
                        group,
                        topology,
                        'dataframe',
                        save_to=save_to_ens_avg
                )
            else:  # example: the matrix principalTMon
                tseries = sort_filenames(
                    observations,
                    fmts=['-' + property_ + species + '.npy']
                )
                if is_segment is True:
                    wholes = whole_from_segment(
                        property_ + species,
                        tseries,
                        parser,
                        geometry,
                        group,
                        topology,
                        'tseries',
                        save_to=save_to_whole
                    )
                else:
                    wholes = whole_from_file(
                        tseries,
                        parser,
                        geometry,
                        group,
                        topology
                    )
                ensembles = ensemble(
                        property_ + species,
                        wholes,
                        parser,
                        geometry,
                        group,
                        topology,
                        'matrix',
                        save_to=save_to_ens
                )
                _ = ensemble_avg(
                        property_ + species,
                        ensembles,
                        geometry,
                        group,
                        'ndarray',
                        save_to=save_to_ens_avg
                )


def stamps(
    observations: List[Tuple[str]],
    group: str,
    geometry: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None],  Union[str, None],  Union[str, None]],
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `HistogramT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `segment` is `True`, `observations` are "segments" and
    are "horizontally" merged to create "wholes".

    Issue
    -----
    HThis function only work for spatial distributions created by
    SpatialDistribution class.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    group : {'bug', 'nucleoid', 'all'}, default cylindrical
        Shape of the simulation box.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    is_segment: bool
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-averages are saved.
    """
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    invalid_keyword(group, ['bug', 'nucleoid', 'all'])
    invalid_keyword(geometry,  ['cylindrical', 'slit', 'cubic'])
    if is_segment is True:
        segments_stamps = children_stamps(
            observations,
            group,
            'segment',  # lineage of the children stamps
            save_to=save_to_whole  # save all the segment stamps in one file
        )
        whole_stamps = parents_stamps(
            segments_stamps,
            geometry,
            group,
            'segment',  # lineage of the children stamps
            save_to=save_to_ens  # save all the whole stamps
        )
    else:
        whole_stamps = children_stamps(
            observations,
            group,
            'whole',  # lineage of the children stamps
            save_to=save_to_ens  # save all the whole stamps
        )
    _ = parents_stamps(
        whole_stamps,
        geometry,
        group,
        'whole',  # lineage of the children stamps
        save_to=save_to_ens_avg  # save all the ensemble-averaged stamps
    )


def analyze_measures(
    input_database: str,
    hierarchy: str,
    parser: ParserT,
    group: str,
    geometry: str,
    topology: str,
    is_segment: bool,
    has_stamp: bool,
    tseries_properties: Optional[List[TimeSeriesT]] = None,
    acf_tseries_properties: Optional[List[TimeSeriesT]] = None,
    hist_properties: Optional[List[HistogramT]] = None,
    hist_properties_no_edge: Optional[List[HistogramT]] = None,
    hist2d_properties: Optional[List[HistogramT]] = None,
    hist2d_edges: Optional[List[EdgeT]] = None,
    rho_phi_hist_properties: Optional[List[HistogramT]] = None,
    nonscalar_hist_t_properties: Optional[List[NonScalarHistT]] = None,
    nonscalar_mat_t_properties: Optional[List[NonScalarMatT]] = None,
    nlags: int = 7000,
    alpha: float = 0.05
) -> None:
    """read in the 'probe' observations of the 'group' particles based on the
    `hierarchy` of directories and files from the `input_database` path to the
    'probe' phase of a 'space' and creates the 'analysis' phase at that parent
    directory of the 'probe' of that 'space', infers 'space' names from
    `input_database` path and creates a 'space' directories at various stages
    in the 'analysis' directory for both 'bug' and 'all' groups.

    `tseries_properties`, `hists_properties`, `rho_phi_hists_properties` are
    list of tuples in which each tuple has three string members. The first
    string is the name of a physical property, the second one is the particle
    type, and the last one is `group` type.

    Parameters
    ----------
    input_database: str
        Path to the input_database; a 'space' directory at a given 'phase'.
    hierarchy: str
        Hierarchy of the directories and files within the `input_database`;
        for instance, "/N*/N*" means files that starts with "N" and are
        located in directories starting with "N".
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    group : {'bug', 'nucleoid', 'all'}, default cylindrical
        Shape of the simulation box.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool
        Whether `observations` are 'segment' (True) or 'whole' (False)
    has_stamp: bool.
        Whether `observations` have 'stamp' files (True) or 'whole' (False).
    tseries_properties: list of TimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all time-series form.
    acf_tseries_properties: list of TimeSeriesT, default None
        A list of tuples where each tuple has three members: the property
        name, species, and group of a 'time-series' property. For
        `cls_tseries_properties`, the auto correlation function (AFC) is
        also computed.
    hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all histogram form.
    hist_properties_no_edge: list of HistogramT default None
        A list of tuples where each tuple has three members: the direction,
        species, and group of a 'histogram' property. This type of histrograms
        does not have an accopanying edge.
    hist2d_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all 2D-histogram form.
    hist2d_edges: list of EdgeT, default None
        A list of tuples in which each tuple has two string members. The
        first string is the name of a physical property, and the second one is
        `group` type. These physical properties are all histogram edge
        form.
    rho_phi_hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all histogram form; however, in contrast to
        `hists_properties`, the local number density and volume fraction of
        `rho_phi_hists_properties` are also calculated.
    nonscalar_hist_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has four string members. The
        first string is the name of a physical property, the second one is
        the particle type, the third one is `group` type, and the last one
        is the axis over which these physical properties are all of
        nonscalar form.
    nonscalar_mat_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    nlags: int, default 7000
        Maximum lag in the auto correlation function (AFC).
    alpha: float, default 0.05
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett's formula.
    """
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    invalid_keyword(group, ['bug', 'nucleoid', 'all'])
    observations = glob.glob(input_database + hierarchy)
    if not observations:
        raise ValueError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    # save_to paths:
    save_to_ens = database_path(
        input_database, 'analysis', stage='ens', group=group
    )
    save_to_ens_avg = database_path(
        input_database, 'analysis', stage='ensAvg', group=group
    )
    save_to_whole = None
    if is_segment is True:
        save_to_whole = database_path(
            input_database, 'analysis', stage='wholeSim', group=group
        )
    # stamps:
    if has_stamp is True:
        stamp_files = sort_filenames(observations, fmts=['-stamps.csv'])
        stamps(
            stamp_files,
            group,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg)
        )
    # physical properties
    if tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            tseries_properties=tseries_properties,
        )
    if acf_tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            acf_tseries_properties=acf_tseries_properties,
            nlags=nlags,
            alpha=alpha
        )
    if hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist_properties=hist_properties
        )
    if hist_properties_no_edge is not None:
        histograms(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist_properties_no_edge=hist_properties_no_edge
        )
    if rho_phi_hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            rho_phi_hist_properties=rho_phi_hist_properties
        )
    if nonscalar_hist_t_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_hist_t_properties=nonscalar_hist_t_properties
        )
    if nonscalar_mat_t_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_mat_t_properties=nonscalar_mat_t_properties
        )
    if hist2d_properties is not None:
        histograms_2d(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist2d_properties=hist2d_properties
        )
    if hist2d_edges is not None:
        histograms_2d(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist2d_edges=hist2d_edges
        )


def analyze_bug_trans_size(
    input_database: str,
    hierarchy: str,
    parser: ParserT,
    geometry: str,
    topology: str,
    is_segment: bool,
    tseries_properties: Optional[List[TimeSeriesT]] = None,
    acf_tseries_properties: Optional[List[TimeSeriesT]] = None,
    hist_properties: Optional[List[HistogramT]] = None,
    rho_phi_hist_properties: Optional[List[HistogramT]] = None,
    nonscalar_hist_t_properties: Optional[List[NonScalarHistT]] = None,
    nonscalar_mat_t_properties: Optional[List[NonScalarMatT]] = None,
    nlags: int = 7000,
    alpha: float = 0.05
) -> None:
    """read in the 'probe' observations of the 'group' particles based on the
    `hierarchy` of directories and files from the `input_database` path to the
    'probe' phase of a 'space' and creates the 'analysis' phase at that parent
    directory of the 'probe' of that 'space', infers 'space' names from
    `input_database` path and creates a 'space' directories at various stages
    in the 'analysis' directory for both 'bug' and 'all' groups.

    `tseries_properties`, `hists_properties`, `rho_phi_hists_properties` are
    list of tuples in which each tuple has three string members. The first
    string is the name of a physical property, the second one is the particle
    type, and the last one is `group` type.

    Parameters
    ----------
    input_database: str
        Path to the input_database; a 'space' directory at a given 'phase'.
    hierarchy: str
        Hierarchy of the directories and files within the `input_database`;
        for instance, "/N*/N*" means files that starts with "N" and are
        located in directories starting with "N".
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    topology:
        Topology of the polymer
    is_segment: bool
        Whether `observations` are 'segment' (True) or 'whole' (False)
    nonscalar_hist_t_properties: list of NonScalarHistT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all of time-varying histogram form.
    nonscalar_mat_t_properties: list of NonScalarHistT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all of time-varying matrix form.
    tseries_properties: list of TimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all time-series form.
    acf_tseries_properties: list of TimeSeriesT, default None
        A list of tuples where each tuple has three members: the property
        name, species, and group of a 'time-series' property. For
        `cls_tseries_properties`, the auto correlation function (AFC) is
        also computed.
    hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all histogram form.
    rho_phi_hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all histogram form; however, in contrast to
        `hists_properties`, the local number density and volume fraction of
        `rho_phi_hists_properties` are also calculated.
    nonscalar_hist_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has four string members. The
        first string is the name of a physical property, the second one is
        the particle type, the third one is `group` type, and the last one
        is the axis over which these physical properties are all of
        nonscalar form.
    nonscalar_mat_t_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particle type, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    nlags: int, default 7000
        Maximum lag in the auto correlation function (AFC).
    alpha: float, default 0.05
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett's formula.
    """
    invalid_keyword(geometry, ['cylindrical', 'slit', 'cubic'])
    observations = glob.glob(input_database + hierarchy)  # observations
    if not observations:
        raise ValueError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    # 'bug' save_to paths:
    save_to_ens = database_path(
        input_database, 'analysis', stage='ens', group='bug'
    )
    save_to_ens_avg = database_path(
        input_database, 'analysis', stage='ensAvg', group='bug'
    )

    if is_segment is True:
        save_to_whole = database_path(
            input_database, 'analysis', stage='wholeSim', group='bug'
        )
    else:
        save_to_whole = None
    # physical properties
    if tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            tseries_properties=tseries_properties,
        )
    if acf_tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            acf_tseries_properties=acf_tseries_properties,
            nlags=nlags,
            alpha=alpha
        )
    if hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist_properties=hist_properties,
        )
    if rho_phi_hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            rho_phi_hist_properties=rho_phi_hist_properties,
        )
    if nonscalar_hist_t_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_hist_t_properties=nonscalar_hist_t_properties,
        )
    if nonscalar_mat_t_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            topology,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_mat_t_properties=nonscalar_mat_t_properties,
        )


def error_calc_block(
    data: np.ndarray,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """
    computes the statistical inefficiency (si) and uncertainty associated with
     a physical quantity calculated in a simulation.

    Using Flybbjerg-Peterson block average method, the plateau should be
     evident after 6-8 transformations.

    Parameters
    ----------
    data: np.array
        Input data
    save_to: bool
        Whether save results to file or not.

    Return
    ------
    block_analysis: pd.dataframe
        A pandas dataframe in which there are several columns of data about
         block-averaging error analysis.

    References
    ----------
    Original python code snippets are from "Computer simulation of liquids",
    Allen MP Tildesley DJ (2017):
    https://github.com/Allen-Tildesley/examples/blob/master/python_examples/error_calc.py

    https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py

    "Error estimates on averages of correlated data",HG Flyvbjerg and
    H Petersen. J Chem Phys, 91, 461 (1989): https://doi.org/10.1063/1.457480
    """
    nframes = len(data)  # size of the data
    var = np.var(data, ddof=1)  # Bias-corrected sample variance
    var_err = np.sqrt(2 / (nframes-1)) * var  # error in bias-corrected var
    sem = np.sqrt(var / nframes)  # correlations are neglected
    sem_err = np.sqrt(1 / (2*(nframes-1))) * sem  # error in SEM
    blocks = data.copy()
    ntransfroms = np.zeros(0, dtype=np.int8)  # number of block transformation
    ntransfroms = np.append(ntransfroms, 0)
    nblocks = np.zeros(0, dtype=np.int8)  # number of blocks
    nblocks = np.append(nblocks, nframes)
    bsize = np.zeros(0, dtype=np.int8)  # size of each block
    bsize = np.append(bsize, 1)
    bvar = np.zeros(0)  # block variances
    bvar = np.append(bvar, var)
    bvar_err = np.zeros(0)  # error in block variances
    bvar_err = np.append(bvar_err, var_err)
    bsem = np.zeros(0)  # block sems
    bsem = np.append(bsem, sem)
    bsem_err = np.zeros(0)  # error in sem
    bsem_err = np.append(bsem_err, sem_err)
    si = np.zeros(0)  # statistical inefficiency (si)
    si_initial = bsize[-1] * bvar[-1] / var  # initial si
    si = np.append(si, si_initial)
    si_err = np.zeros(0)  # error in si
    si_initial_err = np.sqrt(2 / (nframes-1)) * si_initial  # initial si error
    si_err = np.append(si_err, si_initial_err)

    while True:  # Loop over number, and hence length, of blocks
        # halving nblocks and rounding if it is odd:
        if nblocks[-1] <= 3:  # loop counter
            break
        nblocks = np.append(nblocks, nblocks[-1] // 2)
        bsize = np.append(bsize, bsize[-1] * 2)
        blocks[0:nblocks[-1]] = (
            blocks[0:2*nblocks[-1]-1:2] + blocks[1:2*nblocks[-1]:2]
            ) / 2.0  # Blocking transformation, halving the data set
        bvar = np.append(bvar, np.var(blocks[0:nblocks[-1]], ddof=1))
        bvar_err = np.append(
            bvar_err, np.sqrt(2 / (nblocks[-1]-1)) * bvar[-1]
        )
        bsem = np.append(bsem, np.sqrt(bvar[-1] / nblocks[-1]))
        bsem_err = np.append(
            bsem_err, np.sqrt((1 / (2 * (nblocks[-1]-1)))) * bsem[-1]
            )
        si = np.append(si, bsize[-1] * bvar[-1] / var)
        si_err = np.append(
            si_err, np.sqrt((1 / (2 * (nblocks[-1]-1)))) * si[-1]
            )
        ntransfroms = np.append(ntransfroms, ntransfroms[-1]+1)

    cols = [
        "ntransfroms", "bsize", "nblocks", "var", "var_err", "sem", "sem_err",
        "si", "si_err"
        ]
    block_analysis = pd.DataFrame(
        data=np.stack(
            (ntransfroms, bsize, nblocks, bvar, bvar_err, bsem, bsem_err, si,
                si_err),
            axis=1
        ),
        columns=cols
    )
    if save_to is not None:
        block_analysis.to_csv(save_to + '-block_average.csv', index=False)
    return block_analysis


def ensemble_measure(ensemble_db: str, stat_func: Callable) -> pd.DataFrame:
    """
    Apply `stat_func` column-wise (axis=0) on the "whole" columns in an
    "ensemble" dataframe or along a pre-set axis of the numpy "ensemble"
    ndarray of a given property. If the `measure_func` is applied to a numpy's
    ndarray, it must already have the information about the axis along which
    the measurements are done; for example, if the `property_db` refers to a
    ndarray with shape=(8,200,200) which say there are 8 ensembles, where each
    is 200*200 matrix, then `stat_func` should be something like this
    `stat_func=lambda x: np.mean(x, axis=(1,2))` to have 8 averages values for
    8 different ensembles.

    Parameters
    ----------
    ensemble_db: str
        Path to the "ensemble" dataframe or ndarray of a given property,

    stat_func: Callable
        A numpy or numpy-like function. If `ensemble_db` refers to ndarray,
        then `stat_func` must already have included 'axi' information.

    Return
    ------
    ensemble_measures: pd.DataFrame
        A dataframe in which the indexes are "whole" names and there is a
        column with the name of the `measure_func`. The values of this column
        are the measurements.

    Raises
    ------
    TypeError
        If `stat_func` is not a callable function.
    FileNotFoundError
        If the input `ensemble_db` path is invalid or does not exist.
    ValueError
        If the file extension is not '.csv' or '.npy'.
    Requirements
    ------------
    Pandas, Numpy, os and any other package needed for the `measure`.
    """
    file_ext = os.path.splitext(ensemble_db)[1]
    if not os.path.exists(ensemble_db):
        raise FileNotFoundError(f"The file path {ensemble_db} does not exist.")

    if not callable(stat_func):
        raise TypeError(f"The input measure {stat_func} is not callable.")

    if file_ext == '.csv':
        property_ens = pd.read_csv(ensemble_db, header=0)
        ens_measures = property_ens.apply(stat_func, axis=0)
        ens_measures = ens_measures.to_frame(name=stat_func.__name__)
    elif file_ext == '.npy':
        property_ens = np.load(ensemble_db)
        n_ens = property_ens.shape[0]
        # A property can be 0 for a given set of input parameters:
        if property_ens.size != 0:
            ens_measures = stat_func(property_ens)
        else:
            ens_measures = np.zeros(n_ens)
        space_fullname = os.path.basename(ensemble_db).split('-')[0]
        idxs = [f"{space_fullname}ens{i+1}" for i in range(n_ens)]
        ens_measures = pd.DataFrame(
            data=ens_measures,
            index=idxs,
            columns=[stat_func.__name__]
            )
    else:
        raise ValueError(
            f"'{file_ext}' is an invalid file extension. Use '.csv' for"
            "Pandas's 'dataframe' or '.npy' for Numpy's ndaarray."
            )
    return ens_measures


def space_measure(
    property_: str,
    space_db: str,
    stat_func: Callable,
    kind: str = 'dataframe'
) -> pd.DataFrame:
    """Performs `stat_func` on all the "ensembles" of a given physical
    `property_`in a "space" given. The "wholes" in an "ensemble" are time
    series. By performing the `stat_func`, we basically convert a time series
    to a single value such as a mean, standard deviation, or the like.

    It is assumed that a `property_` follow this naming convention:
        path/ensemble-shortnameTspecies.csv
    where *path* is the path to the file, *ensemble* is the name of ensemble,
    *shortname* is the short name of the property, "T" means the file is a
    time series, and "species" is the particle species to which this property
    belongs; for example
        ../ns400nl4al5D20ac1nc0-gyrTMon.csv
    is the time-varying radius of gyration of a polymer composed of "Mon"
    species in "ns400nl4al5D20ac1nc0" ensemble.

    Parameters
    ----------
    property_: str
        The name of physical property.
    space_db: str
        Path to the ensembles of a given property in a given space.
    stat_func: Callable
        A numpy or numpy-like function. If `ensemble_db` refers to ndarray,
        then `stat_func` must already have included 'axi' information.
    kind: {'dataframe', 'array'}, default 'dataframe'
        The kind of the "ensemble" file type.

    Return
    ------
    space_measure: pd.DataFrame
        A dataframe in which the indexes are all the "whole" names in a space
        and the single column are the values of applied 'measure' on that
        property.

    Errors
    ------
    ValueError:
        Raised when an invalid 'kind' parameter is provided.
    FileNotFoundError:
        Raised when no files are found in the specified path or no files
        matching the property pattern are found.
    RuntimeError:
        Raised when no valid ensembles are found for processing.

    Requirements
    ------------
    polyphys, Pandas, Numpy, or any other package needed for the `stat_func`
    """
    if kind not in ["dataframe", "array"]:
        raise ValueError(
            "Invalid 'kind' parameter. Must be 'dataframe' or 'array'."
            )

    suffix = {
        'dataframe': '.csv',
        'array': '.npy',
    }

    property_paths = glob.glob(space_db + '/*')
    if not property_paths:
        raise FileNotFoundError(
            f"No files found in the specified path: {space_db}"
            )

    property_pat = '-' + property_ + suffix[kind]  # pattern of property files.
    property_paths = sort_filenames(
        property_paths, fmts=[property_pat]
        )
    if not property_paths:
        raise FileNotFoundError(
            f"No matching files found for property pattern '{property_pat}'in"
            f"the specified path: {space_db}"
            )

    meas_name = stat_func.__name__
    last_t_index = property_.rfind('T')  # index of the last 'T'
    if last_t_index != -1:
        equil_name = "".join(
            [property_[:last_t_index], property_[last_t_index+1:]]
            )  # new name when measure applied
    else:
        warnings.warn(
            f"The index of last 'T' is '{last_t_index}', meaning "
            f"property '{property_}' is not a timeseries. Therefore, it is "
            "considered as a 'vector' property and its equilibrium name is "
            "found by splitting based on 'Vec' in its name. "
        )
        last_t_index = property_.rfind('Vec')  # index of the last 'Vec'
        equil_name = "".join(
            [property_[:last_t_index], property_[last_t_index+3:]]
            )  # new name when measure applied
    equil_meas_name = equil_name + "-" + meas_name
    spc_measure = []
    for property_db in property_paths:
        try:
            ens_measure = ensemble_measure(property_db[0], stat_func)
            spc_measure.append(ens_measure)
        except Exception as e:
            print(f"Error processing {property_db[0]}: {e}")

    if not spc_measure:
        raise RuntimeError("No valid ensembles found for processing.")

    spc_measure = pd.concat(spc_measure)
    spc_measure.rename(
        columns={meas_name: equil_meas_name}, inplace=True
    )
    return spc_measure


def equilibrium_wholes(
    space: str,
    space_db: str,
    physical_properties: List[str],
    stat_funcs: List[Callable],
    whole_stamps: pd.DataFrame,
    output_type: str = 'dataframe',
    topology: Optional[str] = None,
    output_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Performs a group of `stat_funcs` on a group of physical
    `physical_properties` in a given `space` and merges the resulting dataframe
    with the `whole_stamps` dataset.

    Each statistical function is applied to each "whole" *times series* (a
    column in an "ensemble" data frame) in each "ensemble" of a given physical
    property in a space.

    Parameters
    ----------
    space : str
        The name of the space.
    space_db : str
        The path to all the csv files in the space in globe style; for example:
        "path-to-space/*.csvs"
    physical_properties : List[str]
        The names of physical properties.
    stat_funcs : List[Callable]
        The list of statistics functions.
    whole_stamps : pd.DataFrame
        The dataframe contains the details of each "whole" simulation.
    output_type : {'dataframe', 'array'}, default 'dataframe'
        The kind of the "ensemble" file type.
    topology : str, optional
        The topology of the polymer, if any.
    output_path : str, optional
        Absolute or relative path to which the output is written.

    Returns
    -------
    whole_simulations_properties : pd.DataFrame
        A dataframe of the attributes and physical properties of all the
        "whole" simulations in a given `space`.

    Requirements
    ------------
    Pandas
    """
    equil_properties = []
    for physical_property in physical_properties:
        property_measurements = [
            space_measure(
                physical_property, space_db, stat_func, kind=output_type
            )
            for stat_func in stat_funcs
        ]
        property_measurements = pd.concat(property_measurements, axis=1)
        equil_properties.append(property_measurements)

    equil_properties = pd.concat(equil_properties, axis=1)
    equil_properties.reset_index(inplace=True)
    equil_properties.rename(columns={"index": "whole"}, inplace=True)
    if topology is not None:
        equil_properties["whole"] = equil_properties["whole"] + '.' + topology
    equil_properties = whole_stamps.merge(equil_properties, on="whole")
    if output_path is not None:
        output_filename = '-'.join([space, "whole-equilProps"])
        equil_properties.to_csv(
            output_path + output_filename + ".csv", index=False
            )
    return equil_properties
