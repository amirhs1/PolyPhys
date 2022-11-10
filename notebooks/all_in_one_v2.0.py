from glob import glob
import pandas as pd
import numpy as np
from polyphys.manage.parser import SumRuleCyl, TransFociCyl, TransFociCubic
from polyphys.analyze import measurer
import polyphys.api as api
import warnings
warnings.filterwarnings('ignore')

#  Choose these two before running this script:
project = 'TransFociCubic'  # 'SumRuleCyl', 'TransFociCyl'
analysis_db = './analysis/'
# List of physical properties: Set the project hierarchy
project_details = {
    'SumRuleCyl': {
        'parser': SumRuleCyl,
        'space_pat': 'N*D*ac*',
        'hierarchy': 'N*',
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
        'parser': TransFociCyl,
        'space_pat': 'ns*nl*al*D*ac*',
        'hierarchy': 'eps*',
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
    'TransFociCubic': {
        'parser': TransFociCubic,
        'space_pat': 'ns*nl*al*ac*',
        'hierarchy': 'al*',
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
    }
}
# ensemble-average stamps per project
group = 'bug'
phase = "ensAvg"
space_dbs = glob(analysis_db + project_details[project]['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        group + '-' + phase
    )
]
allInOne_stamps = []
for space_db in ens_avg_space_dbs:
    stamp_path = project_details[project]['space_hierarchy'] + 'stamps*'
    stamp_path = glob(space_db + "/" + stamp_path + '.csv')[0]
    space_stamps = pd.read_csv(stamp_path)
    allInOne_stamps.append(space_stamps)
allInOne_stamps = pd.concat(allInOne_stamps, axis=0)
allInOne_stamps.reset_index(inplace=True, drop=True)
output = analysis_db + "allInOne-" + project + "-stamps-" + phase + ".csv"
allInOne_stamps.to_csv(output, index=False)
# whole stamps per project
group = 'bug'
phase = 'ens'
space_dbs = glob(analysis_db + project_details[project]['space_pat'])
phase_space_dbs = [
    space_db for space_db in space_dbs if space_db.endswith(
        group + '-' + phase
    )
]
allInOne_stamps = []
for space_db in phase_space_dbs:
    stamp_path = project_details[project]['space_hierarchy'] + 'stamps*'
    stamp_path = glob(space_db + "/" + stamp_path + '.csv')[0]
    space_stamps = pd.read_csv(stamp_path)
    allInOne_stamps.append(space_stamps)
allInOne_stamps = pd.concat(allInOne_stamps, axis=0)
allInOne_stamps.reset_index(inplace=True, drop=True)
output = analysis_db + "allInOne-" + project + "-stamps-" + phase + ".csv"
allInOne_stamps.to_csv(output, index=False)
# equilibrium properties
group = 'bug'
spaces = glob(analysis_db + project_details[project]['space_pat'])
spaces = list(set([space.split('/')[-1].split('-')[0] for space in spaces]))
save_space = True
equili_props_wholes = api.allInOne_equil_tseries(
    project,
    analysis_db,
    group,
    spaces,
    project_details[project]['time_varying_props'],
    project_details[project]['equil_measures'],
    save_space=save_space,
    divisor=0.025,
    round_to=3,
    save_to=analysis_db,
)
ens_avg = api.allInOne_equil_tseries_ensAvg(
    project,
    equili_props_wholes,
    group,
    project_details[project]['equil_properties'],
    project_details[project]['equil_attributes'],
    save_to=analysis_db
)
