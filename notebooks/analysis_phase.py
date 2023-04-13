from glob import glob
from polyphys.manage.parser import \
    SumRuleCyl, TransFociCyl, TransFociCub, HnsCub, HnsCyl
from polyphys.analyze import analyzer
import warnings
warnings.filterwarnings("ignore")
bug_details = {
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
        ],
        'hist_properties': [
            # direction, species, group
            ('rflory', 'Mon', 'bug')
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
        ],
        'hist_properties': [
            # direction, species, group
            ('rflory', 'Mon', 'bug')
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
            ('fsdT', 'Mon', 'bug'),
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
        'acf_tseries_properties': [
            ('gyrT', 'Mon', 'nucleoid'),
            ('shapeT', 'Mon', 'nucleoid'),
            ('asphericityT', 'Mon', 'nucleoid')
        ],
        'tseries_properties': [
            ('bondCosineCorrVec', 'Mon', 'nucleoid'),
            ('bondLengthVec', 'Mon', 'nucleoid')
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
            ('nBound', 'HnsPatchT', 'nucleoid'),
            ('nFree', 'HnsPatchT', 'nucleoid'),
            ('nEngaged', 'HnsPatchT', 'nucleoid'),
            ('nFree', 'HnsCoreT', 'nucleoid'),
            ('nBridge', 'HnsCoreT', 'nucleoid'),
            ('nDangle', 'HnsCoreT', 'nucleoid'),
            ('nBound', 'MonT', 'nucleoid'),
            ('nCis', 'MonT', 'nucleoid'),
            ('nTrans', 'MonT', 'nucleoid')
        ], 
        'acf_tseries_properties': [
            # property_, species, group
            ('fsdT', 'Mon', 'nucleoid'),
            ('gyrT', 'Mon', 'nucleoid'),
            ('transSizeT', 'Mon', 'nucleoid'),
            ('shapeT', 'Mon', 'nucleoid'),
            ('asphericityT', 'Mon', 'nucleoid')
        ],
        'hist_properties_no_edge' : [
            # direction, species, group
            ('loopLength', 'Mon', 'nucleoid')
        ]
    }
}

#input_databases = glob("/Users/amirhsi_mini/research_data/probe/N*-probe
# -bugSegment")
#input_databases = glob("/Users/amirhsi_mini/research_data/probe/ns*-probe/")
input_databases = glob("./N*-probe/")
# /N*-nucleoidWhole/")
#print(input_databases)
project = 'HnsCylWhole'
project_details = bug_details[project]
for input_database in input_databases:
    print(input_database)
    analyzer.analyze_measures(
        input_database,
        project_details['hierarchy'],
        project_details['parser'],
        project_details['group'],
        project_details['geometry'],
        project_details['topology'],
        project_details['is_segment'],
        project_details['has_stamp'],
        nonscalar_hist_t_properties=project_details['nonscalar_hist_t_properties'],
        nonscalar_mat_t_properties=project_details['nonscalar_mat_t_properties'],
        acf_tseries_properties=project_details['acf_tseries_properties']
        #tseries_properties=project_details['tseries_properties']
        #nlags=20000
    )
