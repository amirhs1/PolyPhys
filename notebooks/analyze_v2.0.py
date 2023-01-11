from glob import glob
from polyphys.manage.parser import SumRuleCyl, TransFociCyl, TransFociCubic
from polyphys.analyze import analyzer
import warnings
warnings.filterwarnings("ignore")
#  inputs: check the two following before running this script
input_databases = glob("./ns*-bugWhole/")
project = 'TransFociCubicWhole'
#  "bug" files details in various projects
bug_details = {
    'SumRuleCylSegment': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'geometry': 'cylindrical',
        'is_segment': True,
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
        'geometry': 'cylindrical',
        'is_segment': False,
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
        'geometry': 'cylindrical',
        'is_segment': False,
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
        'geometry': 'cubic',
        'is_segment': False,
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
    }
}
#  analyzing "bug" files
project_details = bug_details[project]
for input_database in input_databases:
    print(input_database)
    analyzer.analyze_bug(
        input_database,
        project_details['hierarchy'],
        project_details['parser'],
        project_details['geometry'],
        project_details['is_segment'],
        nonscalar_hist_t_properties=project_details[
            'nonscalar_hist_t_properties'],
        nonscalar_mat_t_properties=project_details[
            'nonscalar_mat_t_properties'],
        acf_tseries_properties=project_details['acf_tseries_properties']
    )
#  "all" files details in different projects
all_details = {
    'SumRuleCylSegment': {
        'hierarchy': 'N*/N*',  # dir/file
        'parser': SumRuleCyl,
        'geometry': 'cylindrical',
        'is_segment': True,
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
        'geometry': 'cylindrical',
        'is_segment': True,
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
        'geometry': 'cubic',
        'is_segment': True,
        'rho_phi_hist_properties': [
            # direction, species, group
            ('r', 'Crd', 'all'),
            ('r', 'Mon', 'all'),
            ('r', 'Foci', 'all'),
        ],
        'hist_properties': [
            # direction, species, group
            ('r', 'Dna', 'all'),
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
            ('z', 'all')
        ]
    }
}
# analyzing "all" files
project_details = all_details[project]
for input_database in input_databases:
    print(input_database)
    analyzer.analyze_all(
        input_database,
        project_details['hierarchy'],
        project_details['parser'],
        project_details['geometry'],
        project_details['is_segment'],
        hist_properties=project_details['hist_properties'],
        hist2d_properties=project_details['hist2d_properties'],
        hist2d_edges=project_details['hist2d_edges'],
        #rho_phi_hist_properties=project_details['rho_phi_hist_properties']
    )
