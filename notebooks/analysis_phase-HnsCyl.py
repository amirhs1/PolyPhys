from glob import glob
from polyphys.analyze import analyzer
import warnings
warnings.filterwarnings("ignore")

input_databases = glob("./N*-probe/")
project = 'HnsCylWhole'
project_details = analyzer.ANALYSIS_DETAILS_NUCLEOID[project]
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
        nonscalar_hist_t_properties=project_details[
            'nonscalar_hist_t_properties'],
        nonscalar_mat_t_properties=project_details[
            'nonscalar_mat_t_properties'],
        acf_tseries_properties=project_details[
            'acf_tseries_properties'],
        tseries_properties=project_details['tseries_properties'],
        # nlags=20000
    )


project_details = analyzer.ANALYSIS_DETAILS_ALL[project]
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
        hist_properties=project_details['hist_properties'],
        hist2d_properties=project_details['hist2d_properties'],
        hist2d_edges=project_details['hist2d_edges'],
        rho_phi_hist_properties=project_details['rho_phi_hist_properties']
    )
