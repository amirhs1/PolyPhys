"""
Does the `analysis` phase in a 'whole' probe directory.
"""
import os
from polyphys.analyze import analyzer

space_database = "./"  # the current directoyu contains the "whole" directories
absolute_path = os.path.abspath(space_database)
space_database = str(absolute_path) + "/"
project = 'HnsCubWhole'
project_details = analyzer.ANALYSIS_DETAILS_NUCLEOID[project]
analyzer.analyze_measures(
    space_database,
    project_details['hierarchy'],
    project_details['parser'],
    project_details['group'],
    project_details['geometry'],
    project_details['topology'],
    project_details['is_segment'],
    project_details['has_stamp'],
    nonscalar_hist_t_properties=project_details['nonscalar_hist_t_properties'],
    nonscalar_mat_t_properties=project_details['nonscalar_mat_t_properties'],
    acf_tseries_properties=project_details['acf_tseries_properties'],
    tseries_properties=project_details['tseries_properties']
)
project_details = analyzer.ANALYSIS_DETAILS_ALL[project]
analyzer.analyze_measures(
    space_database,
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
