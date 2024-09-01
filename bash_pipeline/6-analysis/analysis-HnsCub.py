"""
Does the `analysis` phase in a 'whole' probe directory.
"""
import os
from polyphys.analyze import analyzer

space_database = "./"  # the current directoyu contains the "whole" directories
absolute_path = os.path.abspath(space_database)
space_database = str(absolute_path) + "/"
project = 'HnsCubWhole'
bug_details = analyzer.ANALYSIS_DETAILS_NUCLEOID[project]
analyzer.analyze_measures(
    space_database,
    details_bug['hierarchy'],
    details_bug['parser'],
    details_bug['group'],
    details_bug['geometry'],
    details_bug['topology'],
    details_bug['is_segment'],
    details_bug['has_stamp'],
    nonscalar_mat_t_properties=details_bug[
        'nonscalar_mat_t_properties'],
    acf_tseries_properties=details_bug[
        'acf_tseries_properties'],
    tseries_properties=details_bug['tseries_properties']
)
all_details = analyzer.ANALYSIS_DETAILS_ALL[project]
analyzer.analyze_measures(
    space_database,
    details_all['hierarchy'],
    details_all['parser'],
    details_all['group'],
    details_all['geometry'],
    details_all['topology'],
    details_all['is_segment'],
    details_all['has_stamp'],
    hist_properties=details_all['hist_properties'],
    hist2d_properties=details_all['hist2d_properties'],
    hist2d_edges=details_all['hist2d_edges'],
    rho_phi_hist_properties=details_all['rho_phi_hist_properties']
)
