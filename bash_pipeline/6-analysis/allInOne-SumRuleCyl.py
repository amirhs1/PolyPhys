"""
Does the 'allInOne' step in the 'analysis' directory.
"""
import os
from glob import glob
from itertools import product
import numpy as np
import pandas as pd
from polyphys.manage import organizer
import polyphys.api as api
from polyphys.api import PROJECTS_DETAILS as PSD

project = 'SumRuleCyl'
analysis_db = './'
absolute_path = os.path.abspath(analysis_db)
analysis_db = str(absolute_path) + "/"
project_details = PSD[project]

# Ensemble stamps
print("Analyzing ensemble stamps ...")

space_dbs = glob(analysis_db + project_details['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        project_details['group'] + '-ensAvg'
    )
]
allInOne_stamps = []
for space_db in ens_avg_space_dbs:
    print(space_db)
    stamp_path = project_details['space_hierarchy'] + 'stamps*'
    stamp_path = glob(space_db + "/" + stamp_path + '.csv')[0]
    space_stamps = pd.read_csv(stamp_path)
    allInOne_stamps.append(space_stamps)
allInOne_stamps = pd.concat(allInOne_stamps, axis=0)
allInOne_stamps.reset_index(inplace=True, drop=True)
output = analysis_db + "allInOne-" + project + "-stamps-ensAvg.csv"
allInOne_stamps.to_csv(output, index=False)
# Whole stamps
print("Analyzing whole stamps ...")
space_dbs = glob(analysis_db + project_details['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        project_details['group'] + '-ens'
    )
]
allInOne_stamps = []
for space_db in ens_avg_space_dbs:
    stamp_path = project_details['space_hierarchy'] + 'stamps*'
    stamp_path = glob(space_db + "/" + stamp_path + '.csv')[0]
    space_stamps = pd.read_csv(stamp_path)
    allInOne_stamps.append(space_stamps)
allInOne_stamps = pd.concat(allInOne_stamps, axis=0)
allInOne_stamps.reset_index(inplace=True, drop=True)
output = analysis_db + "allInOne-" + project + "-stamps-ens.csv"
allInOne_stamps.to_csv(output, index=False)

# Auto-correlation functions
print("Analyzing ACFs ...")
phase = 'ensAvg'
space_dbs = glob(analysis_db + project_details['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        project_details['group'] + '-' + phase
    )
]
# list of unique property_measures:
filepath = ens_avg_space_dbs[0] + '*' + project_details['hierarchy'] + '.csv'
_, uniq_props_measures = organizer.unique_property(
    filepath, 2, ["-" + phase], drop_properties=['stamps'])
if project == 'TransFociCyl':
    uniq_props_measures.remove('transSizeTMon-acf')
    uniq_props_measures.remove('transSizeTMon-acfLowerCi')
    uniq_props_measures.remove('transSizeTMon-acfUpperCi')
print(uniq_props_measures)
for ens_avg_space_db in ens_avg_space_dbs:
    ens_avgs = list()
    space = ens_avg_space_db.split('/')[-2].split('-')[0]
    for property_ in uniq_props_measures:
        ens_avg = organizer.space_tseries(
            ens_avg_space_db,
            property_,
            project_details['parser'],
            project_details['hierarchy'],
            project_details['attributes'],
            project_details['group'],
            project_details['geometry'],
            project_details['topology'],
            divisor=project_details['divisor'],
            is_save=False  # if True, save per property per space
        )
        if project in ['HnsCyl', 'HnsCub']:
            ens_avg['phi_c_bulk_round'].replace(0.09, 0.08, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.15, 0.16, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.21, 0.2, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.31, 0.32, inplace=True)
            ens_avg = \
                ens_avg.loc[~ens_avg['phi_c_bulk_round'].isin([0.06, 0.18]), :]
        elif project in ['TransFociCyl', 'TransFociCub']:
            ens_avg = \
                ens_avg.loc[~ens_avg['phi_c_bulk_round'].isin(
                    [0.025, 0.05, 0.075, 0.125, 0.175]), :]
        elif project in ['SumRuleCyl']:
            pass
        else:
            raise ValueError("The 'phi_c drop' condition is not defined for "
                             f"'{project}' project.")
        ens_avgs.append(ens_avg)
        del ens_avg
    ens_avgs = pd.concat(ens_avgs, axis=1)
    # drop duplicated columns:
    ens_avgs = ens_avgs.loc[:, ~ens_avgs.columns.duplicated()]
    output_name = analysis_db + "-".join(
        [space,
         project_details['group'],
         "chainSize-acf.parquet.brotli"
         ]
    )
    ens_avgs.to_parquet(output_name, index=False, compression='brotli')
    del ens_avgs


# Time series
print("Analyzing time series ...")
phase = 'ensAvg'
space_dbs = glob(analysis_db + project_details['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        project_details['group'] + '-' + phase
    )
]
print(ens_avg_space_dbs)
# list of unique property_measures:
filepath = ens_avg_space_dbs[0] + '*' + project_details['hierarchy'] + '.csv'
_, uniq_props_measures = organizer.unique_property(
    filepath, 2, ["-" + phase], drop_properties=['stamps'])
props_tseries = list(
    set(
        [prop.split("-acf")[0] for prop in uniq_props_measures]
    )
)

if project == 'TransFociCyl':
    props_tseries.remove('transSizeTMon')
print(props_tseries)

for ens_avg_space_db in ens_avg_space_dbs:
    ens_avgs = list()
    space = ens_avg_space_db.split('/')[-2].split('-')[0]
    for property_ in props_tseries:
        ens_avg = organizer.space_tseries(
            ens_avg_space_db,
            property_,
            project_details['parser'],
            project_details['hierarchy'],
            project_details['attributes'],
            project_details['group'],
            project_details['geometry'],
            project_details['topology'],
            divisor=project_details['divisor'],
            is_save=False  # if True, save per property per space
        )
        if project in ['HnsCyl', 'HnsCub']:
            ens_avg['phi_c_bulk_round'].replace(0.09, 0.08, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.15, 0.16, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.21, 0.2, inplace=True)
            ens_avg['phi_c_bulk_round'].replace(0.31, 0.32, inplace=True)
            ens_avg = \
                ens_avg.loc[~ens_avg['phi_c_bulk_round'].isin([0.06, 0.18]), :]
        elif project in ['TransFociCyl', 'TransFociCub']:
            ens_avg = \
                ens_avg.loc[~ens_avg['phi_c_bulk_round'].isin(
                    [0.025, 0.05, 0.075, 0.125, 0.175]), :]
        elif project in ['SumRuleCyl']:
            pass
        else:
            raise ValueError("The 'phi_c drop' condition is not defined for "
                             f"'{project}' project.")
        ens_avgs.append(ens_avg)
        del ens_avg
    ens_avgs = pd.concat(ens_avgs, axis=1)
    # drop duplicated columns:
    ens_avgs = ens_avgs.loc[:, ~ens_avgs.columns.duplicated()]
    output_name = analysis_db +  "-".join(
        [space,  project_details['group'], "chainSize.parquet.brotli"]
    )
    ens_avgs.to_parquet(output_name, index=False, compression='brotli')
    del ens_avgs

# Ensemble-averaged time-averaged properties
print("Analyzing ensemble-averaged time-averaged properties ...")
spaces = glob(analysis_db + project_details['space_pat'])
spaces = sorted(
    list(set([space.split('/')[-1].split('-')[0] for space in spaces])))
save_space = True
equili_props_wholes = api.all_in_one_equil_tseries(
    project,
    analysis_db,
    project_details['group'],
    spaces,
    project_details['time_varying_props'],
    project_details['equil_measures'],
    save_space=save_space,
    divisor=project_details['divisor'],
    round_to=3,
    kind='dataframe',
    save_to=analysis_db,
)
ens_avg = api.all_in_one_equil_tseries_ens_avg(
    project,
    equili_props_wholes,
    project_details['group'],
    project_details['equil_properties'],
    project_details['equil_attributes'],
    save_to=analysis_db
)

# Sum-Rule
print("Analyzing spatial distributions ...")
phase = 'ensAvg'
group = 'all'
space_dbs = glob(analysis_db + project_details['space_pat'])
ens_avg_space_dbs = [
    space_db + "/" for space_db in space_dbs if space_db.endswith(
        group + '-' + phase
    )
]
print(ens_avg_space_dbs)
species_dict = project_details['rhosPhisNormalizedScaled']
print('species_dict: ', project_details['rhosPhisNormalizedScaled'])
dir_prop_pairs = list(
    product(project_details['props'],
            project_details['directions'])
)
print('dir_prop_pairs: ', dir_prop_pairs)

for (prop, direction) in dir_prop_pairs:
    all_in_one_list = list()
    for (species, size_attr) in species_dict:
        per_species_list = list()
        for ens_avg_space_db in ens_avg_space_dbs:
            space = ens_avg_space_db.split('/')[-2].split('-')[0]
            per_space = organizer.space_sum_rule(
                ens_avg_space_db,
                prop,
                project_details['parser'],
                project_details['hierarchy'],
                project_details['attributes'],
                species,
                size_attr,
                group,
                project_details['geometry'],
                project_details['topology'],
                direction,
                divisor=project_details['divisor'],
                is_save=False
            )
            if project in ['HnsCyl', 'HnsCub']:
                per_space['phi_c_bulk_round'].replace(0.09, 0.08, inplace=True)
                per_space['phi_c_bulk_round'].replace(0.15, 0.16, inplace=True)
                per_space['phi_c_bulk_round'].replace(0.21, 0.2, inplace=True)
                per_space['phi_c_bulk_round'].replace(0.31, 0.32, inplace=True)
                per_space = \
                    per_space.loc[~per_space['phi_c_bulk_round'].isin(
                        [0.06, 0.18]), :]
            elif project in ['TransFociCyl', 'TransFociCub']:
                per_space = \
                    per_space.loc[~per_space['phi_c_bulk_round'].isin(
                        [0.025, 0.05, 0.075, 0.125, 0.175]), :]
            elif project in ['SumRuleCyl']:
                pass
            else:
                raise ValueError(
                    "The 'phi_c drop' condition is not defined for"
                    f" '{project}' project."
                    )
            per_species_list.append(per_space)
        per_species = pd.concat(per_species_list, axis=0)
        per_species = per_species.loc[:, ~per_species.columns.duplicated()]
        per_species.reset_index(inplace=True, drop=True)
        all_in_one_list.append(per_species)
    all_in_one = pd.concat(all_in_one_list,axis=1)
    all_in_one = all_in_one.loc[:, ~all_in_one.columns.duplicated()]
    all_in_one.reset_index(inplace=True, drop=True)
    output = '-'.join(['allInOne', project, group, direction + prop])
    output += '-NormalizedScaled.parquet.brotli'
    output = analysis_db + output
    all_in_one.to_parquet(output, index=False, compression='brotli')
    del all_in_one, all_in_one_list