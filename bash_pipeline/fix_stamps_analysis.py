import glob
from polyphys.manage.organizer import sort_filenames, database_path

from polyphys.analyze import analyzer

input_databases = glob.glob("./probes/N*-probe/")
project = 'HnsCubWhole'
project_details = analyzer.ANALYSIS_DETAILS_NUCLEOID[project]
hierarchy = project_details['hierarchy']
group = project_details['group']
geometry = project_details['geometry']
is_segment = project_details['is_segment']

for input_database in input_databases:
    observations = glob.glob(input_database + hierarchy)
    save_to_ens = database_path(
        input_database, 'analysis', stage='ens', group=group
    )
    save_to_ens_avg = database_path(
        input_database, 'analysis', stage='ensAvg', group=group
    )
    save_to_whole = None
    stamp_files = sort_filenames(observations, fmts=['-stamps.csv'])
    analyzer.stamps(
        stamp_files,
        group,
        geometry,
        is_segment,
        (save_to_whole, save_to_ens, save_to_ens_avg)
    )