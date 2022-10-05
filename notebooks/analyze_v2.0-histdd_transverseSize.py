from glob import glob
from polyphys.manage import organizer
from polyphys.manage.parser import SumRule, TransFoci
from polyphys.analyze import analyzer
import warnings
warnings.filterwarnings("ignore")

input_databases = glob("./N*-probe-histdd_transverseSize")
print(input_databases)
geometry = 'biaxial'
tseries_properties_bug = [
    # property_, species, group
    ('fsdT', 'Mon', 'bug'),
    ('gyrT', 'Mon', 'bug'),
    ('rfloryT', 'Mon', 'bug'),
    ('shapeT', 'Mon', 'bug'),
    ('asphericityT', 'Mon', 'bug')
]
acf_tseries_properties_bug = [
    # property_, species, group
    #('fsdT', 'Mon', 'bug'),
    #('gyrT', 'Mon', 'bug'),
    ('transSizeT', 'Mon', 'bug'),
    #('rfloryT', 'Mon', 'bug'),
    #('shapeT', 'Mon', 'bug'),
    #('asphericityT', 'Mon', 'bug')
]

#hist_properties_bug = [
    # direction, species, group
#    ('rflory', 'Mon', 'bug')
#]
for input_database in input_databases:
    print(input_database)
    #analyzer.analyze_bug(
    analyzer.analyze_bug_trans_size(
        input_database,
        '/N*/N*',
        SumRule,
        geometry,
        True,
        #nonscalar_hist_properties=nonscalar_hist_properties,
        #tseries_properties=tseries_properties_bug,
        acf_tseries_properties=acf_tseries_properties_bug,
        #hist_properties=hist_properties_bug,
        #,
        #nlags=20000
    )

input_databases = glob("./N*-probe-histdd_transverseSize-bugSegment")
print(input_databases)
geometry = 'biaxial'
tseries_properties_bug = [
    # property_, species, group
    ('fsdT', 'Mon', 'bug'),
    ('gyrT', 'Mon', 'bug'),
    ('rfloryT', 'Mon', 'bug'),
    ('shapeT', 'Mon', 'bug'),
    ('asphericityT', 'Mon', 'bug')
]
acf_tseries_properties_bug = [
    # property_, species, group
    #('fsdT', 'Mon', 'bug'),
    #('gyrT', 'Mon', 'bug'),
    #('rfloryT', 'Mon', 'bug'),
    #('shapeT', 'Mon', 'bug'),
    ('transSizeT', 'Mon', 'bug'),
    #('asphericityT', 'Mon', 'bug')
]

#hist_properties_bug = [
    # direction, species, group
#    ('rflory', 'Mon', 'bug')
#]
for input_database in input_databases:
    print(input_database)
    analyzer.analyze_bug_trans_size(
        input_database,
        '/N*/N*',
        SumRule,
        geometry,
        True,
        #nonscalar_hist_properties=nonscalar_hist_properties,
        #tseries_properties=tseries_properties_bug,
        acf_tseries_properties=acf_tseries_properties_bug,
        #hist_properties=hist_properties_bug,
        #,
        #nlags=20000
    )

input_databases = glob("./N*-probe-*")
print(input_databases)
geometry = 'biaxial'
rho_phi_hist_properties_all = [
    # direction, species, group
    ('r', 'Crd', 'all'),
    ('r', 'Mon', 'all'),
    ('z', 'Crd', 'all'),
    ('z', 'Mon', 'all'),
]
hist_properties_all = [
    # direction, species, group
    ('theta', 'Crd', 'all'),
    ('theta', 'Mon', 'all'),
]
hist2d_properties_all = [
    # direction, species, group
    ('xy', 'Crd', 'all'),
    ('xy', 'Mon', 'all'),
    ('xz', 'Crd', 'all'),
    ('xz', 'Mon', 'all'),
    ('yz', 'Crd', 'all'),
    ('yz', 'Mon', 'all'),
]
hist2d_edges_all = [
    # direction, species, group
    #('theta', 'Crd', 'all'),
    #('theta', 'Mon', 'all'),
    ('x', 'all'),
    ('y', 'all'),
    ('z', 'all'),
]
for input_database in input_databases:
    print(input_database)
    analyzer.analyze_all(
        input_database,
        '/N*/N*',
        SumRule,
        geometry,
        False,
        #hist_properties=hist_properties_all,
        hist2d_properties=hist2d_properties_all,
        hist2d_edges=hist2d_edges_all,
        #rho_phi_hist_properties=rho_phi_hist_properties_all,
    )