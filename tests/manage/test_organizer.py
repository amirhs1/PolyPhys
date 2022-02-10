from glob import glob
from polyphys.analyze import analyzer
from polyphys.manage.organizer import (
    camel_case_split,
    isfloat,
    sort_by_alphanumeric,
    sort_filenames,
    invalid_keyword,
    save_parent,
    whole,
    ensemble,
    ensemble_avg,
    children_stamps,
    parents_stamps,
    time_series,
    histograms,
    database_path,
    analyze_segments
)


def test_camel_case_split():
    pass


def test_isfloat():
    pass


def test_sort_by_alphanumeric():
    pass


def test_sort_filenames():
    pass


def test_invalid_keyword():
    pass


def test_parent_histogram():
    pass


def test_save_parent():
    pass

def test_database_path():
    input_database = '../test_data/probe/N500D10.0ac0.8'
    phase = 'analysis'
    stage = 'ens'
    group = None
    database_path(input_database, phase, stage=stage, group=group)

def test_whole():
    input_database = '../test_data/probe/N500D10.0ac0.8-segment'
    geometry = 'biaxial'
    phase = 'analysis'
    stage = 'wholeSim'
    group = 'bug'
    hierarchy = '/N*/N*'
    lineages = ('segment', 'whole')
    relation = 'time_series'
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise OSError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    observations = sort_filenames(observations, fmts=['bug-gyrTMon.npy'])
    save_to = database_path(input_database, phase=phase, stage=stage, group=group)
    gyr_wholes = whole('gyrTMon', observations, geometry=geometry, group=group, relation=relation, save_to=save_to)

    input_database = '../test_data/probe/N500D10.0ac0.8-segment'
    geometry = 'biaxial'
    phase = 'analysis'
    stage = 'wholeSim'
    group = 'all'
    hierarchy = '/N*/N*'
    lineages = ('segment', 'whole')
    relation = 'histogram'
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise OSError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    observations = sort_filenames(observations, fmts=['-all-rHistMon.npy'])
    save_to = database_path(
        input_database, phase=phase, stage=stage, group=group)
    rhist_wholes = whole(
        'rHistMon', observations, geometry=geometry, group=group,
        relation=relation, save_to=save_to)

    input_database = '../test_data/probe/N500D10.0ac0.8-segment'
    geometry = 'biaxial'
    phase = 'analysis'
    stage = 'wholeSim'
    group = 'all'
    hierarchy = '/N*/N*'
    relation = 'bin_edges'
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise OSError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    observations = sort_filenames(observations, fmts=['-all-rEdgeMon.npy'])
    save_to = database_path(input_database, phase=phase, stage=stage, group=group)
    redge_wholes = whole(
        'rEdgeMon', observations, geometry=geometry,
        group=group, relation=relation, save_to=save_to)


def test_ensmeble():
    group = 'bug'
    input_database = '../test_data/probe/N500D10.0ac0.8-segment'
    phase = 'analysis'
    stage = 'ens'
    save_to = database_path(
        input_database, phase=phase, stage=stage, group=group)
    ensembles = ensemble(
        'gyrTMon',
        gyr_wholes,
        group=group,
        edge_wholes=None,
        save_to=save_to)
    group = 'all'
    input_database = '../test_data/probe/N500D10.0ac0.8-segment'
    phase = 'analysis'
    stage = 'ens'
    save_to = database_path(input_database, phase=phase, stage=stage, group=group)
    ensembles = ensemble(
        'rHistMon',
        rhist_wholes,
        group=group,
        edge_wholes=redge_wholes,
        save_to=save_to)


def test_siblings_statistics():
    pass


def test_children_stamps():
    pass


def test_children_stamps():
    input_database = '../test_data/probe/N500D10.0ac0.8*'
    hierarchy = '/N*/N*'
    lineage = 'segment'
    group = 'bug'
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise OSError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    observations = sort_filenames(observations, fmts=['-stamps.csv'])
    stamps = children_stamps(observations, group=group, lineage=lineage)


def test_parents_stamps():
    input_database = '../test_data/probe/N500D10.0ac0.8*'
    geometry = 'biaxial'
    group = 'bug'
    lineage = 'segment'
    hierarchy = '/N*/N*'
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise OSError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    observations = sort_filenames(observations, fmts=['-stamps.csv'])
    segments_stamps = children_stamps(
        observations, group=group, lineage=lineage
    )
    wholes_stamps = parents_stamps(
        segments_stamps, geometry=geometry, group=group, lineage=lineage
    )


