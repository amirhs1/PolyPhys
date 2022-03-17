from io import TextIOWrapper
from typing import (
    List,
    Dict,
    Union
)
from polyphys.manage.organizer import invalid_keyword
from polyphys.manage.parser import SumRule


def parse_lammps_trj(
    trj_in: TextIOWrapper
) -> Dict[str, Union[float, int]]:
    """parse the header of the first snapshot of a LAMMPS trjectory dump file.

    Issue
    -----
    Is the atom_type alwasy an integer?

    Parameters
    ----------
    trj_in: _io.TextIOWrapper
        Trajectory file opened in the 'read' mode.

    Return
    ------
    snapshot: dict
        a dictionary in which there are the infromation anout the number of
        atoms and the box bounds.
    """
    _ = trj_in.readline()
    _ = trj_in.readline()
    _ = trj_in.readline()
    snapshot = dict()
    snapshot['n_atoms'] = int(trj_in.readline())  # natoms in bug trj
    _ = trj_in.readline()
    line = trj_in.readline()
    words = line.split(' ')
    snapshot['xlo'] = float(words[0])
    snapshot['xhi'] = float(words[1])
    line = trj_in.readline()
    words = line.split(' ')
    snapshot['ylo'] = float(words[0])
    snapshot['yhi'] = float(words[1])
    line = trj_in.readline()
    words = line.split(' ')
    snapshot['zlo'] = float(words[0])
    snapshot['zhi'] = float(words[1])
    _ = trj_in.readline()
    snapshot['Atoms'] = list()
    for line_number in range(snapshot['n_atoms']):
        line = trj_in.readline()
        words = line.split(" ")
        atom_id = int(words[0])
        atom_type = int(words[1])
        x = float(words[2])
        y = float(words[3])
        z = float(words[4])
        snapshot['Atoms'].append(
            [atom_id, atom_type, x, y, z]
        )
    return snapshot


def parse_lammps_data(
    data_in: TextIOWrapper,
    kind: str = 'bond',
) -> Dict[str, Union[float, int]]:
    """parse the headers of a LAMMPS trjectory data
    file of a given `kind`.

    Parameters
    ----------
    data_in: io.TextIOWrapper
        Data file opened in the 'read' mode.
    kind: {'full', 'bond', 'atom'}, default, 'bond'
        The kind of the LAMMMPS data file.

    Return
    ------
    topology: dict
        a dictionary in which there are the infromation anout the number of
        atoms and the box bounds.
    """
    data_kinds = {
        'bond': {  # order of header in bond data file:
            'atom_types': ['Masses', 'Pair_Coeffs'],
            'bond_types': ['Bond_Coeffs'],
            'n_atoms': ['Atoms', 'Velocities'],
            'n_bonds': ['Bonds']
        }
    }
    topology = dict()
    _ = data_in.readline()
    _ = data_in.readline()
    line = data_in.readline()
    words = line.split(' ')
    topology['n_atoms'] = int(words[0])
    line = data_in.readline()
    words = line.split(' ')
    topology['atom_types'] = int(words[0])
    line = data_in.readline()
    words = line.split(' ')
    topology['n_bonds'] = int(words[0])
    line = data_in.readline()
    words = line.split(' ')
    topology['bond_types'] = int(words[0])
    _ = data_in.readline()
    line = data_in.readline()
    words = line.split(' ')
    topology['xlo'] = float(words[0])
    topology['xhi'] = float(words[1])
    line = data_in.readline()
    words = line.split(' ')
    topology['ylo'] = float(words[0])
    topology['yhi'] = float(words[1])
    line = data_in.readline()
    words = line.split(' ')
    topology['zlo'] = float(words[0])
    topology['zhi'] = float(words[1])
    for header_style, header_names in data_kinds[kind].items():
        for header in header_names:
            # ignore header name and the two empty lines above and below it
            _ = data_in.readline()
            topology[header + '_fullname'] = data_in.readline()
            _ = data_in.readline()
            topology[header] = list()
            for _ in range(topology[header_style]):
                line = data_in.readline()
                topology[header].append(line)
    return topology


def write_data_for_trj(
    data_template: str,
    data_out: TextIOWrapper,
    topo_data: dict,
    trj_data: dict,
    kind: str = 'bond',
) -> None:
    """write a data file of a given `kind` to `data_out` file for the
    trajectory file for which the information is given by `trj_data` by
    combining its information with the information `topo_info` of a template
    data file with the same topology.

    Parameters
    ----------
    data_template: str
        Path to the LAMMPS data template file.
    data_out: io.TextIOWrapper
        Data file opened in the 'write' mode.
    topo_data: dict
        Toplogy information
    trj_data: dict
        Trajecory information
    kind: {'full', 'bond', 'atom'}, default, 'bond'
        The kind of the LAMMMPS data file.
    """
    data_kinds = {
        'bond': {  # order of header in bond data file:
            'atom_types': ['Masses', 'Pair_Coeffs'],
            'bond_types': ['Bond_Coeffs'],
            'n_atoms': ['Atoms', 'Velocities'],
            'n_bonds': ['Bonds']
        }
    }
    first_line = "LAMMPS data file via 'data_file_generator',"
    first_line_cont = " based on this LAMMPS template data: "
    first_line_cont_to = f"'{data_template}' \n"
    data_out.write(first_line + first_line_cont + first_line_cont_to)
    data_out.write('\n')
    data_out.write(f"{trj_data['n_atoms']} atoms\n")
    data_out.write(f"{topo_data['atom_types']} atom types\n")
    data_out.write(f"{trj_data['n_atoms'] - 1} atoms\n")
    data_out.write(f"{topo_data['bond_types']} atom types\n")
    data_out.write('\n')
    data_out.write(f"{trj_data['xlo']} {trj_data['xhi']} xlo xhi\n")
    data_out.write(f"{trj_data['ylo']} {trj_data['yhi']} ylo yhi\n")
    data_out.write(f"{trj_data['zlo']} {trj_data['zhi']} zlo zhi\n")
    for header_style, header_names in data_kinds[kind].items():
        for header in header_names:
            # ignore header name and the two empty lines above and below it
            data_out.write('\n')
            data_out.write(topo_data[header + '_fullname'])
            data_out.write('\n')
            if header == 'Atoms':
                for idx in range(topo_data[header_style]):
                    trj_line = trj_data['Atoms'][idx]
                    atom_id, atom_type, x, y, z = trj_line
                    topo_line = topo_data['Atoms'][idx]
                    words = topo_line.split(" ")
                    molecule_tag = words[1]
                    data_out.write(f"{idx + 1} {molecule_tag} {atom_type}"
                                   f" {x} {y} {z} 0 0 0\n")
            elif header == 'Velocities':
                for idx in range(topo_data[header_style]):
                    data_out.write(f"{idx + 1} 0.0 0.0 0.0\n")
            else:
                for idx in range(topo_data[header_style]):
                    data_out.write(topo_data[header][idx])


def bug_data_file_generator(
    data_template: str,
    trjs: List[str],
    geometry: str = 'biaixal',
    lineage: str = 'whole',
    save_to: str = './'
) -> None:
    """generate conjugate LAMMPS data files based on `data_template` for
    all the LAMMPS trajectory files `bug_trjs` of `lineage` type in a 'bug'
    particle group. The input and generated data files are both based of
    'bond' type; that is, they have only 'Pair Coeffs' and 'Bond Coeffs'
    and there is not 'charge' information in the system.

    To-do List
    ----------
    1. How can this method generalized to any data and trj file?
    2. Convert the parts of code about parsing trj and data_templet into
    two separate functions.

    Parameters
    ----------
    data_template: str
        Path to the LAMMPS data template file.
    trjs: list of str
        Pathes to the LAMMPS trajectory files for whih the conjugate
        data file are generated.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box
    lineage: {'segment', 'whole'}, default 'whole'
        Lineage type of children' stamps
    save_to : str, default './'
        Absolute or relative path of a directory to which outputs are saved.

    Notes
    ----
    Since the generated data file is used in visulaization and trajectory
    analysis in MDAnalysis, the toplogy information in the generated data file
    is only important. Such a topology file cannot be directly used in LAMMPS
    via the 'read_data' command since the velocities are all zero. If the
    generated is read, then there should be some command like 'velocity' in the
    LAMMPS input script that corrects the velocity information.

    The all the sections/headers above the 'Atoms' section are copied from
    the `data_template`; however, the number of atoms, the number of bonds,
    and the box bounds are inferred from the trajectory files, and replace
    their corresponding values in the data file. The 'Atoms' section is
    written based on the first TIMESTEP in the trajectory files. The
    'Velocities' section is set to zero. The 'Bonds' section is read from the
    `data_template`. There are no other sections that these three sections.

    By copying atom postions from the first TIMESTEP in the trajectory file to
    the Atoms section in the generated data file, it is ensured that the
    positions of atoms are always within the box bounds of the system inffered
    from the trajectory file itself.

    References
    ----------
    See LAMMPS website for more details about the strcuture of data and
    lammps files:
    https://docs.lammps.org/2001/data_format.html
    https://docs.lammps.org/Modify_dump.html
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(lineage, ['segment', 'whole'])
    if not data_template.endswith('data'):
        raise ValueError(
            f"The '{data_template}'"
            " file is not LAMMPS data file."
        )
    for trj in trjs:
        if not trj.endswith('lammpstrj'):
            raise ValueError(
                f"The '{trj}'"
                " is not LAMMPS trajectroy file."
            )
        trj_info = SumRule(
            trj,
            geometry=geometry,
            group='bug',
            lineage=lineage
        )
        bug_data_name = save_to + trj_info.whole + '.bug.data'
        with open(data_template, 'r') as data_in,\
                open(trj, 'r') as trj_in,\
                open(bug_data_name, 'w') as data_out:
            # ----------------------------------------------
            # infer n_atoms and the box bounds from trj file:
            # ----------------------------------------------
            trj_data = parse_lammps_trj(trj_in)
            # ------------------------------------------------------------
            # infer atom and bond types and information from data template
            # ------------------------------------------------------------
            topo_data = parse_lammps_data(data_in, kind='bond')
            # The bond data are used below.
            # ---------------------
            # write trj data header
            # ---------------------
            if topo_data['n_atoms'] != trj_data['n_atoms']:
                raise ValueError("The topology and trajectory data do not"
                                 "have the the same number of atoms.")
            if topo_data['n_bonds'] != (trj_data['n_atoms'] - 1):
                raise ValueError("The topology and trajectory data do not"
                                 "have the the same number of bonds.")
            write_data_for_trj(
                data_template,
                data_out,
                topo_data,
                trj_data,
                kind='bond'
            )
