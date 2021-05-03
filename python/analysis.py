"""
This file sets-up and arranges the data to be used
"""
import time
import numpy as np
from classes import Atom, Molecule


def count_atoms(file: str = "conf.gro") -> int:
    """
    Reads the second line of the file which contains the number of atoms
    :param file: str, name of the file to read (conf.gro by default)
    :return: nrow: int, number of atoms in the file
    """
    nrow = np.loadtxt(file, skiprows=1, max_rows=1, usecols=0, dtype=int)

    return int(nrow)


def load_atoms(nrow: int, file: str = "conf.gro") -> np.ndarray:
    """
    Reads all the atoms and their coordinates in a given .gro file, then sorts
    them in the Atom dataclass
    :param nrow: int, number of atoms contained in the file
    :param file: str, name of the file to read (conf.gro by default)
    :return: a ndarray containing all the Atom dataclasses
    """
    names = np.loadtxt(file, skiprows=2, usecols=1, max_rows=nrow,
                       dtype=str)
    ids = np.loadtxt(file, skiprows=2, usecols=2, max_rows=nrow,
                     dtype=str)
    x = np.loadtxt(file, skiprows=2, usecols=3, max_rows=nrow)
    y = np.loadtxt(file, skiprows=2, usecols=4, max_rows=nrow)
    z = np.loadtxt(file, skiprows=2, usecols=5, max_rows=nrow)

    atoms_list = []

    for i in range(nrow):
        atoms_list.append(Atom(names[i], ids[i], x[i], y[i], z[i]))

    return np.array(atoms_list, dtype=Atom)


def compute_molecules(nrow: int,
                      atoms_list: np.ndarray,
                      file: str = "conf.gro") -> np.ndarray:
    """
    Reads the first column of a given .gro file containing the names of the
    molecules and sorts them in the Molecule dataclass along with the
    corresponding Atoms
    :param nrow: int, number of atoms contained in the file
    :param atoms_list: np.ndarray, contains the Atoms to be sorted in Molecules
    :param file: str, name of the file to read (conf.gro by default)
    :return: a ndarray containing the Molecule dataclasses
    """
    mol = np.loadtxt(file, skiprows=2, usecols=0, max_rows=nrow,
                     dtype=str)

    mols_list = []

    counter = 0
    wat_count = 1
    met_count = 1

    while True:
        if "SOL" in mol[counter]:
            mols_list.append(Molecule("WAT-" + "{:0>4}".format(wat_count),
                                      counter // 4,
                                      (atoms_list[counter],
                                       atoms_list[counter + 1],
                                       atoms_list[counter + 2],
                                       atoms_list[counter + 3])))
            wat_count += 1
            counter += 4

        elif "CH4" in mol[counter]:
            mols_list.append(Molecule("MET-" + "{:0>4}".format(met_count),
                                      counter // 4,
                                      (atoms_list[counter],
                                       atoms_list[counter + 1],
                                       atoms_list[counter + 2],
                                       atoms_list[counter + 3],
                                       atoms_list[counter + 4],)))
            met_count += 1
            counter += 5

        if counter + 4 >= nrow:
            break

    return np.array(mols_list, dtype=Molecule)


def print_atom(atom: Atom):
    """
    Pretty printer for the Atom dataclass
    :param atom: the Atom dataclass to print
    :return: None
    """
    print("name: {}\n"
          "id: {:0>4}\n"
          "coordinates:\n"
          "    x: {}\n"
          "    y: {}\n"
          "    z: {}\n".format(atom.name, atom.id, atom.x, atom.y, atom.z))


def print_mol(mol: Molecule):
    """
    Pretty printer for the Molecule dataclass
    :param mol: the Molecule to print
    :return: None
    """
    atoms_name = ''
    atoms_id = ''
    atoms_x = ''
    atoms_y = ''
    atoms_z = ''

    for i in range(len(mol.contains)):
        atoms_name += "{:>2}atom name: {:<5}".format('', mol.contains[i].name)
        atoms_id += "{0:>2}atom id: {1:0>4}{0:>3}".format('',
                                                          mol.contains[i].id)
        atoms_x += "{0:<4}x: {1:>6.3f}{0:>6}".format('', mol.contains[i].x)
        atoms_y += "{0:<4}y: {1:>6.3f}{0:>6}".format('', mol.contains[i].y)
        atoms_z += "{0:<4}z: {1:>6.3f}{0:>6}".format('', mol.contains[i].z)

    coordinates = "{0:>2}{1:<18}{1:<18}{1:<18}{1:<18}".format('',
                                                              'coordinates:')

    print("mol name: {}\n"
          "id: {:0>4}\n"
          "contains:".format(mol.name, mol.id))
    print(atoms_name.rstrip() + '\n' +
          atoms_id.rstrip() + '\n' +
          coordinates.rstrip() + '\n' +
          atoms_x.rstrip() + '\n' +
          atoms_y.rstrip() + '\n' +
          atoms_z.rstrip())

def trajectories_processing(atoms: np.ndarray):

    # t = [time.time()]
    # with open('dump_traj') as file:
    #     for lines, l in enumerate(file):
    #         print(lines)
    #     lines += 1
    # t.append(time.time())
    # print("Lines counted in {:.3f} seconds".format(t[1] - t[0]))

    # Initiates the trajectory
    for i in range(atoms.size):
        atoms[i].x_traj = np.zeros(atoms.size, dtype=np.float32)
        atoms[i].y_traj = np.zeros(atoms.size, dtype=np.float32)
        atoms[i].z_traj = np.zeros(atoms.size, dtype=np.float32)

    current_line = 1
    step = 0
    while True:

        # Loads data
        data_str = np.loadtxt('dump_traj', delimiter=',', dtype=str,
                              skiprows=(7 * (step + 1) + atoms.size * step),
                              max_rows=atoms.size)
        data = np.zeros((atoms.size, 3), dtype=np.float32)

        print("Loading data done")
        print(np.char.strip(data_str[:, 0], ''))
        print(np.char.strip(data_str[:, 1], ' '))
        print(np.char.strip(data_str[:, 2], ' }'))
        # Truncates data and converts it to floats
        for i in range(data_str.size // 3):
            data_str[i, 0] = data_str[i, 0][17:]
            data_str[i, 1] = data_str[i, 1][2:]
            data_str[i, 2] = data_str[i, 2][2:-1]
        data[:, 0] = data_str[:, 0].astype(np.float32)
        data[:, 1] = data_str[:, 1].astype(np.float32)
        data[:, 2] = data_str[:, 2].astype(np.float32)

        # Loads the positions in the Atoms
        for i in range(atoms.size):
            atoms[i].x_traj[step] = data[i, 0]
            atoms[i].y_traj[step] = data[i, 1]
            atoms[i].z_traj[step] = data[i, 2]

        step += 1
        t[2] = time.time()
        if step > 2:
            break
    pass
