"""
This file sets-up and arranges the data to be used
"""
import time
import numpy as np
from classes import Atom, Molecule
from parameters import radius


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
                      file: str = "conf.gro") -> dict:
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

    met_list = []
    wat_list = []
    mols_dic = {}

    counter = 0
    wat_count = 1
    met_count = 1

    while True:
        if "SOL" in mol[counter]:
            wat_list.append(Molecule("WAT-" + "{:0>4}".format(wat_count),
                                      counter // 4,
                                      (atoms_list[counter],
                                       atoms_list[counter + 1],
                                       atoms_list[counter + 2],
                                       atoms_list[counter + 3])))
            wat_count += 1
            counter += 4

        elif "CH4" in mol[counter]:
            met_list.append(Molecule("MET-" + "{:0>4}".format(met_count),
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

    mols_dic['MET'] = met_list
    mols_dic['WAT'] = wat_list

    return mols_dic


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

def load_frame(atoms: np.ndarray, frame):

    # Loads data
    data = np.loadtxt('trimmed_data', delimiter=',', dtype=np.float32,
                      skiprows=atoms.size * frame,
                      max_rows=atoms.size)

    # Loads the positions in the Atoms
    for i in range(atoms.size):
        atoms[i].x = data[i, 0]
        atoms[i].y = data[i, 1]
        atoms[i].z = data[i, 2]

    pass

def compute_rdf(mols: dict, frame: int, met_rdf: dict):

    t_init = time.time()
    for met in mols['MET']:
        t0 = time.time()
        met_rdf[met.name] = []
        met_rdf[met.name + "-dst"] = []
        for wat in mols['WAT']:
            dst0 = np.sqrt((wat.contains[0].x - met.contains[0].x) ** 2
                           + (wat.contains[0].y - met.contains[0].y) ** 2
                           + (wat.contains[0].z - met.contains[0].z) ** 2)
            dst1 = np.sqrt((wat.contains[1].x - met.contains[0].x) ** 2
                           + (wat.contains[1].y - met.contains[0].y) ** 2
                           + (wat.contains[1].z - met.contains[0].z) ** 2)
            dst2 = np.sqrt((wat.contains[2].x - met.contains[0].x) ** 2
                           + (wat.contains[2].y - met.contains[0].y) ** 2
                           + (wat.contains[2].z - met.contains[0].z) ** 2)
            dst3 = np.sqrt((wat.contains[3].x - met.contains[0].x) ** 2
                           + (wat.contains[3].y - met.contains[0].y) ** 2
                           + (wat.contains[3].z - met.contains[0].z) ** 2)
            if (dst0 and dst1 and dst2 and dst3) < radius:
                met_rdf[met.name].append(wat)
                met_rdf[met.name + "-dst"].append(dst1)
        t1 = time.time()
        print("{:0>2} neighbors in the vicinity of {},"
              " computed in {:.3f} seconds"
              .format(len(met_rdf[met.name]), met.name, t1 - t0))
    t_end = time.time()
    print("All neighbors computed in {:.3f} seconds".format(t_end - t_init))
    return met_rdf
