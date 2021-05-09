"""
This file sets-up and arranges the data to be used
"""
import time
import numpy as np
from scipy.spatial import cKDTree
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
        atoms_list.append(Atom(names[i], int(ids[i]) - 1, x[i], y[i], z[i]))

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
    wat_count = 0
    met_count = 0

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

    mols_dic['WAT'] = wat_list
    mols_dic['MET'] = met_list

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

def compute_aop(atoms: np.ndarray):
    oxygen = []
    neighbors = []
    atoms_aop = []

    # Counts the oxygen atoms
    for a in atoms:
        if 'OW' in a.name:
            oxygen.append(a)

    # Iterates on all the oxygen atoms
    for ox in oxygen:
        center = [ox.x, ox.y, ox.z]
        others = []

        # Stores all the atoms except the center
        for outer in oxygen:
            if outer.id != ox.id:
                others.append([outer.x, outer.y, outer.z])

        # Computes and stores its closest neighbours (r <= 0.35 nm)
        tree = cKDTree(others)
        count = tree.query_ball_point(center, r=0.35)
        positions = []
        for i in count:
            positions.append(others[i])
        neighbors.append((center, positions))
        print("Added neighbours for oxygen {}".format(ox.id))

        # Computes the angles
        angles = []
        while True:
            for i in range(len(positions) - 1):
                v1 = (center[0] - positions[0][0],
                      center[1] - positions[0][1],
                      center[2] - positions[0][2])
                v2 = (center[0] - positions[i + 1][0],
                      center[1] - positions[i + 1][1],
                      center[2] - positions[i + 1][2])
                theta = np.arccos(np.dot(v1, v2)
                                  / (np.linalg.norm(v1) * np.linalg.norm(v2)))
                angles.append(theta)

            positions = positions[1:]

            if len(positions) <= 1:
                break
        print("Computed angles for oxygen {}".format(ox.id))

        # Computes the AOP for the oxygen
        aop = 0
        for i in range(len(angles)):
            aop += (np.abs(np.cos(angles[i]))
                    * np.cos(angles[i])
                    + np.cos(np.radians(109.47)) ** 2) ** 2

        print("Computed aop for oxygen {}".format(ox.id))
        atoms_aop.append((ox, aop))

    return atoms_aop

def save_aop(atoms_aop: list):
    output = []
    output_moy = []
    for i in atoms_aop:
        output.append((i[0].x, i[1]))

    output_moy.append(((atoms_aop[0][0].x + atoms_aop[1][0].x) / 2,
                      atoms_aop[0][1]))
    for i in range(len(atoms_aop) - 1):
        output_moy.append(((atoms_aop[i - 1][0].x
                            + atoms_aop[i][0].x
                            + atoms_aop[i + 1][0].x) / 3, atoms_aop[i][1]))
    output_moy.append(((atoms_aop[len(atoms_aop) - 2][0].x
                       + atoms_aop[len(atoms_aop) - 1][0].x) / 2,
                      atoms_aop[len(atoms_aop) - 1][1]))

    dtype = [('x', float), ('aop', float)]
    output = np.array(output, dtype=dtype)
    output = np.sort(output, order='x')
    output_moy = np.array(output_moy, dtype=dtype)
    output_moy = np.sort(output_moy, order='x')
    np.savetxt('aop.dat', output_moy)


def compute_rdf(mols: dict, met_rdf: dict):

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
