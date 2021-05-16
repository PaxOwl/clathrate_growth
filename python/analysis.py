"""
This file sets-up and arranges the data to be used
"""
import time
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from MDAnalysis.lib.pkdtree import PeriodicKDTree
from parameters import aop_radius


def count_atoms(file: str = "conf.gro") -> int:
    """
    Reads the second line of the file which contains the number of atoms
    :param file: str, name of the file to read (conf.gro by default)
    :return: nrow: int, number of atoms in the file
    """
    nrow = np.loadtxt(file, skiprows=1, max_rows=1, usecols=0, dtype=int)

    return int(nrow)


def load_atoms(nrow: int, file: str = "conf.gro") -> pd.DataFrame:
    """
    Reads all the atoms and their coordinates in a given .gro file, then sorts
    them in the Atom dataclass
    :param nrow: int, number of atoms contained in the file
    :param file: str, name of the file to read (conf.gro by default)
    :return: a ndarray containing all the Atom dataclasses
    """

    data = pd.read_csv(file, sep=' ', nrows=nrow,
                       names=['mol', 'atom', 'x', 'y', 'z'],
                       usecols=[0, 1, 3, 4, 5], skiprows=[0, 1],
                       dtype={'mol': str, 'atom': str,
                              'x': np.float32,
                              'y': np.float32,
                              'z': np.float32})

    return data


def load_frame(file: str, atoms: pd.DataFrame, frame: int, nrow: int):

    # Loads data
    data = pd.read_csv(file, nrows=nrow, names=['x', 'y', 'z'],
                       dtype={'x': np.float32,
                              'y': np.float32,
                              'z': np.float32})

    # Replaces the old data with the positions from the loaded frame
    atoms.x = data.x
    atoms.y = data.y
    atoms.z = data.z

    return atoms


def load_box(file: str, frame: int) -> np.ndarray:
    return np.loadtxt(file, max_rows=3, skiprows=3 * frame)


def filter_data(data: pd.DataFrame, keep: list):

    lst = []
    for element in keep:
        lst.append(data[(data.atom == element)])
    output = pd.concat(lst)

    return output.sort_index()


def distance(p1: pd.Series, p2: pd.Series) -> tuple:
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    dz = p1.z - p2.z

    dst = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return dst, (dx, dy, dz)


def nearest_neighbours(data: pd.DataFrame, center: pd.Series,
                       radius: float, periodic: bool = False) -> pd.DataFrame:
    neighbours = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z',
                                       'x_dst', 'y_dst', 'z_dst'])
    for index, i in data.iterrows():
        if i.mol == center.mol:
            pass
        else:
            dst, details = distance(center, i)
            if dst <= radius:
                neighbours = neighbours.append(i)
                neighbours.loc[index, 'x_dst'] = details[0]
                neighbours.loc[index, 'y_dst'] = details[1]
                neighbours.loc[index, 'z_dst'] = details[2]
                print('appended')

    if periodic:
        pass

    return neighbours


def compute_aop(center: pd.DataFrame, neighbours: pd.DataFrame):

    # Computes the angles between the central oxygen and any other pair of
    # oxygen atoms within the range of the nearest neighbours
    angles = []
    while True:
        for i in range(len(neighbours) - 1):
            v1 = (neighbours.iloc[0].x_dst,
                  neighbours.iloc[0].y_dst,
                  neighbours.iloc[0].z_dst)
            v2 = (neighbours.iloc[i + 1].x_dst,
                  neighbours.iloc[i + 1].y_dst,
                  neighbours.iloc[i + 1].z_dst)
            theta = np.arccos(np.dot(v1, v2)
                              / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            angles.append(theta)

        neighbours = neighbours.drop(neighbours.index[0])

        if len(neighbours) <= 1:
            break

    # Computes the AOP for the oxygen
    aop = 0
    for i in range(len(angles)):
        aop += (np.abs(np.cos(angles[i]))
                * np.cos(angles[i])
                + np.cos(np.radians(109.47)) ** 2) ** 2

    return aop

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
    np.savetxt('aop.dat', output)
    np.savetxt('aop_moy.dat', output_moy)


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
            if (dst0 and dst1 and dst2 and dst3) < aop_radius:
                met_rdf[met.name].append(wat)
                met_rdf[met.name + "-dst"].append(dst1)
        t1 = time.time()
        print("{:0>2} neighbors in the vicinity of {},"
              " computed in {:.3f} seconds"
              .format(len(met_rdf[met.name]), met.name, t1 - t0))
    t_end = time.time()
    print("All neighbors computed in {:.3f} seconds".format(t_end - t_init))
    return met_rdf
