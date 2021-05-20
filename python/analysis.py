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


def distance(p1: pd.Series, p2: pd.Series, box: np.ndarray,
             periodic: bool) -> tuple:
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    dz = p1.z - p2.z

    if periodic:
        dx = dx - int(round(dx / box[0])) * box[0]
        dy = dy - int(round(dy / box[1])) * box[1]
        dz = dz - int(round(dz / box[2])) * box[2]
    dst = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return dst, (dx, dy, dz)


def angle(p1: pd.Series, p2: pd.Series, p3: pd.Series) -> np.float64:
    v1 = (p1.x - p2.x, p1.y - p2.y, p1.z, p2.z)
    v1 /= np.linalg.norm(v1)
    v2 = (p1.x - p3.x, p1.y - p3.y, p1.z, p3.z)
    v2 /= np.linalg.norm(v2)
    theta = np.arccos(np.dot(v1, v2))

    return np.degrees(theta)


def nearest_neighbours(data: pd.DataFrame, center: pd.Series,
                       radius: float, box: np.ndarray,
                       periodic: bool) -> pd.DataFrame:
    neighbours = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z',
                                       'x_dst', 'y_dst', 'z_dst'])
    for index, i in data.iterrows():
        if i.mol == center.mol:
            pass
        else:
            if periodic:
                dst, details = distance(center, i, box, True)
            else:
                dst, details = distance(center, i, box, False)
            if dst <= radius:
                neighbours = neighbours.append(i)
                neighbours.loc[index, 'x_dst'] = details[0]
                neighbours.loc[index, 'y_dst'] = details[1]
                neighbours.loc[index, 'z_dst'] = details[2]

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

    print("AOP computed for molecule {}".format(center.mol))

    return aop

def save_aop(aop: np.ndarray, oxygen: pd.DataFrame, periodic: bool):
    output = np.zeros((aop.shape[0], 2))
    output[:, 1] = aop
    for i in range(aop.shape[0]):
        output[i, 0] = oxygen.iloc[i].x

    output = output[output[:, 0].argsort()]

    if periodic:
        np.savetxt('aop_periodic.dat', output)
    else:
        np.savetxt('aop.dat', output)

def hydrogen_bonds():
    data = pd.read_csv('test_hbonds/out', sep=',',
                       usecols=[0, 1, 3, 4, 5],
                       names=['mol', 'atom', 'x', 'y', 'z'])
    box = np.array([3.12913551, 2.94142906, 3.61460741], dtype=np.double)
    oxygen = filter_data(data, ['OW'])
    hydrogen1 = filter_data(data, ['HW1'])
    hydrogen2 = filter_data(data, ['HW2'])
    pairs = []
    for index, i in oxygen.iterrows():
        center = i
        for jdex, j in oxygen.iterrows():
            if center.mol == j.mol:
                pass
            else:
                dst, _ = distance(center, j, box, True)
                if 0.25 <= dst <= 0.35:
                    # Finds the closer hydrogen
                    d1, _ = distance(center,
                                  hydrogen1.loc[hydrogen1.mol
                                                == center.mol].squeeze(),
                                  box, False)
                    d2, _ = distance(j,
                                  hydrogen1.loc[hydrogen1.mol
                                                == center.mol].squeeze(),
                                  box, False)
                    d3 = d1 + d2
                    d4, _ = distance(center,
                                  hydrogen2.loc[hydrogen2.mol
                                                == center.mol].squeeze(),
                                  box, False)
                    d5, _ = distance(j,
                                  hydrogen2.loc[hydrogen2.mol
                                                == center.mol].squeeze(),
                                  box, False)
                    d6 = d4 + d5

                    if d3 < d6:
                        hydrogen_atom = hydrogen1.loc[hydrogen1.mol
                                                      == center.mol].squeeze()
                    else:
                        hydrogen_atom = hydrogen2.loc[hydrogen2.mol
                                                      == center.mol].squeeze()
                    # Compute the angle
                    theta = angle(hydrogen_atom, center, j)
                    print(theta)
                    if 90 < theta < 180:
                        print("{} accepts from {}".format(center.mol, j.mol))
                        pairs.append((center.mol, j.mol))
    with open("pairs.dat", 'w') as file:
        for element in pairs:
            file.write("{}\n".format(element))
