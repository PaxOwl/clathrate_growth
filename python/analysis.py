"""
This file sets-up and arranges the data to be used
"""
import time
import numpy as np
import pandas as pd
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

def periodic_conditions(d: list, box: np.ndarray):
    d[0] = d[0] - int(round(d[0] / box[0])) * box[0]
    d[1] = d[1] - int(round(d[1] / box[1])) * box[1]
    d[2] = d[2] - int(round(d[2] / box[2])) * box[2]

    return d

def distance(p1: pd.Series, p2: pd.Series, box: np.ndarray,
             periodic: bool) -> tuple:
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    dz = p2.z - p1.z

    if periodic:
        dx, dy, dz = periodic_conditions([dx, dy, dz], box)
    dst = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return dst, (dx, dy, dz)


def norm_vec(p1: pd.Series, p2: pd.Series, box) -> tuple:
    _, v = distance(p1, p2, box, True)
    v /= np.linalg.norm(v)

    return v

def angle(p1: pd.Series, p2: pd.Series, p3: pd.Series, box) -> np.float64:
    v1 = norm_vec(p1, p2, box)
    v2 = norm_vec(p1, p3, box)
    theta = np.arccos(np.dot(v1, v2))

    return np.degrees(theta)


def closest_atom(ox1: pd.Series, ox2: pd.Series,
                 hy1: pd.Series, hy2: pd.Series, box):
    d1 = distance(ox1, hy1, box, True)[0]\
        + distance(ox2, hy1, box, True)[0]
    d2 = distance(ox1, hy2, box, True)[0]\
        + distance(ox2, hy2, box, True)[0]

    if d1 > d2:
        return hy2, d2
    else:
        return hy1, d1


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


def hydrogen_bonds(center: pd.Series, oxygen: pd.DataFrame,
                   hydrogen1: pd.Series, hydrogen2: pd.Series,
                   box: np.ndarray):
    bonds = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z', 'ox_mol'])
    # i is the oxygen atom of the neigbour molecule
    for index, i in oxygen.iterrows():
        if center.mol == i.mol:
            pass
        else:
            dst, _ = distance(center, i, box, True)
            if 0.25 <= dst <= 0.35:
                # Finds the closest hydrogen
                hydrogen_atom = closest_atom(center, i,
                                             hydrogen1, hydrogen2, box)
                # Compute the angle
                theta = angle(hydrogen_atom[0], center, i, box)
                if 90 < theta < 180:
                    bonds = bonds.append(hydrogen_atom[0],
                                         ignore_index=True)
                    bonds.loc[bonds.shape[0] - 1, 'ox_mol'] = i.mol

    return bonds


def hydrogen_bonds_test():
    data = pd.read_csv('test_hbonds/out', sep=',',
                       usecols=[0, 1, 3, 4, 5],
                       names=['mol', 'atom', 'x', 'y', 'z'])
    box = np.array([3.12913551, 2.94142906, 3.61460741])
    oxygen = filter_data(data, ['OW'])
    hydrogen1 = filter_data(data, ['HW1'])
    hydrogen2 = filter_data(data, ['HW2'])
    bonds = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z', 'ox_mol'])
    # i is the oxygen atom of the considered molecule
    for index, i in oxygen.iterrows():
        # j is the oxygen atom of the neigbour molecule
        for jdex, j in oxygen.iterrows():
            if i.mol == j.mol:
                pass
            else:
                dst, _ = distance(i, j, box, True)
                if 0.25 <= dst <= 0.35:
                    # Finds the closest hydrogen
                    hy1 = hydrogen1.loc[hydrogen1.mol == i.mol].squeeze()
                    hy2 = hydrogen2.loc[hydrogen2.mol == i.mol].squeeze()
                    hydrogen_atom = closest_atom(i, j, hy1, hy2, box)
                    # Compute the angle
                    theta = angle(hydrogen_atom[0], i, j, box)
                    if 90 < theta < 180:
                        bonds = bonds.append(hydrogen_atom[0],
                                             ignore_index=True)
                        bonds.loc[bonds.shape[0] - 1, 'ox_mol'] = j.mol
                        print("Bond added between {} and {}"
                              .format(i.mol, j.mol))
        if index >= 50:
            break
