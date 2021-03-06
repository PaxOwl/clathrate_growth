import ctypes
from ctypes import *
import pandas as pd
import numpy as np
from numpy.ctypeslib import ndpointer


# Imports the shared library from C
so_file = "c-functions/utils.so"
utils = cdll.LoadLibrary(so_file)

# Sets some types to use with the C library
_1ddoublepp = ndpointer(dtype=float, ndim=1, flags='C')
_2ddoublepp = ndpointer(dtype=float, ndim=2, flags='C')
_1dlongpp = ndpointer(dtype=np.int64, ndim=1, flags='C')
_2dlongpp = ndpointer(dtype=np.int64, ndim=2, flags='C')

# Sets argtypes and restype for used C functions
utils.neighbours.argtypes = [_2ddoublepp, _2ddoublepp,
                             _1ddoublepp, ctypes.c_double, ctypes.c_double,
                             ctypes.c_int, ctypes.c_int,
                             _1ddoublepp, _2dlongpp]
utils.neighbours.restype = None

utils.distance.argtypes = [_1ddoublepp, _1ddoublepp, _1ddoublepp, _1ddoublepp]
utils.distance.restype = None
utils.angle.argtypes = [_1ddoublepp, _1ddoublepp, _1ddoublepp]
utils.angle.restype = None
utils.aop.argtypes = [_1ddoublepp, _2ddoublepp, _1ddoublepp,
                      ctypes.c_int, _1ddoublepp, _1ddoublepp,
                      _1ddoublepp, _1ddoublepp, _1ddoublepp]
utils.aop.restype = None
utils.hydrogen_bonds.argtypes = [_1ddoublepp, _2ddoublepp, _2ddoublepp,
                                 _1ddoublepp, _1ddoublepp, _1ddoublepp,
                                 _1ddoublepp, ctypes.c_int, _1ddoublepp,
                                 _1ddoublepp, _1ddoublepp]
utils.hydrogen_bonds.restype = None
utils.clath_size.argtypes = [_1ddoublepp, _1ddoublepp, ctypes.c_int,
                             ctypes.c_int, _1ddoublepp, _1ddoublepp]
utils.clath_size.restype = None


def neighbours(df_centers: pd.DataFrame, df_neigh: pd.DataFrame,
               box: np.ndarray, l_lim: float, h_lim: float):
    """
    Search the nearest neighbours in a sphere/circle of width h_lim - l_lim
    around each centers
    :param df_centers: DataFrame, atoms to which the neighbours are to be
    searched
    :param df_neigh: DataFrame, all potential neighbours
    :param box: ndarray, size of the box
    :param l_lim: float, lower bound of the search
    :param h_lim: float, upper bound of the search
    :return: a DataFrame containing the neighbours for each center
    """

    # Conversions from DataFrame to ndarray for use in C code
    np_center = np.ascontiguousarray(np.delete(df_centers.to_numpy(),
                                               (0, 1), 1).astype(float))
    np_neigh = np.ascontiguousarray(np.delete(df_neigh.to_numpy(),
                                              (0, 1), 1).astype(float))
    output = np.ascontiguousarray(np.zeros((df_centers.shape[0],
                                            df_neigh.shape[0]))).astype(int)
    dst = np.ascontiguousarray(np.zeros((df_centers.shape[0],
                                         df_neigh.shape[0])))
    output.fill(-1)
    dst.fill(-1)
    vec = np.ascontiguousarray(np.zeros(3, dtype=float))

    # Call of the C code
    utils.neighbours(np_center, np_neigh, box, h_lim, l_lim,
                     np_center.shape[0], np_neigh.shape[0], vec, output)

    # Formatting the data for use in DataFrames
    outdf = pd.DataFrame(output)
    outdf = outdf[outdf[1] != -1]
    outdf = outdf.replace([-1], np.nan)
    outdf.dropna(how='all', axis=1, inplace=True)

    return outdf


def caop(df_center, neigh, box):
    """
    Computes the AOP for a molecule of water
    :param df_center: Dataframe, atom to chich the aop must be computed
    :param neigh: DataFrame, contains the neighbours
    :param box: ndarray, size of the box
    :return: the AOP of the center
    """

    # Conversions from DataFrame to ndarray for use in C code
    vec1 = np.ascontiguousarray(np.zeros(3, dtype=float))
    vec2 = np.ascontiguousarray(np.zeros(3, dtype=float))
    np_center = np.ascontiguousarray(np.array([df_center.x,
                                               df_center.y,
                                               df_center.z]), dtype=float)
    np_neigh = np.ascontiguousarray(np.delete(neigh.to_numpy(),
                                              (0, 1), 1)).astype(float)
    angles = np.zeros((neigh.shape[0] - 1) * neigh.shape[0] // 2,
                      dtype=float)
    angles = np.ascontiguousarray(angles)
    angles.fill(-1)
    theta = np.ascontiguousarray(np.zeros(1, dtype=float))
    aop = np.ascontiguousarray(np.zeros(1, dtype=float))

    # Call of the C code
    utils.aop(np_center, np_neigh, box, np_neigh.shape[0],
              vec1, vec2, theta, angles, aop)

    return aop[0]


def hbonds(sr_center: pd.Series, df_oxygens: pd.DataFrame,
           sr_hydrogen1: pd.Series, sr_hydrogen2: pd.Series,
           box: np.ndarray):
    """
    Computes the hydrogen bonds for a given water molecule
    :param sr_center: Series, atom to which the hydrogen bounds must be found
    :param df_oxygens: DataFrame, all the oxygens
    :param sr_hydrogen1: DataFrame, first hydrogen of each water molecule
    :param sr_hydrogen2: DataFrame, second hydrogen of each water molecule
    :param box: ndarray, size of the box
    :return: a DataFrame containing the molecules and their hydrogen bonds
    """

    # Conversions from DataFrame to ndarray for use in C code
    np_center = np.ascontiguousarray(np.array([sr_center.x,
                                               sr_center.y,
                                               sr_center.z]), dtype=float)
    np_oxygens = np.ascontiguousarray(np.delete(df_oxygens.to_numpy(),
                                                (0, 1), 1)).astype(float)
    np_hydrogens = np.ascontiguousarray(np.array(([sr_hydrogen1.x,
                                                   sr_hydrogen1.y,
                                                   sr_hydrogen1.z],
                                                  [sr_hydrogen2.x,
                                                   sr_hydrogen2.y,
                                                   sr_hydrogen2.z])),
                                        dtype=float)
    vec1 = np.ascontiguousarray(np.zeros(3, dtype=float))
    vec2 = np.ascontiguousarray(np.zeros(3, dtype=float))
    vec3 = np.ascontiguousarray(np.zeros(3, dtype=float))
    closest = np.ascontiguousarray(np.zeros(3, dtype=float))
    theta = np.ascontiguousarray(np.zeros(1, dtype=float))
    bonds = np.ascontiguousarray(np.zeros(np_oxygens.shape[0], dtype=float))

    # Call of the C code
    utils.hydrogen_bonds(np_center, np_oxygens, np_hydrogens, box,
                         vec1, vec2, vec3, np_oxygens.shape[0],
                         closest, theta, bonds)

    # Formatting the data for use in DataFrames
    df_bonds = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z', 'wat'])
    for index, i in enumerate(bonds):
        if i == 1:
            df_bonds = df_bonds.append(df_oxygens.iloc[index])
    df_bonds.loc[:, 'wat'] = sr_center.mol

    return df_bonds


def clath_phase(df_small: pd.DataFrame, df_large: pd.DataFrame,
                box: np.ndarray):
    """
    Computes the size of the clathrate phase
    :param df_small: DataFrame, contains the molecules of the small cages
    :param df_large: DataFrame, contains the molecules of the large cages
    :param box: ndarray, size of the box
    :return: the size of the clathrate phase
    """

    # Conversions from DataFrame to ndarray for use in C code
    xmin = np.ascontiguousarray(np.zeros(1, dtype=float))
    xmin[0] = box[0]
    xmax = np.ascontiguousarray(np.zeros(1, dtype=float))
    np_small = np.ascontiguousarray(np.delete(df_small.to_numpy(),
                                              (0, 1, 3, 4), 1)).astype(float)
    np_large = np.ascontiguousarray(np.delete(df_large.to_numpy(),
                                              (0, 1, 3, 4), 1)).astype(float)
    np_small = np_small.reshape(df_small.shape[0])
    np_large = np_large.reshape(df_large.shape[0])

    # Call of the C code
    utils.clath_size(np_small, np_large,
                     np_small.shape[0], np_large.shape[0],
                     xmin, xmax)

    return xmax[0] - xmin[0]
