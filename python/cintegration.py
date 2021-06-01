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


def neighbours(df_centers: pd.DataFrame, df_neigh: pd.DataFrame,
               box: np.ndarray, l_lim: float, h_lim: float):

    utils.neighbours.argtypes = [_2ddoublepp, _2ddoublepp,
                                 _1ddoublepp, ctypes.c_double, ctypes.c_double,
                                 ctypes.c_int, ctypes.c_int,
                                 _1ddoublepp, _2dlongpp]
    utils.neighbours.restype = None
    np_center = np.ascontiguousarray(np.delete(df_centers.to_numpy(),
                                               (0, 1), 1).astype(float))
    np_neigh = np.ascontiguousarray(np.delete(df_neigh.to_numpy(),
                                                   (0, 1), 1).astype(float))
    output = np.ascontiguousarray(np.zeros((df_centers.shape[0],
                                            df_neigh.shape[0]))).astype(int)
    output.fill(-1)

    vec = np.ascontiguousarray(np.zeros(3, dtype=float))
    utils.neighbours(np_center, np_neigh, box, h_lim, l_lim,
                     np_center.shape[0], np_neigh.shape[0], vec, output)

    outdf = pd.DataFrame(output)
    outdf = outdf[outdf[1] != -1]
    outdf = outdf.replace([-1], np.nan)
    outdf.dropna(how='all', axis=1, inplace=True)
    for index, i in outdf.iterrows():
        for j in range(len(i)):
            if np.isnan(i[j]):
                break
            print(int(i[j]))
    return outdf

def aop(center, neigh):

    return None
