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


def neighbours(centers: pd.DataFrame, neighbours: pd.DataFrame,
               box: np.ndarray, h_lim: float, l_lim: float):

    utils.neighbours.argtypes = [_2ddoublepp, _2ddoublepp,
                                 _1ddoublepp, ctypes.c_double, ctypes.c_double,
                                 ctypes.c_int, ctypes.c_int,
                                 _1ddoublepp, _2dlongpp]
    utils.neighbours.restype = None
    np_center = np.ascontiguousarray(np.delete(centers.to_numpy(),
                                               (0, 1), 1).astype(float))
    np_neighbours = np.ascontiguousarray(np.delete(neighbours.to_numpy(),
                                                   (0, 1), 1).astype(float))
    output = np.ascontiguousarray(np.zeros((centers.shape[0],
                                            neighbours.shape[0]))).astype(int)
    output.fill(-1)

    vec = np.ascontiguousarray(np.zeros(3, dtype=float))
    utils.neighbours(np_center, np_neighbours, box, h_lim, l_lim,
                     np_center.shape[0], np_neighbours.shape[0], vec, output)

    outdf = pd.DataFrame(output)
    outdf = outdf.replace([-1], np.nan)

    return outdf
