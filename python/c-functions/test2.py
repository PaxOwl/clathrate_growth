from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

so_file = "./test.so"
utils = cdll.LoadLibrary(so_file)

doublepp = ndpointer(dtype=float, ndim=2, flags='C')
size = 4
matrix = np.zeros((size, size), dtype=float)

utils.test.argtypes = [c_int, doublepp]
utils.test.restype = None

for i in range(size):
    for j in range(size):
        matrix[i, j] = i + j
matrix = np.ascontiguousarray(matrix)

utils.test(size, matrix)
