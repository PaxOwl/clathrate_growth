from ctypes import *
import numpy as np
import timeit
so_file = "./utils.so"
utils = cdll.LoadLibrary(so_file)

box = np.zeros(3)
distance = np.zeros(3)
for i in range(3):
    box[i] = 3
    distance[i] = i

print(timeit.timeit('utils.periodic_conditions(c_void_p(distance.ctypes.data),'
                                              'c_void_p(box.ctypes.data),'
                                              'c_void_p(distance.ctypes.data))',
                    'from test import utils, distance, box, c_void_p',
                    number=1000000))


def periodic_conditions(d: np.ndarray, box: np.ndarray):
    for i in range(len(d)):
        d[i] = d[i] - int(round(d[i] / box[i])) * box[i]
    return d


print(timeit.timeit('periodic_conditions(distance, box)', 'from test import periodic_conditions, distance, box', number=1000000))
