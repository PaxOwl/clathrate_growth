from ctypes import *
import numpy as np
so_file = "./utils.so"
utils = cdll.LoadLibrary(so_file)

box = np.zeros(3)
distance = np.zeros(3)
for i in range(3):
    box[i] = 3
    distance[i] = i

output = utils.periodic_conditions(c_void_p(distance.ctypes.data),
                                   c_void_p(box.ctypes.data))
print(output)
