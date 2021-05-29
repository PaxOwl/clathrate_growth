from ctypes import *
import numpy as np
import timeit
so_file = "./utils.so"
utils = cdll.LoadLibrary(so_file)

box = np.zeros(3)
distance = np.zeros(3)
for i in range(3):
    box[i] = 20
    distance[i] = i

p1 = np.zeros(3)
p2 = np.zeros(3)
p3 = np.zeros(3)

p1[0] = 2
p1[1] = 4
p1[2] = 6
p2[0] = 5
p2[1] = 8
p2[2] = 9
p3[0] = 4
p3[1] = 7
p3[2] = 1
period = True
v1 = np.zeros(3)
v2 = np.zeros(3)
utils.distance(c_void_p(p1.ctypes.data), c_void_p(p2.ctypes.data),
               c_void_p(box.ctypes.data), c_void_p(v1.ctypes.data))
utils.distance(c_void_p(p1.ctypes.data), c_void_p(p3.ctypes.data),
               c_void_p(box.ctypes.data), c_void_p(v2.ctypes.data))

print(v1)
print(v2)

utils.norm_vec(c_void_p(v1.ctypes.data))
utils.norm_vec(c_void_p(v2.ctypes.data))

print(np.sqrt(v1 ** 2))

theta = np.zeros(1)
utils.angle(c_void_p(v1.ctypes.data), c_void_p(v2.ctypes.data),
            c_void_p(theta.ctypes.data))
print(theta[0])
