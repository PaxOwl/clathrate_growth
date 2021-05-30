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

names = np.zeros(2)
names[0] = 6
names[1] = 19
oxygens = np.zeros((2, 3))
oxygens[0, 0] = 1
oxygens[0, 1] = 2
oxygens[0, 2] = 3
oxygens[1, 0] = 1
oxygens[1, 1] = 4
oxygens[1, 2] = 5

hydrogens = np.zeros((2, 3))
hydrogens[0, 0] = 1
hydrogens[0, 1] = 2
hydrogens[0, 2] = 3
hydrogens[1, 0] = 1
hydrogens[1, 1] = 4
hydrogens[1, 2] = 5

name = np.zeros(1)
print(name)
utils.closest_atom(c_void_p(names.ctypes.data), c_void_p(oxygens.ctypes.data),
                   c_void_p(hydrogens.ctypes.data), c_void_p(box.ctypes.data),
                   c_void_p(name.ctypes.data))
print(name)
