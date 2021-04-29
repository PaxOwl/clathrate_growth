import numpy as np
from numpy.linalg import norm


theta = 107.8 * np.pi / 180
c = np.zeros(3)
h1 = np.copy(c)
h1[1] = 0.109

rx = np.zeros((3, 3))
rx[0, 0] = 1
rx[1, 1] = np.cos(theta)
rx[1, 2] = np.sin(theta)
rx[2, 1] = -np.sin(theta)
rx[2, 2] = np.cos(theta)

ry = np.zeros((3, 3))
ry[0, 0] = np.cos(theta)
ry[0, 2] = -np.sin(theta)
ry[1, 1] = 1
ry[2, 0] = np.sin(theta)
ry[2, 2] = np.cos(theta)

h2 = np.matmul(rx, h1)
h3 = np.matmul(rx, h2)
h4 = np.matmul(ry, h3)

print(np.arccos(np.dot(h2, h3) / (norm(h2) * norm(h3))) * 180 / np.pi)
print(norm(h1))
print(norm(h2))
print(norm(h3))
print(norm(h4))
