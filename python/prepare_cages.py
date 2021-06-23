import numpy as np


INPUT = "data/10frames30-300k"
OUTPUT = INPUT + "-cages.dat"
SMALL = INPUT + "-small.dat"
LARGE = INPUT + "-large.dat"
INTER = INPUT + "-inter.dat"
IRREG = INPUT + "-irreg.dat"

width = np.loadtxt("data/60melt310-width.dat")
smallin = np.loadtxt(SMALL)
largein = np.loadtxt(LARGE)
interin = np.loadtxt(INTER)
irregin = np.loadtxt(IRREG)

cages = np.zeros((width.shape[0], 5))
cages[:, 0] = width[:, 0]
cages[:, 1] = smallin
cages[:, 2] = largein
cages[:, 3] = interin
cages[:, 4] = irregin

np.savetxt(OUTPUT, cages)
