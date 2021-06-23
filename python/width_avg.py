import numpy as np

def width_average():
    width = np.loadtxt("data/growth-width.dat")

    linreg = np.zeros((width.shape[0], 2))
    p = np.polyfit(width[:, 0], width[:, 1], 1)
    linreg[:, 0] = width[:, 0]
    linreg[:, 1] = width[:, 0] * p[0] + p[1]

    np.savetxt("data/growth-linreg.dat", linreg)
    print(p)
    print((linreg[width.shape[0] - 1, 1] - linreg[0, 1]) / 40.)
