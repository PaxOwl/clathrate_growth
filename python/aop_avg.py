import numpy as np


def aop_compute_avg():
    aop = np.loadtxt("data/100melt310-aop.dat")
    avg = np.zeros((100, 2), dtype=float)
    for i in range(100):
        tmp = []
        for yindex, y in enumerate(aop[:, 0]):
            if y > i / 10.:
                break
            else:
                tmp.append(aop[yindex, 1])
        avg[i, 1] = np.average(tmp)
        avg[i, 0] = i / 10.
        avg[0, 1] = 0

    np.savetxt("data/TESTAOP.DAT", avg)
