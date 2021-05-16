from scipy import interpolate, signal
import numpy as np
import pandas as pd

def sample():
    data = pd.read_csv('aop.dat', sep=' ', names=['x', 'y'])

    linterp = interpolate.interp1d(data.x, data.y)
    xnew = np.arange(min(data.x), max(data.x), 0.001)

    ynew = linterp(xnew)

    sos = signal.butter(2, 0.5, 'lowpass', fs=1000, output='sos')

    filtered = signal.sosfilt(sos, ynew)

    new_sig = np.array((xnew, filtered), dtype=np.float).transpose()

    np.savetxt('rs_aop.dat', new_sig)

    return None
