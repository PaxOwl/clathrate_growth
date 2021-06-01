#!/usr/bin/env python
"""
Core part of the program
"""
import ctypes
import time
import sys
from analysis import *
from parameters import *
from cintegration import *
from ctypes import *
from numpy.ctypeslib import ndpointer


so_file = "c-functions/utils.so"
utils = cdll.LoadLibrary(so_file)

frame = 0

if __name__ == "__main__":

    # hydrogen_bonds_test()
    # sys.exit("Done")

    # Read the number of atoms (rows)
    nrows = count_atoms(filename)

    # Loads the atoms in a pandas DataFrame
    atoms = load_atoms(nrows, filename)

    # Load the data of the selected frame in the DataFrame
    load_frame(trimmed_data, atoms, frame, nrows)

    # Load the size of the box
    box = load_box(box_file, frame)

    # Retrieves only the oxygen atoms
    oxygen = filter_data(atoms, ['OW'])
    carbon = filter_data(atoms, ['C'])
    # np_oxygen = np.ascontiguousarray(np.delete(oxygen.to_numpy(),
    #                                            (0, 1), 1).astype(float))
    # np_carbon = np.ascontiguousarray(np.delete(carbon.to_numpy(),
    #                                            (0, 1), 1).astype(float))
    # _output = np.ascontiguousarray(np.zeros((oxygen.shape[0],
    #                                          oxygen.shape[0]))).astype(int)
    # _output.fill(-1)
    # limit = 0.35
    # vec = np.ascontiguousarray(np.zeros(3, dtype=float))

    # _1ddoublepp = ndpointer(dtype=float, ndim=1, flags='C')
    # _2ddoublepp = ndpointer(dtype=float, ndim=2, flags='C')
    # _1dlongpp = ndpointer(dtype=np.int64, ndim=1, flags='C')
    # _2dlongpp = ndpointer(dtype=np.int64, ndim=2, flags='C')
    # utils.neighbours.argtypes = [_2ddoublepp, _2ddoublepp,
    #                              _1ddoublepp, ctypes.c_double,
    #                              ctypes.c_int, ctypes.c_int,
    #                              _1ddoublepp, _2dlongpp]
    # utils.neighbours.restype = None
    # utils.nearest_neighbours.argtypes = [_1ddoublepp, _2ddoublepp,
    #                                      _1ddoublepp, ctypes.c_double,
    #                                      ctypes.c_int,
    #                                      _1ddoublepp, _1dlongpp]
    # utils.nearest_neighbours.restype = None
    t1c = time.time()
    # utils.neighbours(np_oxygen, np_oxygen, box, limit,
    #                  np_oxygen.shape[0], np_oxygen.shape[0], vec, _output)

    # _outdf = pd.DataFrame(_output)
    outdf = neighbours(carbon, oxygen, box, 0.35, 0.0)
    # _outdf = _outdf[_outdf[1] != -1]
    # _outdf = _outdf.replace([-1], np.nan)
    # _outdf.dropna(how='all', axis=1, inplace=True)
    t2c = time.time()
    # utils.nearest_neighbours(np_oxygen[47], np_oxygen, box,
    #                          limit, oxygen.shape[0],
    #                          vec, _output[0])

    # sys.exit()
    # Initiates the array storing the AOP numbers
    aop = oxygen.copy()
    aop.loc[:, 'aop'] = 0.

    t1 = time.time()
    for index, i in outdf.iterrows():
        # Select an atom of oxygen
        center = oxygen.iloc[index]
        for j in i:
            pass
        print("oui")
        # Finds the nearest neighbours under a certain distance of a
        # given atom
        print("oui")

        # Computes the AOP for the selecter atom
        # aop.iat[i, 5] = compute_aop(center, neighbours)
    t2 = time.time()
    print("C : {:.4f} s".format(t2c - t1c))
    print("Python : {:.4f} s".format(t2 - t1))
    sys.exit()
    print("Elapsed time: {:.4f} s".format(t2 - t1))
    save_aop(aop.aop.values, oxygen, periodic)

    carbon = filter_data(atoms, ['C'])
    t_metinit = time.time()
    for i in range(carbon.shape[0]):
        t_met_un = time.time()
        # Select an atom of methane
        center = carbon.iloc[i]

        # Finds the nearest neighbours under a certain distance of a
        # given atom
        neighbours = nearest_neighbours(oxygen, center, 0.55, box, periodic)
        print("{} neighbours found".format(neighbours.shape[0]))
        low_aop = neighbours.copy()

        # Filter the neighbours to keep the one with an aop < 0.4
        for index, j in low_aop.iterrows():
            if aop.loc[aop.mol == j.mol].aop.values > 0.4:
                low_aop.drop(index)

        # Compute the hydrogen bonding of the returned molecules
        bonds = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z', 'ox_mol'])
        for index, j in low_aop.iterrows():
            hy1 = atoms.loc[(atoms.atom == 'HW1') &
                            (atoms.mol == low_aop.loc[index].mol)].squeeze()
            hy2 = atoms.loc[(atoms.atom == 'HW2') &
                            (atoms.mol == low_aop.loc[index].mol)].squeeze()
            bonds = bonds.append(hydrogen_bonds(j, oxygen, hy1, hy2, box))
        print("mol = {}, x = {}"
              .format(carbon.iloc[i].mol, carbon.iloc[i].x))
        print("nw = {}, nh = {}, nb = {}"
              .format(neighbours.shape[0], low_aop.shape[0], bonds.shape[0]))
        t_met_done = time.time()
        print("Elapsed time {}".format(t_met_done - t_met_un))
    print("Total time {}".format(time.time() - t_metinit))
    print('done')