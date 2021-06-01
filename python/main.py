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

    # Finds the nearest neighbours under a certain distance of a
    # given atom
    ox_neigh_ids = neighbours(oxygen, oxygen, box, 0.0, 0.35)

    # Initiates the array storing the AOP numbers
    aop = oxygen.copy()
    aop.loc[:, 'aop'] = 0.

    t1 = time.time()
    for index, i in ox_neigh_ids.iterrows():

        neigh = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
        # Select an atom of oxygen
        center = oxygen.iloc[index]
        for j in i:
            if np.isnan(j):
                break
            neigh = neigh.append(oxygen.iloc[int(j)])


        # Computes the AOP for the selecter atom
        # aop.iat[i, 5] = compute_aop(center, neighbours)
    t2 = time.time()
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