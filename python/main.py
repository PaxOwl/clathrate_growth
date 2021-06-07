#!/usr/bin/env python
"""
Core part of the program
"""
import time
import sys
from analysis import *
from parameters import *
from cintegration import *


frame = 0

if __name__ == "__main__":
    # TEST HBONDS
    # atoms = pd.read_csv('test_hbonds/out', sep=',',
    #                    usecols=[0, 1, 3, 4, 5],
    #                    names=['mol', 'atom', 'x', 'y', 'z'])
    # box = np.array([3.12913551, 2.94142906, 3.61460741])
    # oxygen = filter_data(atoms, ['OW'])
    # hydrogen1 = filter_data(atoms, ['HW1'])
    # hydrogen2 = filter_data(atoms, ['HW2'])
    # count = 0
    # for index, i in oxygen.iterrows():
    #     hy1 = atoms.loc[(atoms.atom == 'HW1') & (atoms.mol == i.mol)].squeeze()
    #     hy2 = atoms.loc[(atoms.atom == 'HW2') & (atoms.mol == i.mol)].squeeze()
    #     count += hbonds(i, oxygen, hy1, hy2, box)
    # print("\n", count)
    # sys.exit()
    # END OF TEST HBONDS

    # Read the number of atoms (rows)
    nrows = count_atoms(filename)

    # Loads the atoms in a pandas DataFrame
    atoms = load_atoms(nrows, filename)

    # Load the data of the selected frame in the DataFrame
    # load_frame(trimmed_data, atoms, frame, nrows)

    # Load the size of the box
    # box = load_box(box_file, frame)
    box = np.array([1.20034886, 1.20034886, 1.20034886])

    # Retrieves the oxygen and carbon atoms
    oxygen = filter_data(atoms, ['OW'])
    carbon = filter_data(atoms, ['C'])

    # Finds the nearest neighbours under a certain distance of a
    # given atom
    ox_neigh_ids = neighbours(oxygen, oxygen, box, 0.0, 0.35)

    # Initiates the array storing the AOP numbers
    aop_values = oxygen.copy()
    aop_values.loc[:, 'aop'] = 0.

    # Computes and stores the aop of the oxygen atoms
    t1 = time.time()
    for index, i in ox_neigh_ids.iterrows():
        neigh_ox = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])

        # Select an atom of oxygen
        center = oxygen.iloc[index]
        for j in i:
            if np.isnan(j):
                break
            neigh_ox = neigh_ox.append(oxygen.iloc[int(j)])

        # Computes the AOP for the selecter atom
        aop_values.iat[index, 5] = caop(center, neigh_ox, box)
        print(index)
    t2 = time.time()
    print("Elapsed time: {:.4f} s".format(t2 - t1))
    save_aop(aop_values.aop.values, oxygen, periodic)

    ca_neigh_ids = neighbours(carbon, oxygen, box, 0., 0.55)
    t_metinit = time.time()
    for index, i in ca_neigh_ids.iterrows():
        neigh_ca = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
        # Select an atom of carbon
        center = carbon.iloc[index]

        for j in i:
            if np.isnan(j):
                break
            neigh_ca = neigh_ca.append(oxygen.iloc[int(j)])

        low_aop = neigh_ca.copy()
        # Filter the neighbours to keep the one with an aop < 0.4
        for jindex, j in low_aop.iterrows():
            if aop_values.loc[aop_values.mol == j.mol].aop.values[0] > 0.4:
                low_aop = low_aop.drop(jindex)

        # Compute the hydrogen bonding of the returned molecules
        bonds = 0
        for jindex, j in low_aop.iterrows():
            hy1 = atoms.loc[(atoms.atom == 'HW1') &
                            (atoms.mol == low_aop.loc[jindex].mol)].squeeze()
            hy2 = atoms.loc[(atoms.atom == 'HW2') &
                            (atoms.mol == low_aop.loc[jindex].mol)].squeeze()
            bonds += hbonds(j, oxygen, hy1, hy2, box)
        print("mol = {}, x = {:.3f}"
              .format(carbon.iloc[index].mol, carbon.iloc[index].x))
        print("nw = {}, nh = {}, nb = {}"
              .format(neigh_ca.shape[0], low_aop.shape[0], bonds))
    print("Total time {:.4f} s".format(time.time() - t_metinit))
    print('done')
