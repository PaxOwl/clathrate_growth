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
    # Read the number of atoms (rows)
    nrows = count_atoms(filename)

    # Loads the atoms in a pandas DataFrame
    atoms = load_atoms(nrows, filename)

    # Load the data of the selected frame in the DataFrame
    load_frame(trimmed_data, atoms, frame, nrows)

    # Load the size of the box
    box = load_box(box_file, frame)

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
    cages_small = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
    cages_large = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
    cages_inter = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
    cages_irreg = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
    alone_met = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z'])
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
        bonds = pd.DataFrame(columns=['mol', 'atom', 'x', 'y', 'z', 'wat'])
        for jindex, j in low_aop.iterrows():
            hy1 = atoms.loc[(atoms.atom == 'HW1') &
                            (atoms.mol == low_aop.loc[jindex].mol)].squeeze()
            hy2 = atoms.loc[(atoms.atom == 'HW2') &
                            (atoms.mol == low_aop.loc[jindex].mol)].squeeze()
            bonds = bonds.append(hbonds(j, low_aop, hy1, hy2, box),
                                 ignore_index=True)
        sort_bonds = bonds.copy()
        for jindex, j in bonds.iterrows():
            neigh = bonds[bonds.wat == j.mol]
            if neigh.shape[0] == 0:
                sort_bonds = sort_bonds.drop(jindex)

        if low_aop.shape[0] == 20:
            cages_small = cages_small.append(center)
        elif low_aop.shape[0] == 24:
            cages_large = cages_large.append(center)
        elif 10 <= low_aop.shape[0] <= 13:
            cages_inter = cages_inter.append(center)
        elif (low_aop.shape[0] >= 14
              & low_aop.shape[0] != 20
              & low_aop.shape[0] != 24):
            cages_irreg = cages_irreg.append(center)
        else:
            alone_met = alone_met.append(center)

        print("mol = {}, x = {:.3f}"
              .format(carbon.iloc[index].mol, carbon.iloc[index].x))
        print("nw = {}, nh = {}, nb = {}"
              .format(neigh_ca.shape[0], low_aop.shape[0],
                      sort_bonds.shape[0]))
    size = clath_phase(cages_small, cages_large)
    print("Size of the clathrate phase: {} A".format(size))
    print("Total time {:.4f} s".format(time.time() - t_metinit))
    print('done')
