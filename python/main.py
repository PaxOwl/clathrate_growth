#!/usr/bin/env python
"""
Core part of the program
"""
import time
import sys
from analysis import *
from parameters import *


frame = 0

if __name__ == "__main__":
    hydrogen_bonds()
    sys.exit("Done")

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

    # Initiates the array storing the AOP numbers
    aop = oxygen.copy()
    aop.loc[:, 'aop'] = 0.

    t1 = time.time()
    for i in range(oxygen.shape[0] // 30):
        # Select an atom of oxygen
        center = oxygen.iloc[i]

        # Finds the nearest neighbours under a certain distance of a
        # given atom
        neighbours = nearest_neighbours(oxygen, center, 0.35,
                                        box, periodic)

        # Computes the AOP for the selecter atom
        aop.iat[i, 5] = compute_aop(center, neighbours)
    t2 = time.time()
    print("Elapsed time: {:.4f} s".format(t2 - t1))
    save_aop(aop.aop.values, oxygen, periodic)

    methane = filter_data(atoms, ['C'])
    for i in range(methane.shape[0]):
        # Select an atom of methane
        center = methane.iloc[i]

        # Finds the nearest neighbours under a certain distance of a
        # given atom
        neighbours = nearest_neighbours(oxygen, center, 0.55, box, periodic)
        print("{} neighbours found".format(neighbours.shape[0]))
        low_aop = neighbours.copy()
        for index, j in low_aop.iterrows():
            if aop[aop.mol == j.mol].aop > 0.4:
                low_aop.drop(index)
        print("done")
