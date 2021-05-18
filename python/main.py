#!/usr/bin/env python
"""
Core part of the program
"""
import time
from analysis import *
from parameters import *


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

    # Retrieves only the oxygen atoms
    oxygen = filter_data(atoms, ['OW'])

    # Initiates the array storing the AOP numbers
    aop = np.zeros((oxygen.shape[0], 2), dtype=np.float32)

    t1 = time.time()
    for i in range(oxygen.shape[0]):
        # Select an atom of oxygen
        center = oxygen.iloc[i]

        # Finds the nearest neighbours under a certain distance of a given atom
        neighbours = nearest_neighbours(oxygen, center, 0.35, box, periodic)

        # Computes the AOP for the selecter atom
        aop[i, 1] = compute_aop(center, neighbours)
    t2 = time.time()
    print("Elapsed time: {:.4f} s".format(t2 - t1))
    save_aop(aop, oxygen, periodic)
    # compute_rdf(mols, frame, met_rdf)
