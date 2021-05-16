#!/usr/bin/env python
"""
Core part of the program
"""
from analysis import *
from parameters import *
from resampling import sample
import sys


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

    # Select an atom of oxygen
    center = oxygen[oxygen.mol == '1SOL'].loc[0]

    # Finds the nearest neighbours under a certain distance of a given atom
    neighbours = nearest_neighbours(oxygen, center, 0.35)
    aop = compute_aop(center, neighbours)
    save_aop(aop)
    # compute_rdf(mols, frame, met_rdf)
