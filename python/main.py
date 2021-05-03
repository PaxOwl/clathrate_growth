#!/usr/bin/env python
"""
Core part of the program
"""
from analysis import *
from parameters import *

nrows = count_atoms(filename)
atoms = load_atoms(nrows, filename)
mols = compute_molecules(nrows, atoms, filename)


if __name__ == "__main__":
    trajectories_processing(atoms)
