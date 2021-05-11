#!/usr/bin/env python
"""
Core part of the program
"""
from analysis import *
from parameters import *

nrows = count_atoms(filename)
atoms = load_atoms(nrows, filename)

met_rdf = {}
frame = 0

if __name__ == "__main__":
    load_frame(trimmed_data, atoms, frame, nrows)
    aop = compute_aop(atoms)
    save_aop(aop)
    # compute_rdf(mols, frame, met_rdf)
