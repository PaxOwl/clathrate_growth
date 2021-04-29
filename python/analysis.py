"""
This file sets-up and arranges the data to be used
"""

import numpy as np
from classes import Atom, Molecule
from parameters import *


def count_rows(file):
    with open(file, "r") as gro:
        for i, l in enumerate(gro):
            pass
    nrow = i - 2

    return nrow

def load_atoms(file, nrow):

    names = np.loadtxt(file, skiprows=2, usecols=1, max_rows=nrow,
                       dtype=str)
    ids = np.loadtxt(file, skiprows=2, usecols=2, max_rows=nrow,
                     dtype=str)
    x = np.loadtxt(file, skiprows=2, usecols=3, max_rows=nrow)
    y = np.loadtxt(file, skiprows=2, usecols=4, max_rows=nrow)
    z = np.loadtxt(file, skiprows=2, usecols=5, max_rows=nrow)

    # atoms = np.empty(nrow, dtype=Atom)
    atoms_list = []

    for i in range(nrow):
        atoms_list.append(Atom(names[i], ids[i], x[i], y[i], z[i]))

    return atoms_list


def compute_molecules(file, nrow, atoms_list):

    mol = np.loadtxt(file, skiprows=2, usecols=0, max_rows=nrow,
                     dtype=str)

    nmols = nrow // 4
    mols_list = []

    counter = 0

    while True:
        if mol[counter] == "1SOL":
            mols_list.append(Molecule("WAT",
                                      counter // 4,
                                      (atoms_list[counter],
                                       atoms_list[counter + 1],
                                       atoms_list[counter + 2],
                                       atoms_list[counter + 3])))
        elif mol[counter] == "1CH4":
            mols_list.append(Molecule("MET",
                                      counter // 4,
                                      (atoms_list[counter],
                                       atoms_list[counter + 1],
                                       atoms_list[counter + 2],
                                       atoms_list[counter + 3])))
        counter += 4

        if (counter + 3) // 4 >= nmols:
            break

    return mols_list

def print_atom(atom):
    print("name: {}\n"
          "id: {}\n"
          "coordinates:\n"
          "    x: {}\n"
          "    y: {}\n"
          "    z: {}\n".format(atom.name, atom.id, atom.x, atom.y, atom.z))

def print_mol(mol):
    atoms_name = ""
    atoms_id = ""
    atoms_x = ""
    atoms_y = ""
    atoms_z = ""

    for i in range(4):
        atoms_name += "  atom name: {:<5}".format(mol.contains[i].name)
        atoms_id += "  atom id: {:<7}".format(mol.contains[i].id)
        atoms_x += "    x: {:.3f}      ".format(mol.contains[i].x)
        atoms_y += "    y: {:.3f}      ".format(mol.contains[i].y)
        atoms_z += "    z: {:.3f}      ".format(mol.contains[i].z)
    coordinates = "  {0:<18}{0:<18}{0:<18}{0:<18}".format('coordinates:')

    print("mol name: {}\n"
          "id: {}\n"
          "contains:".format(mol.name, mol.id))
    print(atoms_name+ "\n" +
          atoms_id + "\n" +
          coordinates + "\n" +
          atoms_x + "\n" +
          atoms_y + "\n" +
          atoms_z)


nrows = count_rows(filename)
atoms = load_atoms(filename, nrows)
mols = compute_molecules(filename, nrows, atoms)
