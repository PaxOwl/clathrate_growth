"""
This file contains the classes needed to store the data
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class Atom:
    name: str
    id: int
    x: float
    y: float
    z: float


@dataclass
class Molecule:
    name: str
    id: int
    contains: tuple
    rdf = np.ndarray


class PrettyPrint:
    @staticmethod
    def print_atom(atom: Atom):
        """
        Pretty printer for the Atom dataclass
        :param atom: the Atom dataclass to print
        :return: None
        """
        print("name: {}\n"
              "id: {:0>4}\n"
              "coordinates:\n"
              "    x: {}\n"
              "    y: {}\n"
              "    z: {}\n".format(atom.name, atom.id, atom.x, atom.y, atom.z))

    @staticmethod
    def print_mol(mol: Molecule):
        """
        Pretty printer for the Molecule dataclass
        :param mol: the Molecule to print
        :return: None
        """
        atoms_name = ''
        atoms_id = ''
        atoms_x = ''
        atoms_y = ''
        atoms_z = ''

        for i in range(len(mol.contains)):
            atoms_name += "{:>2}atom name: {:<5}".format('',
                                                         mol.contains[i].name)
            atoms_id += "{0:>2}atom id: {1:0>4}{0:>3}".format('',
                                                              mol.contains[
                                                                  i].id)
            atoms_x += "{0:<4}x: {1:>6.3f}{0:>6}".format('', mol.contains[i].x)
            atoms_y += "{0:<4}y: {1:>6.3f}{0:>6}".format('', mol.contains[i].y)
            atoms_z += "{0:<4}z: {1:>6.3f}{0:>6}".format('', mol.contains[i].z)

        coordinates = "{0:>2}{1:<18}{1:<18}{1:<18}{1:<18}"\
                      .format('', 'coordinates:')

        print("mol name: {}\n"
              "id: {:0>4}\n"
              "contains:".format(mol.name, mol.id))
        print(atoms_name.rstrip() + '\n' +
              atoms_id.rstrip() + '\n' +
              coordinates.rstrip() + '\n' +
              atoms_x.rstrip() + '\n' +
              atoms_y.rstrip() + '\n' +
              atoms_z.rstrip())
