"""
This file prints all the molecules and their respective atoms (dump purpose)
"""
from main import mols
from classes import PrettyPrint


print("------------------ DUMPED DATA ------------------")
print("Dumped {} molecules, {} WATER and {} METHANE".format(len(mols['WAT'])
                                                            + len(mols['MET']),
                                                            len(mols['WAT']),
                                                            len(mols['MET'])))
print("-------------------------------------------------\n")

for key in mols:
    for i in range(len(mols[key])):
        PrettyPrint.print_mol(mols[key][i])
        print("")
