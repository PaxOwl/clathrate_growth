"""
This file prints all the molecules and their respective atoms (dump purpose)
"""
from main import mols, print_mol


print("------------------ DUMPED DATA ------------------")
print("Dumped {} molecules, {} WATER and {} METHANE".format(len(mols),
                                                            len(mols['WAT']),
                                                            len(mols['MET'])))
print("-------------------------------------------------\n")

for key in mols:
    for i in range(len(mols[key])):
        print_mol(mols[key][i])
        print("")
