"""
This file prints all the molecules and their respective atoms (dump purpose)
"""
from main import mols, print_mol

wat = 0
met = 0
for key in mols:
    for i in range(len(mols[key])):
        if "WAT" in mols[key][i].name:
            wat += 1
        elif "MET" in mols[key][i].name:
            met += 1

print("------------------ DUMPED DATA ------------------")
print("Dumped {} molecules, {} WATER and {} METHANE".format(len(mols),
                                                            wat, met))
print("-------------------------------------------------\n")

for key in mols:
    for i in range(len(mols[key])):
        print_mol(mols[key][i])
        print("")
