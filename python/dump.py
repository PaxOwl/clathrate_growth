"""
This file prints all the molecules and their respective atoms (dump purpose)
"""
from main import mols, print_mol

wat = 0
met = 0
for i in range(len(mols)):
    if "WAT" in mols[i].name:
        wat += 1
    elif "MET" in mols[i].name:
        met += 1

print("------------------ DUMPED DATA ------------------")
print("Dumped {} molecules, {} WATER and {} METHANE".format(len(mols),
                                                            wat, met))
print("-------------------------------------------------\n")

for i in range(len(mols)):
    print_mol(mols[i])
    print("")
