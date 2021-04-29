"""
Core part of the program
"""


from analysis import atoms, mols, print_mol


wat = 0
met = 0
for i in range(len(mols)):
    if mols[i].name == "WAT":
        wat += 1
    elif mols[i].name == "MET":
        met += 1

print("Water: {}".format(wat))
print("Methane: {}\n".format(met))

print_mol(mols[0])
