from rdkit import Chem
from rdkit.Chem import Draw

mol = Chem.MolFromSmiles("OCCC(F)CC(CCF)CC(O)(CCCCCl)O")
mols = [mol, mol]

img_SVG_false_PNG_false = Draw.MolsToGridImage(mols, useSVG=False, returnPNG=False)
img_SVG_false_PNG_false.save("img_SVG_false_PNG_false.PNG")
print(f"{img_SVG_false_PNG_false=}")
# print()

# img_SVG_false_PNG_true = Draw.MolsToGridImage(mols, useSVG=False, returnPNG=True)
# print(f"{img_SVG_false_PNG_true=}")
# img_SVG_false_PNG_true.save("img_SVG_false_PNG_true.PNG")

# img_SVG_true_PNG_false = Draw.MolsToGridImage(mols, useSVG=True, returnPNG=False)
# print(f"{img_SVG_true_PNG_false=}")
# img_SVG_true_PNG_false.save("img_SVG_true_PNG_false.SVG")

# img_SVG_true_PNG_true = Draw.MolsToGridImage(mols, useSVG=True, returnPNG=True)
# print(f"{img_SVG_true_PNG_true=}")
# img_SVG_true_PNG_true.save("img_SVG_true_PNG_true.SVG")

a = Chem.MolToSmiles(mol, doRandom = True)
print(a)
b = Chem.MolToSmiles(mol, doRandom = True)
print(b)
