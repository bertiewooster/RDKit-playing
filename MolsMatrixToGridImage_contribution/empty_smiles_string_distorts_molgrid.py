from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

salt = Chem.MolFromSmiles("[Na][Cl]")

mols_None = [salt, None, salt]
img_None = Chem.Draw.MolsToGridImage(mols_None)
img_None.save("img_None.png")

smiles = ["[Na][Cl]", "", "[Na][Cl]"]
mols_smiles = [Chem.MolFromSmiles(string) for string in smiles]
# print(f"{mols_smiles=}")
img_smiles = Chem.Draw.MolsToGridImage(mols_smiles)
img_smiles.save("img_smiles.png")

empty_smiles = Chem.MolFromSmiles("")
print(f"{empty_smiles=}")

bad = Chem.MolFromSmiles("alskjfowijfjsdfo")
print(f"{bad=}")
