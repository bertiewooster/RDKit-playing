from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw, rdFMCS

mol = Chem.MolFromSmiles("C")
print(mol)

row1 = [mol] * 10000
Draw.MolsToGridImage(row1)
