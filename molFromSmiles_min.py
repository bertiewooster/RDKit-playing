from rdkit import Chem

sml = "CO"
smls = sml * 1_000
mols = [Chem.MolFromSmiles(sml) for sml in smls]
print(mols[0])
