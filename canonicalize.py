from rdkit import Chem

smile:str = "C1=CC=CC=C1"
m:Chem.Mol = Chem.MolFromSmiles(smile)
print(f"{m=}")
m2:Chem.rdchem.Mol = Chem.MolFromSmiles(smile)
print(f"{m2=}")
canon_smile:str = Chem.MolToSmiles(m)
print(f"{canon_smile=}")
canon_smile2:str = Chem.MolToSmiles(m2)
print(f"{canon_smile2=}")

# smls = ["C1=CC=CC=C1", "C1C=CC=CC=1"]
# mols = [Chem.MolFromSmiles(sml) for sml in smls]
# canon_smls = [Chem.MolToSmiles(mol) for mol in mols]
# print(canon_smls)
