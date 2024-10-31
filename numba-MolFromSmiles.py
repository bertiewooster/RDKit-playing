import datetime

from numba import jit
from rdkit import Chem


@jit
def mols_from_smiles(smiles):
    mols = []
    for sml in smiles:
        mols.append(Chem.MolFromSmiles(sml))
    return mols

# count = 1_000
count = 1_000_000
smls = ["CC"] * count

start = datetime.datetime.now()
mols_from_smiles(smls)
print(f"Converting {count} sml to mol took {datetime.datetime.now() - start}s")

