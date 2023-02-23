import polars as pl
import pyreadr
from rdkit import Chem
from rdkit.Chem import Descriptors

print(Chem.MolFromSmiles("CCCC"))

# Data from https://github.com/cran/rcdk/blob/master/data/bpdata.RData
# "A dataset containing the structures and associated boiling points for 277 molecules, primarily alkanes and substituted alkanes.""
# Actually has an additional initial column, molecule name
result = pyreadr.read_r('data/bpdata.RData')

# print(result)
df1 = result["bpdata"]
# print(df1)
# df1.to_csv('data/bp.csv')

df = pl.from_pandas(df1, include_index=True)
# print(df)
df = df.with_columns([
    # pl.col('BP').abs().alias('AbsBP'),
    pl.col('SMILES').apply(lambda s: Chem.MolFromSmiles(s)).alias('mol'),
])
print(df)

def wiener_index(m):
    res = 0
    amat = Chem.GetDistanceMatrix(m)
    num_atoms = m.GetNumAtoms()
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            res += amat[i][j]
    return res

df = df.with_columns([
    pl.col('mol').apply(lambda m: Chem.MolToSmiles(m)).alias('CanonicalSMILES'),
    pl.col('mol').apply(lambda m: Descriptors.MolWt(m)).alias('MolWt'),
    pl.col('mol').apply(lambda m: wiener_index(m)).alias('Wiener_Index'),
])
print(df)

print(df["SMILES"])