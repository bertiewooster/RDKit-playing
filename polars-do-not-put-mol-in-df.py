import sys
from functools import wraps
from time import time

import polars as pl
from rdkit import Chem


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print('func:%r args:[%r, %r] took: %2.4f sec' % \
          (f.__name__, args, kw, te-ts))
        return result
    return wrap

count = 1_000
# count = 1_000_000
# count = 5
smls = ["CC"] * count

def wiener_index(m: Chem.Mol):
    """
    From https://sourceforge.net/p/rdkit/mailman/message/36802142/ by Greg Landrum
    :returns: Wiener index, aka path number
    :rtype: int
    :param m: RDKit molecule
    """
    res = 0
    amat = Chem.GetDistanceMatrix(m)
    num_atoms = m.GetNumAtoms()
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            res += amat[i][j]
    return int(res)

def CalculatePolarityNumber(mol: Chem.Mol) -> int:
    """
    #################################################################
    Copyright BSD 3-Clause "New" or "Revised" License
    Author : gadsbyfly
    https://codesuche.com/view-source/python/gadsbyfly/PyBioMed/

    Calculation of Polarity number.
    It is the number of pairs of vertexes at
    distance matrix equal to 3
    ---->Pol
    Usage:
        result=CalculatePolarityNumber(mol)
        Input: mol is a molecule object
        Output: result is a numeric value
    #################################################################
    """
    Distance = Chem.GetDistanceMatrix(mol)
    res = int(1.0 / 2 * sum(sum(Distance == 3)))

    return res

def mol_props(sml):
    """
    Convert SMILES to an RDKit molecule, then calculate various properties of it
    :returns: dictionary of molecular properties
    :param sml: SMILES to convert to a molecule
    """
    mol = Chem.MolFromSmiles(sml)
    CanonicalSMILES = Chem.MolToSmiles(mol)
    omega = wiener_index(mol)
    p = CalculatePolarityNumber(mol)
    n = mol.GetNumAtoms()
    return dict(
        CanonicalSMILES=CanonicalSMILES,
        omega=omega,
        p=p,
        n=n,
        )

@timing
def add_mol_props(df):
    return df.with_columns(
    molecular_props = pl.col('SMILES').map_elements(mol_props)
    ).unnest('molecular_props')

df = pl.DataFrame({'SMILES': smls})
df = add_mol_props(df)
print(df)
print(f"{df.estimated_size()=}")

@timing
def add_mols_and_props(df):
    return df.with_columns(
        [
            pl.col("SMILES").map_elements(lambda s: Chem.MolFromSmiles(s)).alias("mol"),
        ]
    ).with_columns(
    [
        pl.col("mol").map_elements(lambda m: Chem.MolToSmiles(m)).alias("CanonicalSMILES"),
        pl.col("mol").map_elements(lambda m: wiener_index(m)).alias("omega"),
        pl.col("mol").map_elements(lambda m: CalculatePolarityNumber(m)).alias("p"),
        pl.col("mol").map_elements(lambda m: m.GetNumAtoms()).alias("n"),
    ]
)

df2 = pl.DataFrame({'SMILES': smls})
df2 = add_mols_and_props(df2)
print(df2)
print(f"{df2.estimated_size()=}")
