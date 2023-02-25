import asyncio
import re
import ssl
import time

import aiohttp
import polars as pl
import pubchempy as pcp
from pypdf import PdfReader
from rdkit import Chem
from rdkit.Chem import Descriptors

ssl._create_default_https_context = ssl._create_unverified_context

reader = PdfReader("Journal articles/Wiener index ja01193a005.pdf")
number_of_pages = len(reader.pages)
page = reader.pages[2]
text = page.extract_text()

with open('data/wiener_page_2.txt', 'w') as f:
    f.write(text)

def wiener_index(m):
    res = 0
    amat = Chem.GetDistanceMatrix(m)
    num_atoms = m.GetNumAtoms()
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            res += amat[i][j]
    return res

def CalculatePolarityNumber(mol):
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
    res = 1./2*sum(sum(Distance==3))
    
    return res

# Calculate delta t using Wiener equation 4
def calc_delta_t(n, delta_omega, delta_p):
    return (98/(n**2) * delta_omega) + (5.5 * delta_p)

# delta_t = calc_delta_t(4, 1, 1) # should give 11.625; does
# delta_t = calc_delta_t(8, 26, -4) # should give 17.8; does
delta_t = calc_delta_t(8, 26, -4) # should give 17.8; does

def get_canonical_smiles(name):
    for compound in pcp.get_compounds(name, 'name'):
        print(f"{name} {compound.canonical_smiles}")
        # if compound.canonical_smiles == "null":
        #     print(f"  {name} returned Null for canonical_smiles")
        # time.sleep(0.5)
        return compound.canonical_smiles

# Fix data from tables
# Convert these characters
#   « -> n
#   ^ -> 2
# Delete any line starting with a decimal point or comma

molecules = []
tobss = []

with open("data/wiener_table_III_edited.txt") as f:
    content = f.readlines()

# Show the file contents line by line.
# We added the comma to print single newlines and not double newlines.
# This is because the lines contain the newline character '\n'.
ignore_line_chars = (".", ",")
# for line in content:

# for line in content:
# Temporarily cutting down dataset size during coding
for line in content[0:3]:
    if line[0] not in ignore_line_chars:
        end_marker = "ane "
        end_of_molecule = line.find(end_marker) + len(end_marker)
        no_spaces_in_molecule = line[:end_of_molecule].replace(" ", "")
        # print(f"{no_spaces_in_molecule=}")
        molecule_clean = no_spaces_in_molecule.replace("«", "n").replace("^", "2").replace("!", "l").replace("Ihexane", "lhexane").replace("Ioctane", "loctane").replace("Iheptane", "lheptane")

        words = line[end_of_molecule:].split()
        tobs = words[0]
        # print(f"{molecule_clean} {tobs}")
        molecules.append(molecule_clean)
        tobss.append(tobs)

df = pl.DataFrame({"molecules": molecules,
                          "tobss": tobss})

print(df)

# print(df["molecules"])

# for molecule in df["molecules"]:
#     s = get_canonical_smiles(molecule)
#     if s == "null":
#         print(f"{molecule} {s}")

df = df.with_columns([
    pl.col('molecules').apply(lambda s: get_canonical_smiles(s)).alias('SMILES'),
])

# print(df)

df = df.with_columns([
    pl.col('SMILES').apply(lambda s: Chem.MolFromSmiles(s)).alias('mol'),
])

df = df.with_columns([
    pl.col('mol').apply(lambda m: Chem.MolToSmiles(m)).alias('CanonicalSMILES'),
    pl.col('mol').apply(lambda m: Descriptors.MolWt(m)).alias('MolWt'),
    pl.col('mol').apply(lambda m: wiener_index(m)).alias('omega'),
    pl.col('mol').apply(lambda m: CalculatePolarityNumber(m)).alias('p'),
    pl.col('mol').apply(lambda m: m.GetNumAtoms()).alias('n'),
])

print(df)

linear_alkanes = pl.DataFrame({"compound": ["n-Butane", "n-Pentane", "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Decane", "n-Undecane", "n-Dodecane"], 
                          "t0": [-0.5, 36.1, 68.7, 98.4, 125.7, 150.8, 174.0, 195.8, 216.2],
                          "SMILES": ["CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", "CCCCCCCCCCC", "CCCCCCCCCCCC"],
                          "n": [4, 5, 6, 7, 8, 9, 10, 11, 12],
                          "omega0": [10, 20, 35, 56, 84, 120, 165, 220, 286],
                          "p0": [1, 2, 3, 4, 5, 6, 7, 8, 9],
                          })

print(linear_alkanes)

# Join to copy in values from corresponding straight-chain alkane
df = df.join(linear_alkanes, on="n", how="inner", suffix="_lin_alkane")

# Calculate delta omega, p values
df = df.with_columns([
    pl.struct(["omega", "omega0"]).apply(lambda x: x["omega0"] - x["omega"]).alias("delta_omega"),
    pl.struct(["p", "p0"]).apply(lambda x: x["p0"] - x["p"]).alias("delta_p"),
])

# Calculate delta t
# df = df.with_columns([
#     pl.struct(["n", "delta_omega", "delta_p"]).apply(lambda x: calc_delta_t(x["n"], x["delta_omega"], x["delta_p"])).alias("delta_t"),
# ])

# df = df.with_columns([
#     pl.col('n').apply(lambda s: Chem.MolFromSmiles(s)).alias('mol'),
# ])

# Calculate delta t values
# df = df.with_columns([
#     pl.struct(["n", "delta_omega", "delta_p"]).apply(lambda x: calc_delta_t(x["n"], x["delta_omega"], x["delta_p"])).alias("delta_t"),
# ])

print(df)


# smiless = []
# mols = []
# wieners = []
# polarities = []

# for mol in mols:
#     wieners = 