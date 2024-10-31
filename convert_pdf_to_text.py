import asyncio
import math
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
    return int(res)

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
    res = int(1./2*sum(sum(Distance==3)))
    
    return res

# Calculate delta t using Wiener equation 4
def calc_delta_t(n, delta_omega, delta_p):
    return (98/(n**2) * delta_omega) + (5.5 * delta_p)

def egloff(n):
  return 745.42 * math.log10(n + 4.4) - 689.4

def get_canonical_smiles(name):
    print(name)
    for compound in pcp.get_compounds(name, 'name'):
        print(f"  {name} {compound.canonical_smiles}")
        # if compound.canonical_smiles == "null":
        #     print(f"  {name} returned Null for canonical_smiles")
        time.sleep(0.01) # 0.01 second delay seems sufficient for PUG REST to not give Server Busy
        return compound.canonical_smiles

# Fix data from tables
# Convert these characters
#   « -> n
#   ^ -> 2
# Delete any line starting with a decimal point or comma

molecules = []
tobss = []

# with open("data/wiener_table_III_edited.txt") as f:
with open("data/wiener_table_II_edited.txt") as f:
    content = f.readlines()

# Show the file contents line by line.
# We added the comma to print single newlines and not double newlines.
# This is because the lines contain the newline character '\n'.
ignore_line_chars = (".", ",")

for line in content:
# Temporarily cutting down dataset size during coding
# for line in content[3:7]:
    if line[0] not in ignore_line_chars:
        if "2,2" in line:
            pass
        line_clean = line.replace("«", "n").replace("^", "2").replace("!", "l").replace("Ihexane", "lhexane").replace("Ioctane", "loctane").replace("Iheptane", "lheptane").replace("ro", "n").replace("pnpane", "propane").replace("pentaue","pentane").replace("raethyl", "methyl")
        end_marker = "ane "
        end_of_molecule = line_clean.find(end_marker) + len(end_marker)
        no_spaces_in_molecule = line_clean[:end_of_molecule].replace(" ", "")
        # print(f"{no_spaces_in_molecule=}")
        words = line_clean[end_of_molecule:].split()
        tobs = words[0]
        # print(f"{molecule_clean} {tobs}")
        molecules.append(no_spaces_in_molecule)
        tobss.append(float(tobs))

df = pl.DataFrame({"molecules": molecules,
                        #   "tobss": tobss, # For table III
                          "delta_t_obs": tobss, # For table II
                          })

# print(df)

# Debugging: Print molecule names
# print(df["molecules"])

# Debugging: Print molecule names and SMILES
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

print(df)

df = df.with_columns([
    pl.col('mol').apply(lambda m: Chem.MolToSmiles(m)).alias('CanonicalSMILES'),
    pl.col('mol').apply(lambda m: Descriptors.MolWt(m)).alias('MolWt'),
    pl.col('mol').apply(lambda m: wiener_index(m)).alias('omega'),
    pl.col('mol').apply(lambda m: CalculatePolarityNumber(m)).alias('p'),
    pl.col('mol').apply(lambda m: m.GetNumAtoms()).alias('n'),
])

print(df)

linear_alkanes = pl.DataFrame({"compound": ["n-Butane", "n-Pentane", "n-Hexane", "n-Heptane", "n-Octane", "n-Nonane", "n-Decane", "n-Undecane", "n-Dodecane"], 
                          "t0_obs": [-0.5, 36.1, 68.7, 98.4, 125.7, 150.8, 174.0, 195.8, 216.2],
                          "SMILES": ["CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", "CCCCCCCCCCC", "CCCCCCCCCCCC"],
                          "n": [4, 5, 6, 7, 8, 9, 10, 11, 12],
                          "omega0": [10, 20, 35, 56, 84, 120, 165, 220, 286],
                          "p0": [1, 2, 3, 4, 5, 6, 7, 8, 9],
                          })

linear_alkanes = linear_alkanes.with_columns([
    pl.col('n').apply(lambda n: egloff(n)).alias('t0_calc'),
    # pl.col('mol').apply(lambda m: wiener_index(m)).alias('omega0'),
    # pl.col('mol').apply(lambda m: CalculatePolarityNumber(m)).alias('p0'),
])

print(linear_alkanes)

# Join to copy in values from corresponding straight-chain alkane
df = df.join(linear_alkanes, on="n", how="inner", suffix="_lin_alkane")

# Calculate delta omega, p values
df = df.with_columns([
    pl.struct(["omega", "omega0"]).apply(lambda x: x["omega0"] - x["omega"]).alias("delta_omega"),
    pl.struct(["p", "p0"]).apply(lambda x: x["p0"] - x["p"]).alias("delta_p"),
])

# Calculate delta t
df = df.with_columns([
    pl.struct(["n", "delta_omega", "delta_p"]).apply(lambda x: calc_delta_t(x["n"], x["delta_omega"], x["delta_p"])).alias("delta_t_calc"),
])

df = df.with_columns([
    # Calculate deviation in delta t: obs - calc
    pl.struct(["delta_t_obs", "delta_t_calc"]).apply(lambda x: x["delta_t_obs"] - x["delta_t_calc"]).alias("Dev"),
    # Calculate t_calc
    pl.struct(["t0_calc", "delta_t_calc"]).apply(lambda x: x["t0_calc"] - x["delta_t_calc"]).alias("t_calc"),
    # Calculate t_obs
    pl.struct(["t0_obs", "delta_t_obs"]).apply(lambda x: x["t0_obs"] - x["delta_t_obs"]).alias("t_obs"),
])

print(df)

# Debugging: Print key t values for one molecule
row = df.filter(pl.col("molecules") == "2,3,4-Trimethylpentane")
row_t = row.select(["molecules", "SMILES", "t0_obs", "t0_calc", "delta_t_obs", "delta_t_calc", "t_obs", "t_calc"])
print(row_t)

