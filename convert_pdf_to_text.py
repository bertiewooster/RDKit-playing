import asyncio
import re
import ssl
import time

import aiohttp
import polars as pl
import pubchempy as pcp
from pypdf import PdfReader
from rdkit import Chem

ssl._create_default_https_context = ssl._create_unverified_context

# reader = PdfReader("Journal articles/Wiener index ja01193a005.pdf")
# number_of_pages = len(reader.pages)
# page = reader.pages[1]
# text = page.extract_text()

def get_canonical_smiles(name):
    for compound in pcp.get_compounds(name, 'name'):
        print(f"{molecule} {compound.canonical_smiles}")
        if compound.canonical_smiles == "null":
            print(f"  {name} returned Null for canonical_smiles")
        time.sleep(0.5)
        return compound.canonical_smiles

# Fix data from tables
# Convert these characters
#   « -> n
#   ^ -> 2
# Delete any line starting with a decimal point or comma

molecules = []
tobss = []

with open("data/wiener_table_III.txt") as f:
    content = f.readlines()

# Show the file contents line by line.
# We added the comma to print single newlines and not double newlines.
# This is because the lines contain the newline character '\n'.
ignore_line_chars = (".", ",")
for line in content:
    if line[0] not in ignore_line_chars:
        end_marker = "ane "
        end_of_molecule = line.find(end_marker) + len(end_marker)
        no_spaces_in_molecule = line[:end_of_molecule].replace(" ", "")
        # print(f"{no_spaces_in_molecule=}")
        molecule_clean = no_spaces_in_molecule.replace("«", "n").replace("^", "2").replace("!", "l").replace("Ihexane", "lhexane").replace("Ioctane", "loctane").replace("Iheptane", "lheptane")
        print(f"{molecule_clean=}")

        words = line[end_of_molecule:].split()
        tobs = words[0]
        molecules.append(molecule_clean)
        tobss.append(tobs)

df = pl.DataFrame({"molecules": molecules,
                          "tobss": tobss})

print(df)

print(df["molecules"])

for molecule in df["molecules"]:
    s = get_canonical_smiles(molecule)
    if s == "null":
        print(f"{molecule} {s}")

df = df.with_columns([
    pl.col('molecules').apply(lambda s: get_canonical_smiles(s)).alias('SMILES'),
])

# print(df)

# df = df.with_columns([
#     pl.col('SMILES').apply(lambda s: Chem.MolFromSmiles(s)).alias('mol'),
# ])

# print(df)


# smiless = []
# mols = []
# wieners = []
# polarities = []

# for mol in mols:
#     wieners = 