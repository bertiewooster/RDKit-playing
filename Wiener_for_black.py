import math

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn
from py2opsin import py2opsin
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import IPythonConsole


# Cheminformatic utilities
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
    for j in range(i+1, num_atoms):
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
    res = int(1./2*sum(sum(Distance==3)))
    
    return res

def calc_Δt(n: int, Δomega: int, Δp: int) -> float:  
  """
  Calculate Δt using Wiener equation 4
  https://pubs.acs.org/doi/10.1021/ja01193a005
  :returns: Δt, difference in boiling point between alkane and its structral 
  isomer of a linear alkane
  :param n: number of carbon atoms
  :param Δomega: Wiener index of straight-chain alkane minus this molecule
  :param Δp: polarity number of straight-chain alkane minus this molecule
  """
  return (98/(n**2) * Δomega) + (5.5 * Δp)

def egloff(n: int) -> float:
  """
  Calculate boiling point for linear alkane using Egloff's equation
  https://pubs.acs.org/doi/pdf/10.1021/j150402a006
  :returns: Δt, difference in boiling point between alkane and its structral 
  isomer of a linear alkane
  :param n: number of carbon atoms
  """
  return 745.42 * math.log10(n + 4.4) - 689.4



  molecules = []
tables = []
ts_read_in = []
ts_which = []

# Dictionary of typos and corrections. Italicized "n" in "n-" is particularly difficult for OCR.
replace_typos = {
    "w-": "n-", 
    "ro-": "n-",
    "«-": "n-",
    "^": "2",
    "!": "l",
    "thyI": "thyl",
    "Methyt": "Methyl",
    "raethyl": "methyl",
    "pentaue": "pentane",
}

for table_num in ('II', 'III'):
  table_file = f"/content/drive/MyDrive/data/wiener_table_{table_num}_edited.txt"
  with open(table_file) as f:
    content = f.readlines()

    # Ignore lines that start with period or comma--these were incorrectly split across two lines
    ignore_line_chars = (".", ",")

    for line in content:
      if line[0] not in ignore_line_chars:
        line_clean = line
        for typo, correct in replace_typos.items():
          line_clean = line_clean.replace(typo, correct)
        end_marker = "ane "
        end_of_molecule = line_clean.find(end_marker) + len(end_marker)
        no_spaces_in_molecule = line_clean[:end_of_molecule].replace(" ", "")
        words = line_clean[end_of_molecule:].split()
        t_read_in = words[0]

        # Some table entries have no observed data; we process only molecules with observed data
        if t_read_in != "Null":
          molecules.append(no_spaces_in_molecule)
          tables.append(table_num)
          ts_read_in.append(float(t_read_in))
          if table_num == 'II':
            ts_which.append("Δt")
          elif table_num == 'III':
            ts_which.append("t")



smiles = py2opsin(molecules)



df = pl.DataFrame({"molecule": molecules,
                          "table": tables,
                          "Smiles": smiles,
                          "t_read_in": ts_read_in,
                          "t_which": ts_which,
                          })
df = df.with_row_count(name="Compound_Id", offset=1)

print(df)



df = df.with_columns([
  pl.col('Smiles').apply(lambda s: Chem.MolFromSmiles(s)).alias('mol'),
])