import pandas as pd
from rdkit.Chem import PandasTools

PandasTools.InstallPandasTools() # <- only necessary during testing, you don't need to do this
import os

from rdkit import RDConfig

antibiotics = pd.DataFrame(columns=['Name','Smiles'])
antibiotics = pd.concat([antibiotics, pd.DataFrame.from_records([
    {'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
  'Name':'Penicilline G',
  'Integer': 1,
  'Float': 1.01,
  }])], ignore_index=True) #Penicilline G
# antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
#   'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O',
#   'Name':'Tetracycline'}])], ignore_index=True) #Tetracycline
# antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
#   'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C',
#   'Name':'Ampicilline'}])], ignore_index=True) #Ampicilline
PandasTools.AddMoleculeColumnToFrame(
    antibiotics,
    smilesCol="Smiles"
)
print([str(x) for x in  antibiotics.columns])
print(antibiotics)
print(antibiotics.to_html())