from rdkit import Chem
from rdkit.Chem import rdqueries

m = Chem.MolFromSmiles('F[P-](F)(F)(F)(F)F.CN(C)C(F)=[N+](C)C',sanitize=False)
m
# Build a query for the P
q = rdqueries.AtomNumEqualsQueryAtom(15)

# Select the first and only P
phosphorus = m.GetAtomsMatchingQuery(q)[0]

print(phosphorus.GetHybridization())

# Regenerate computed properties like implicit valence and ring information
m.UpdatePropertyCache(strict=False)

# Apply several sanitization rules
Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
m

print(phosphorus.GetHybridization())
