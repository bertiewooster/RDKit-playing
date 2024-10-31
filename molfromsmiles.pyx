# distutils: language = c++
# distutils: include_dirs = venv/lib/python3.11/site-packages/rdkit/include
# distutils: libraries = RDKitChemTransforms RDKitChemSmiles RDKitChemGraph RDKitSubstructMatch RDKitFingerprints RDKitDescriptors RDKitDataStructs RDKitRDGeometryLib RDKitRDGeneral RDKitRDBoost

from rdkit cimport Chem
from rdkit.ROMol cimport Mol

cpdef Chem.Mol[::1] mols_from_smiles(list smiles):
    cdef list mols = []
    for sml in smiles:
        mols.append(Chem.MolFromSmiles(sml))
    return mols
