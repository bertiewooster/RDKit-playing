import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

# Define the Cython extension module
ext_modules = [
    Extension(
        "molfromsmiles",
        sources=["molfromsmiles.pyx"],
        include_dirs=[np.get_include(), "venv/lib/python3.11/site-packages/rdkit/include"],
        libraries=["RDKitChemTransforms", "RDKitChemSmiles", "RDKitChemGraph", "RDKitSubstructMatch", "RDKitFingerprints", "RDKitDescriptors", "RDKitDataStructs", "RDKitRDGeometryLib", "RDKitRDGeneral", "RDKitRDBoost"],
        library_dirs=["venv/lib/python3.11/site-packages/rdkit/lib"],
        language="c++"
    )
]

setup(
    ext_modules=cythonize(ext_modules)
)
