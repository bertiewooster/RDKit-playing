from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

salt = Chem.MolFromSmiles("[Na][Cl]")
mols_matrix = [[salt, salt], [salt, None, salt]]

# Minimal example: Only mols_matrix is supplied,
# result will be a drawing containing (where each row contains molecules)
# NaCl NaCl
# NaCl      NaCl
img = Chem.Draw.MolsMatrixToGridImage(mols_matrix)
# img is a PIL object for a PNG image file:
# <PIL.PngImagePlugin.PngImageFile image mode=RGB size=600x200 at 0x1648CC390>

# Exhaustive example: All parameters are supplied,
# result will be a drawing containing (where each row of molecules is followed by a row of legends):
# 0NaCl1                  0NaCl1
# no highlighting         sodium
# 0NaCl1                            0NaCl1
# chloride                          sodium and bond
legends_matrix = [["no highlighting", "sodium"], ["chloride", None, "sodium and bond"]]
highlightAtomLists_matrix = [[[],[]], [[0], None, [1]]]
highlightBondLists_matrix = [[[],[0]], [[], None, [0]]]

dopts = rdMolDraw2D.MolDrawOptions()
dopts.addAtomIndices = True

img_file = Chem.Draw.MolsMatrixToGridImage(mols_matrix=mols_matrix, subImgSize=(300, 400), legends_matrix=legends_matrix, highlightAtomLists_matrix=highlightAtomLists_matrix, highlightBondLists_matrix=highlightBondLists_matrix, useSVG=False, returnPNG=True, drawOptions=dopts)
img_file.save("salt.png")
# Drawing will be saved as PNG file salt.png