        from rdkit import Chem
        from rdkit.Chem.Draw import MolsMatrixToGridImage, rdMolDraw2D
        FCl = Chem.MolFromSmiles("FCl")
        mols_matrix = [[FCl, FCl], [FCl, None, FCl]]

        # Minimal example: Only mols_matrix is supplied,
        # result will be a drawing containing (where each row contains molecules):
        # Cl-F    Cl-F
        # Cl-F             Cl-F
        img = MolsMatrixToGridImage(mols_matrix)
        img.save("MolsMatrixToGridImage_minimal.png")
        # img is a PIL object for a PNG image file like:
        # <PIL.PngImagePlugin.PngImageFile image mode=RGB size=600x200 at 0x1648CC390>
        # Drawing will be saved as PNG file MolsMatrixToGridImage_minimal.png

        # Exhaustive example: All parameters are supplied,
        # result will be a drawing containing (where each row of molecules is followed by a row of legends):
        # 1 Cl-F 0             1 Cl-F 0
        # no highlighting       bond highlighted         
        # 1 Cl-F 0                                1 Cl-F 0
        # sodium highlighted                       chloride and bond highlighted
        legends_matrix = [["no highlighting", "bond highlighted"], 
        ["F highlighted", "", "Cl and bond highlighted"]]
        highlightAtomLists_matrix = [[[],[]], [[0], None, [1]]]
        highlightBondLists_matrix = [[[],[0]], [[], None, [0]]]

        dopts = rdMolDraw2D.MolDrawOptions()
        dopts.addAtomIndices = True

        img_file = MolsMatrixToGridImage(mols_matrix=mols_matrix, subImgSize=(300, 400), 
        legends_matrix=legends_matrix, highlightAtomLists_matrix=highlightAtomLists_matrix, 
        highlightBondLists_matrix=highlightBondLists_matrix, useSVG=False, returnPNG=True, drawOptions=dopts)
        img_file.save("MolsMatrixToGridImage_exhaustive.png")
        # Drawing will be saved as PNG file MolsMatrixToGridImage_exhaustive.png
