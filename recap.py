from rdkit import Chem
from rdkit.Chem import Draw, Recap

# cisapride = Chem.MolFromSmiles("Clc1cc(c(OC)cc1N)C(=O)NC3CCN(CCCOc2ccc(F)cc2)CC3OC")
# # Draw.MolToImage(cisapride)
# hierarch = Recap.RecapDecompose(cisapride)
# ks = hierarch.children.keys()
# ks = sorted(ks)
# # print(f"{ks=}")

# # all_children = hierarch.GetAllChildren()
# # print(f"{all_children=}")

# generations = [hierarch.smiles]

# level = hierarch
# level_counter = 0


# def get_children(base_node, level, cols=0):
#     for smiles, node in base_node.children.items():
#         spaces = " " * level
#         # print(f"{spaces}{smiles=}")
#         children = node.children.keys()
#         children = sorted(children)
#         if len(children) > 0:
#             get_children(node, level + 1, cols)
#             # print(f"{children=}")
#         else:
#             cols += 1
#             # print(f"{cols=}")
#             # return cols


# cols = get_children(hierarch, 0, 0)
# # print(f"{cols=}")


class NonBinTree:
    """Adapted from https://stackoverflow.com/questions/60579330/non-binary-tree-data-structure-in-python#60579464"""

    def __init__(self, val):
        self.val = val
        self.nodes = []

    def add_node(self, val):
        self.nodes.append(NonBinTree(val))
        return self.nodes[-1]

    def __repr__(self):
        return f"NonBinTree({self.val}): {self.nodes}"

    def get_ncols(self):
        self.ncols = 0
        if len(self.nodes) > 0:
            # If there are nodes under this one, call get_ncols on them recursively
            for node in self.nodes:
                self.ncols += node.get_ncols()
        else:
            # If there are no nodes under this one, add 1 for this node
            self.ncols += 1
        return self.ncols

    def get_max_depth(self):
        max_depth = 0
        if len(self.nodes) > 0:
            for node in self.nodes:
                this_depth = node.get_max_depth()
                max_depth = max(this_depth + 1, max_depth)
        else:
            max_depth = max(1, max_depth)
        self.max_depth = max_depth
        return self.max_depth

    def get_grid(self):
        self.get_ncols()
        # Create top row: Molecule, then the rest of columns are blank (empty strings)
        grid = [self.val] + [""] * (self.ncols - 1)
        return grid


root = NonBinTree("cisapride")
f1 = root.add_node("f1")
f2 = root.add_node("f2")
f21 = f2.add_node("f21")
f22 = f2.add_node("f22")
f23 = f2.add_node("f23")
f26 = f2.add_node("f24")
f3 = root.add_node("f3")
f4 = root.add_node("f4")
f41 = f4.add_node("f41")
f42 = f4.add_node("f42")
f5 = root.add_node("f5")
f6 = root.add_node("f6")
f61 = f6.add_node("f61")
f62 = f6.add_node("f62")

# Adding for debugging--not real for cisapride
f7 = root.add_node("f7")
f71 = f7.add_node("f71")
f72 = f7.add_node("f71")
f711 = f71.add_node("f711")
f712 = f71.add_node("f712")
f7111 = f711.add_node("f7111")
f7112 = f712.add_node("f7112")

print(root)
# print(f"{f1.get_ncols()=}")
# print(f"{f2.get_ncols()=}")
# print(f"{root.get_ncols()=}")
# grid = root.get_grid()
# print(f"{grid=}, {len(grid)=}")
# print(f"{root.get_max_depth()=}, should be 3")
print(f"{root.get_max_depth()=}, should be 5")
print(f"{f1.get_max_depth()=}, should be 1")
print(f"{f4.get_max_depth()=}, should be 2")
print(f"{f41.get_max_depth()=}, should be 1")
print(f"{f42.get_max_depth()=}, should be 1")
print(f"{f7111.get_max_depth()=}, should be 1")
print(f"{f711.get_max_depth()=}, should be 2")
print(f"{f71.get_max_depth()=}, should be 3")
print(f"{f72.get_max_depth()=}, should be 1")
print(f"{f7.get_max_depth()=}, should be 4")

