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
    def __init__(self, val):
        self.val = val
        self.nodes = []

    def add_node(self, val):
        self.nodes.append(NonBinTree(val))
        return self.nodes[-1]

    def __repr__(self):
        return f"NonBinTree({self.val}): {self.nodes}"

    def get_ncols(self):
        ncols = 0
        if len(self.nodes) > 0:
            # for node_counter in range(len(self.nodes)):
            #     ncols += self.get_ncols(self.nodes[node_counter])
            for node in self.nodes:
                ncols += node.get_ncols()
        else:
            ncols += 1
        return ncols


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
print(root)
print(f"{f1.get_ncols()=}")
print(f"{f2.get_ncols()=}")
print(f"{root.get_ncols()=}")
