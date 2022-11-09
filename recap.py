from itertools import chain

from rdkit import Chem
from rdkit.Chem import Draw, Recap

# def flatten_to_list(nested_list):
#     if nested_list == []:
#         return []
#     mylist = nested_list.copy()
#     new_list = []
#     for item in nested_list:
#         while isinstance(item, (list, tuple)):
#             for subitem in item:
#                 while isinstance(subitem, (list, tuple)):
#                     subitem = list(chain(*subitem))
#                 new_list += subitem
#             # mylist[0] = list(chain(*mylist))
#         # new_list += item
#     return mylist


a = [1, 2]
b = [[1], [2]]
c = [[[1], [2], [3]]]

# print(flatten_to_list(a))
# print(flatten_to_list(b))
# print(flatten_to_list(c))

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

def concat(grid1, grid2):
    combined = []
    for row_counter in range(len(grid1)):
        combined += [grid1[row_counter] + grid2[row_counter]]
    return combined

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
        self_cols = self.get_ncols()
        self_rows = self.get_max_depth()

        # Create top row: Node value, then the rest of columns are blank (empty strings)
        grid = [[self.val] + [""] * (self.ncols - 1)]
        # for rows_to_pad in range(self_rows - 1):
        #     grid += [[]]

        n_nodes = len(self.nodes)
        if n_nodes > 0:
            # nodes_grid = [[""]] * n_nodes
            nodes_grid = [[]]
            for node_counter, node in enumerate(self.nodes):
                # nodes_grid.append(node.get_grid())
                node_grid = node.get_grid()

                if self.val == "CCNN":
                    ...
                # Add spacer rows if needed
                node_grid_rows = len(node_grid)
                rows_padding = self_rows - node_grid_rows - 1
                for padding in range(rows_padding):
                    node_grid.append(["" * rows_padding])
                # nodes_grid[node_counter] = node_grid

                nodes_grid = concat(nodes_grid, node_grid)

            # combined = []
            # for row_counter in range(node_grid_rows):
            #     for node_number in range(len(self.nodes)):
            #         combined += [nodes_grid[node_number][row_counter]]
            grid += nodes_grid
            # grid[-1] = nodes_grid

        return grid


nested = [
    ((["CCN"], ["CCNN", ""]), ["CCNNN"]),
    (([""], [(["CCNNO"], ["CCNNOO"])]), [""]),
]

root = NonBinTree("Fe")
f1 = root.add_node("C")
f2 = root.add_node("CC")
f21 = f2.add_node("CCN")
f22 = f2.add_node("CCNN")
f221 = f22.add_node("CCNNO")
f222 = f22.add_node("CCNNOO")
f23 = f2.add_node("CCNNN")
f3 = root.add_node("CCC")
f4 = root.add_node("CCCC")
f41 = f4.add_node("CCCCN")
f42 = f4.add_node("CCCCNN")

# print(f"{f1.get_grid()=}")
# print(f"{f21.get_grid()=}")
print(f"{f22.get_grid()=}")
# print(f"{f221.get_grid()=}")
# print(f"{f23.get_grid()=}")
# print(f"{f2.get_grid()=}")
# print(f"{f3.get_grid()=}")
# print(f"{f4.get_grid()=}")

# root = NonBinTree("cisapride")
# f1 = root.add_node("f1")
# f2 = root.add_node("f2")
# f21 = f2.add_node("f21")
# f22 = f2.add_node("f22")
# f23 = f2.add_node("f23")
# f24 = f2.add_node("f24")
# f3 = root.add_node("f3")
# f4 = root.add_node("f4")
# f41 = f4.add_node("f41")
# f42 = f4.add_node("f42")
# f5 = root.add_node("f5")
# f6 = root.add_node("f6")
# f61 = f6.add_node("f61")
# f62 = f6.add_node("f62")

# Adding for debugging--not real for cisapride
# f7 = root.add_node("f7")
# f71 = f7.add_node("f71")
# f72 = f7.add_node("f71")
# f711 = f71.add_node("f711")
# f712 = f71.add_node("f712")
# f7111 = f711.add_node("f7111")
# f7112 = f712.add_node("f7112")

# def combine_grids(*grids):
#     grid = list(zip(*grids))
#     grid_flat = []
#     for row in grid:
#         # Remove tuple, flatten into single list
#         flat_row = [item for item in row]
#         grid_flat.append(flat_row)
#         # print(f"{flat_row=}")


# print(root)
# print(f"{f1.get_ncols()=}")
# print(f"{f2.get_ncols()=}")
# print(f"{root.get_ncols()=}")
# grid = root.get_grid()
# print(f"{grid=}, {len(grid)=}")
# print(f"{root.get_max_depth()=}, should be 3")
# print(f"{root.get_max_depth()=}, should be 5")
# print(f"{f1.get_max_depth()=}, should be 1")
# print(f"{f4.get_max_depth()=}, should be 2")
# print(f"{f41.get_max_depth()=}, should be 1")
# print(f"{f42.get_max_depth()=}, should be 1")
# print(f"{f7111.get_max_depth()=}, should be 1")
# print(f"{f711.get_max_depth()=}, should be 2")
# print(f"{f71.get_max_depth()=}, should be 3")
# print(f"{f72.get_max_depth()=}, should be 1")
# print(f"{f7.get_max_depth()=}, should be 4")

# molecule_grid = [["molecule"] + [""] * 8]
# molecule_grid = [["O"] + [""] * 8]
# print(f"{molecule_grid=}")

# f1_grid = [["f1"], [""], [""]]
# f2_grid = [["f2", "", "", ""], ["f21", "f22", "", "f23"], ["", "f221", "f222", ""]]
# f3_grid = [["f3"], [""], [""]]
# f4_grid = [["f4", ""], ["f41", "f42"], ["", "", ""]]

# f1_grid = [["C"], [""], [""]]
# f21_grid = ["CCN"]

# f221_grid = ["CCNNO"]
# f222_grid = ["CCNNOO"]
# f22_grid = [["CCNN", ""], f221_grid + f222_grid]

# f22_sub_grids = combine_grids(f221_grid, f222_grid)
# print(f"{f22_sub_grids=}")

# f23_grid = ["CCNNN"]

# f2_grid =

# f2_grid = [
#     ["CC", "", "", ""],
#     ["CCN", "CCNN", "", "CCNNN"],
#     ["", "CCNNO", "CCNNOO", ""],
# ]
# f3_grid = [["CCC"], [""], [""]]
# f4_grid = [["CCCC", ""], ["CCCCN", "CCCCNN"], ["", "", ""]]

# grid = list(zip(f1_grid, f2_grid, f3_grid, f4_grid))
# print(f"{grid=}")

# print("grid rows:")
# for row in grid:
#     # Remove tuple, flatten into single list
#     flat_row = [item for sublist in row for item in sublist]
#     molecule_grid.append(flat_row)
#     # print(f"{flat_row=}")

# # molecule_grid = molecule_grid + grid

# molecule_grid_flat = [item for sublist in molecule_grid for item in sublist]

# print(f"{molecule_grid_flat=}")
