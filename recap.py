from rdkit import Chem
from rdkit.Chem import Draw, Recap

cisapride = Chem.MolFromSmiles("Clc1cc(c(OC)cc1N)C(=O)NC3CCN(CCCOc2ccc(F)cc2)CC3OC")
# Draw.MolToImage(cisapride)
hierarch = Recap.RecapDecompose(cisapride)
ks = hierarch.children.keys()
ks = sorted(ks)
# print(f"{ks=}")

# all_children = hierarch.GetAllChildren()
# print(f"{all_children=}")

generations = [hierarch.smiles]

level = hierarch
level_counter = 0


def get_children(base_node, level, cols=0):
    for smiles, node in base_node.children.items():
        spaces = " " * level
        print(f"{spaces}{smiles=}")
        children = node.children.keys()
        children = sorted(children)
        if len(children) > 0:
            get_children(node, level + 1, cols)
            # print(f"{children=}")
        else:
            cols += 1
            # print(f"{cols=}")
            # return cols


cols = get_children(hierarch, 0, 0)
# print(f"{cols=}")


class Node:
    def __init__(self, data):

        self.left = None
        self.right = None
        self.data = data

    # Insert Node
    def insert(self, data):

        if self.data:
            # if data < self.data:
            if self.left is None:
                self.left = Node(data)
            else:
                self.left.insert(data)
        # elif data > self.data:
        #     if self.right is None:
        #         self.right = Node(data)
        #     else:
        #         self.right.insert(data)
        else:
            self.data = data
        return self.left

    # Print the Tree
    def PrintTree(self):
        if self.left:
            self.left.PrintTree()
        print(self.data),
        if self.right:
            self.right.PrintTree()

    # Preorder traversal
    # Root -> Left -> Right
    def PreorderTraversal(self, root):
        res = []
        if root:
            res.append(root.data)
            res = res + self.PreorderTraversal(root.left)
            res = res + self.PreorderTraversal(root.right)
        return res


root = Node("cisapride")
f1 = root.insert("f1")
f2 = root.insert("f2")
f21 = f2.insert("f21")
f22 = f2.insert("f22")
f23 = f2.insert("f23")
f26 = f2.insert("f24")
f3 = root.insert("f3")
f4 = root.insert("f4")
f41 = f4.insert("f41")
f42 = f4.insert("f42")
f5 = root.insert("f5")
f6 = root.insert("f6")
f61 = f6.insert("f61")
f62 = f6.insert("f62")

print(root.PreorderTraversal(root))
