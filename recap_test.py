import pytest

from recap import NonBinTree

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


@pytest.mark.parametrize(
    "name, tree, grid",
    [
        ("root", root, [['Fe', '', '', '', '', '', '', ''], ["C", "CC", "", "", "", "CCC", "CCCC", ""], ["", 'CCN', 'CCNN', '', 'CCNNN', "", "CCCCN", "CCCCNN"], ["", '', 'CCNNO', 'CCNNOO', '', '', '', '']]),
        ("f1", f1, [['C']]),
        ("f2", f2, [["CC", "", "", ""], ['CCN', 'CCNN', '', 'CCNNN'], ['', 'CCNNO', 'CCNNOO', '']]),
        ("f21", f21, [['CCN']]),
        ("f22", f22, [['CCNN', ''], ['CCNNO', 'CCNNOO']]),
        ("f221", f221, [['CCNNO']]),
        ("f23", f23, [['CCNNN']]),
        ("f3", f3, [['CCC']]),
        ("f4", f4, [['CCCC', ''], ['CCCCN', 'CCCCNN']]),
    ])
def test_constraint_string(
    name, tree, grid
):
    a = name
    assert grid == tree.get_grid()
