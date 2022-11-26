grid = [['parent', '', ''],
 ['child1', 'child2', '',],
 ['', 'grandchild1', 'grandchild2']]

for col1, col2, col3 in grid:
    print (f"{col1:<20}{col2:<20}{col3:<20}")

f_string = ''
for col_index in range(len(grid[0])):
    f_string += "{row[" + str(col_index) + "]:<20}"

print(f"{f_string=}")

tup = ()
tup += (())

print(f"{tup=}")


for tup in grid:
    print(tup)

for row_index, row in enumerate(grid):
    print(row)
    # print (f"{f_string[row_index]}")
