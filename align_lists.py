def align_lists(
        list1,
        list2,
        filler = "",
):
    """Align list2 under list1, using filler for items in list1 not in list2, 
    for example:
    
    input:
    list1: a,  b,  c
    list2: a,  d,  c,  e,  f

    output:
    list2: a, "",  c,  d,  e,  f
    """
    # Determine index for each element in list2:
    list_2_indices = []
    list2_index_max = len(list1) - 1
    for _, list2_item in enumerate(list2):
        try:
            list2_index = list1.index(list2_item)
        except:
            list2_index_max += 1
            list2_index = list2_index_max
        print(f"{list2_item} will be placed in column {list2_index}")
        list_2_indices.append(list2_index)
    # Create sorted (with blanks) list2 by that ordering list
    list2_aligned = [filler] * (max(list_2_indices) + 1)
    for i, list_2_col in enumerate(list_2_indices):
        list2_aligned[list_2_col] = list2[i]
    return list2_aligned

list1 = ["a", "c", "e"]
list2 = ["a", "d", "c", "e", "f"]
align_lists(list1, list2)