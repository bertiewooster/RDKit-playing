import time

import polars as pl

n = 5_000_000
# n = 1_000
# n = 3
a_values = list(range(n))
b_values = list(range(n))

df = pl.DataFrame(
    {
        "a": a_values,
        "b": b_values,
    }
)

start_time = time.time()
series = df.to_struct("nums")
print(f"to_struct took {time.time() - start_time} seconds")
# print(series)

def struct_to_dict(series):
    names = series.struct.fields

    my_dict = {a:b for a, b in zip(series.struct.field(names[0]), series.struct.field(names[1]))}
    return my_dict

start_time = time.time()
result_dict_via_series = struct_to_dict(series)
print(f"two series method took {time.time() - start_time} seconds")

names = series.struct.fields

# Convert series to dictionary
start_time = time.time()
result_dict = dict((row[names[0]], row[names[1]]) for row in series)
print(f"struct method took {time.time() - start_time} seconds")

# print(result_dict)

start_time = time.time()
result_dict_iter_rows = {}
for row in df.iter_rows():
    result_dict_iter_rows[row[0]] = row[1]
print(f"iter_rows method took {time.time() - start_time} seconds")
