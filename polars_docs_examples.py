# From https://docs.pola.rs/user-guide/io/multiple/#dealing-with-multiple-files

import glob
import os

import polars as pl

df = pl.DataFrame({"foo": [1, 2, 3], "bar": [None, "ham", "spam"]})

os.makedirs("docs/data", exist_ok=True)

for i in range(5):
    df.write_csv(f"docs/data/my_many_files_{i}.csv")

df = pl.read_csv("docs/data/my_many_files_*.csv")
print(df)

# pl.scan_csv("docs/data/my_many_files_*.csv").show_graph()

# queries = []
# for file in glob.glob("docs/data/my_many_files_*.csv"):
#     q = pl.scan_csv(file).group_by("bar").agg(pl.len(), pl.sum("foo"))
#     queries.append(q)

# dataframes = pl.collect_all(queries)
# print(dataframes)

import numpy as np
import pandas as pd
import polars as pl

polars_data = pl.DataFrame({
    "A": [1, 2, 3, 4, 5],
    "B": [6, 7, 8, 9, 10]
})

pandas_data = pd.DataFrame({
    "A": [1, 2, 3, 4, 5],
    "B": [6, 7, 8, 9, 10]
})

numpy_data = np.array([
    [1, 2, 3, 4, 5],
    [6, 7, 8, 9, 10]
]).T

fp = pl.from_pandas(pandas_data)
print(fp)

fn = pl.from_numpy(numpy_data, schema={"A": pl.Int64, "B": pl.Int64})
print(fn)

tp = polars_data.to_pandas()
print(f"{tp=}")

tn = polars_data.to_numpy()
print(f"{tn=}")
