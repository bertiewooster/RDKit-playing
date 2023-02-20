import pandas
import polars
import pyreadr

# Data from https://github.com/cran/rcdk/blob/master/data/bpdata.RData
# "A dataset containing the structures and associated boiling points for 277 molecules, primarily alkanes and substituted alkanes.""
# Actually has an additional initial column, molecule name
result = pyreadr.read_r('data/bpdata.RData')

# print(result)
df1 = result["bpdata"]
# print(df1)
# df1.to_csv('data/bp.csv')

df = polars.from_pandas(df1, include_index=True)
print(df)
