import pandas
import pyreadr

result = pyreadr.read_r('data/bpdata.RData')

# print(result)
df1 = result["bpdata"]
print(df1)
df1.to_csv('data/bp.csv')