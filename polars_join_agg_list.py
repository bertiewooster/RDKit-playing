import polars as pl

# Create DataFrame df1
df1 = pl.DataFrame({
    'id': [1, 2, 3]
})

print(df1)

# Create DataFrame df2
df2 = pl.DataFrame({
    'id': [1, 1],
    'name': ['Jim', 'Mary']
})

print(df2)

# Perform a left join on 'id' column and aggregate 'name' column into a list
result = df1.join(df2, on='id', how='left').groupby('id').agg('name')

print(result)
