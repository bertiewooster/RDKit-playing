import polars as pl

df = pl.DataFrame({
    'id': [1, 1],
    'name': [["Bob"], ["Mary", "Sue"]],
    'name_nested': [[["Bob"], ["Mary", "Sue"]], None],
})
print(df)

# df_desired = pl.DataFrame({
#     'id': [1],
#     'name': [["Bob", "Mary", "Sue"]],
# })
# print(df_desired)

# Perform the aggregation
result_df = df.groupby('id').agg(pl.col("name").explode())

print(result_df)
