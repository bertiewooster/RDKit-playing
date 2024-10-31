import polars as pl

df = pl.DataFrame(
    {
        "a": [1, 2, 3],
        "b": ["x", "y", "z"],
    }
)
df = df.with_columns(pl.all().reverse().name.prefix("flat_"))
print(df)

df = pl.DataFrame({
    'x': [11, 22],
})

def uses_object(x):
    r = list(range(0, x))
    c10 = r.count(10)
    c12 = r.count(12)
    return dict(count_of_10=c10, count_of_12=c12)

df = df.with_columns(
   count = pl.col('x').map_elements(uses_object)
).unnest('count')


print(df)

df = pl.DataFrame({
    'x': [11, 22],
})

df = df.with_columns(
    cnt=pl.col('x').map_elements(uses_object)
).unnest('cnt')

print(df)

# print(df)
