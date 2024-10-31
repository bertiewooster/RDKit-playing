import polars as pl

s = pl.Series([1,2,3])
print(s)

s = s.append(pl.Series([99, 11]))
print(f"{s=}")

print(f"{s.get_chunks()=}")
