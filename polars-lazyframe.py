import numpy as np
import polars as pl

print(pl.__version__)

num_rows = 5000
rng = np.random.default_rng(seed=7)

buildings = {
     "sqft": rng.exponential(scale=1000, size=num_rows),
     "price": rng.exponential(scale=100_000, size=num_rows),
     "year": rng.integers(low=1995, high=2023, size=num_rows),
     "building_type": rng.choice(["A", "B", "C"], size=num_rows),
  }
buildings_lazy = pl.LazyFrame(buildings)
# print(buildings_lazy)

lazy_query = (
    buildings_lazy
    .with_columns(
        (pl.col("price") / pl.col("sqft")).alias("price_per_sqft")
    )
    .filter(pl.col("price_per_sqft") > 100)
    .filter(pl.col("year") < 2010)
 )
# print(lazy_query)
# lazy_query.show_graph()
# print(lazy_query.explain())

lazy_query = (
    buildings_lazy
    .with_columns(
        (pl.col("price") / pl.col("sqft")).alias("price_per_sqft")
    )
    .filter(pl.col("price_per_sqft") > 100)
    .filter(pl.col("year") < 2010)
 )

r = (
    lazy_query
    .collect()
    .select(pl.col(["price_per_sqft", "year"]))
)
# print(r)
# print(r.describe())

lazy_car_data = pl.scan_csv("electric_cars.csv")
print(lazy_car_data)
print(dict(lazy_car_data.schema))

lazy_car_query = (
    lazy_car_data
    .filter((pl.col("Model Year") >= 2018))
    .filter(
        pl.col("Electric Vehicle Type") == "Battery Electric Vehicle (BEV)"
    )
    .groupby(["State", "Make"])
    .agg(
        pl.mean("Electric Range").alias("Average Electric Range"),
        pl.min("Model Year").alias("Oldest Model Year"),
        pl.count().alias("Number of Cars"),
    )
    .filter(pl.col("Average Electric Range") > 0)
    .filter(pl.col("Number of Cars") > 5)
    .sort(pl.col("Number of Cars"), descending=True)
)

r = lazy_car_query.collect()
print(r)

data = pl.DataFrame({
    "A": [1, 2, 3, 4, 5],
    "B": [6, 7, 8, 9, 10],
})

data.write_csv("data.csv")
data.write_ndjson("data.json")
data.write_parquet("data.parquet")

data_csv = pl.read_csv("data.csv")
data_csv_lazy = pl.scan_csv("data.csv")
print(data_csv_lazy.schema)

data_json = pl.read_ndjson("data.json")
data_json_lazy = pl.scan_ndjson("data.json")
print(data_json_lazy.schema)


data_parquet = pl.read_parquet("data.parquet")
data_parquet_lazy = pl.scan_parquet("data.parquet")
print(dict(data_parquet_lazy.schema))


