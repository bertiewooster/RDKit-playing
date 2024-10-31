import polars as pl
from rdkit import Chem

smls = ["CC"] * 10

data = {"sml": smls}

df = pl.DataFrame(data)

# df = df.with_columns(
#     [
#         pl.col("sml").map_elements(Chem.MolFromSmiles).alias("mol_map_elements"),
#     ]
# )

df = df.with_columns(
    [
        pl.col("sml"). map_elements(lambda s: list(Chem.MolFromSmiles(s))).alias("mol_map_elements"),
    ]
)
pass