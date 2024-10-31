import datetime
from functools import wraps
from time import time

import polars as pl
from rdkit import Chem
from rdkit.Chem import PandasTools

from polars_xlsx import xlsx_from_polarsdf, xlsx_from_polarsdf_idiomatic


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f'Function {f.__name__} took {te-ts:2.4f} seconds') 
        return result
    return wrap

sml = 'CO'
mol = Chem.MolFromSmiles(sml)
df = pl.DataFrame({
    'sml': sml,
    'mol': mol,
    })

mult = 5
smls_many = [sml] * mult
mols_many = [mol] * mult
ints_many = [1] * mult
floats_many = [1.01] * mult
times_many = [datetime.datetime.now()] * mult

df_many = pl.DataFrame({
    'sml': smls_many,
    'mol': mols_many,
    'int': ints_many,
    'float': floats_many,
    'datetime': times_many,
    })

df_many_html = df_many._repr_html_()

with open('polars_output.html', 'w') as f:
    f.write(df_many_html)


# xlsx_from_polarsdf(
#     df=df_many,
#     outFile="xlsx_from_polars.xlsx",
#     molCol="mol"
# )

xlsx_from_polarsdf_idiomatic(
    df_many,
    outFile="xlsx_from_polars_idiomatic.xlsx",
    molCol='mol',
)

# Pandas

df_many_pd = df_many.to_pandas()

@timing
def run_pandas():
    PandasTools.SaveXlsxFromFrame(
        frame = df_many_pd,
        outFile="xlsx_from_pandas.xlsx",
        molCol="mol"
    )

with open('pandas_output.html', 'w') as f:
    f.write(df_many_pd.to_html())

# run_pandas()

