import datetime
from multiprocessing import Pool, cpu_count

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import PandasTools

# count = 1_000_000
# smls = ["CC"] * count

# start = datetime.datetime.now()

# mols = [Chem.MolFromSmiles(sml) for sml in smls]

# print(f"Converting {count} sml to mol took {datetime.datetime.now() - start}s")

def process_smiles(sml):
    return Chem.MolFromSmiles(sml)

if __name__ == '__main__':
    # count = 1_000_000
    count = 2

    smls = np.full(count, "CC")

    start = datetime.datetime.now()

    # Convert the array of SMILES strings into a pandas DataFrame
    df = PandasTools.SmilesToMol(smls, smilesCol='SMILES', nameCol='Molecule')

    end = datetime.datetime.now()
    print("Time taken using np.frompyfunc:", end - start)

    # smls = ["CC"] * count

    # start = datetime.datetime.now()

    # Determine the number of CPU cores available
    # num_cores = cpu_count()

    # Create a Pool of worker processes
    # with Pool(num_cores) as pool:
        # Distribute the workload across the worker processes
        # mols = pool.map(process_smiles, smls)

    # end = datetime.datetime.now()
    # print("Time taken:", end - start)

# data = {"sml": smls}
# df = pl.DataFrame(data)

# start = datetime.datetime.now()
# df = df.with_columns(
#     [
#         pl.col("sml").map_elements(lambda s: Chem.MolFromSmiles(s)).alias("mol_map_elements"),
#     ]
# )
# print(f"Converting {count} sml to mol took {datetime.datetime.now() - start}s")

# start = datetime.datetime.now()
# df = df.with_columns(
#     [
#         pl.col("sml").str.len_bytes().alias("len"),
#     ]
# )
# print(f"len(sml) for {count} sml took {datetime.datetime.now() - start}s")

# df = df.with_columns(
#     [
#         pl.col("sml"). map_batches(lambda s: Chem.MolFromSmiles(s)).alias("mol_map_batches"),
#     ]
# )

# print(df)
