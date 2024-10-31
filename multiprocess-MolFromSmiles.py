import datetime
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem


def process_smiles(sml):
    return Chem.MolFromSmiles(sml)

with ProcessPoolExecutor() as executor:
    count = 2
    # count = 1_000_000
    smls = ["CC"] * count

    start = datetime.datetime.now()
    mols = list(executor.map(process_smiles, smls))
    print(f"Converting {count} sml to mol took {datetime.datetime.now() - start}s")

    # multiprocessing.freeze_support()

    # pool = multiprocessing.Pool()
    # mols = pool.map(process_smiles, smls)
    # pool.close()
    # pool.join()
