import datetime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem


def process_smiles(sml):
    return Chem.MolFromSmiles(sml)

if __name__ == "__main__":
    multiprocessing.freeze_support()

    # count = 2
    count = 100_000
    smls = ["CC"] * count

    start = datetime.datetime.now()

    with ProcessPoolExecutor() as executor:
        mols = list(executor.map(process_smiles, smls))
    print(f"Converting {count} sml to mol took {datetime.datetime.now() - start}s")
