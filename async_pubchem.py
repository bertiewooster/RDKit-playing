# Source: https://realpython.com/python-async-features/#asynchronous-non-blocking-http-calls
import asyncio

import aiohttp
from codetiming import Timer
from rdkit import Chem


class Reactant():
    """Define a reactant's commercial availability."""

    def __init__(self, smiles):
        """
        Construct a Reactant object to store the commercial availability of a reactant

        :param input_file: Name of the file to read the Constraint from (string)
        :param output_file: Name of the file to write the valid points to
        :param n_results_str: Number of valid points to find, as a string until determine that it's an integer
        """
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)

        @property
        def in_pubchem(self):
            return self._in_pubchem

        @in_pubchem.setter
        def in_pubchem(self, value):
            self._in_pubchem = value

        @property
        def pubchem_cid(self):
            return self._pubchem_cid

        @pubchem_cid.setter
        def pubchem_cid(self, value):
            self._pubchem_cid = value

        @property
        def commercially_available(self):
            return self._commercially_available

        @commercially_available.setter
        def commercially_available(self, value):
            self._commercially_available = value

        @property
        def pubchem_page(self):
            return self._pubchem_page

        @pubchem_page.setter
        def pubchem_page(self, value):
            self._pubchem_page = value

        def __repr__(self):
            str_print = f"SMILES {self.smiles}"
            str_print += f", in_pubchem = {self.in_pubchem}"
            str_print += f", pubchem_cid = {self.pubchem_cid}"
            str_print += f", in_pubchem = {self.in_pubchem}"
            str_print += f", commercially_available = {self.commercially_available}"
            str_print += f", pubchem_page = {self.pubchem_page}"
            return str_print

        def __str__(self):
            return __repr__(self)

async def is_commercially_available(smiles, work_queue):
    async with aiohttp.ClientSession() as session:
        while not work_queue.empty():
            smiles_again = await work_queue.get()
            print(f"{smiles=} {smiles_again=}")
            # Create Reactant object, which will be populated during this function
            reactant = Reactant(smiles)

            timer = Timer(text=f"Task {smiles} elapsed time: {{:.2f}}")

            print(f"Task {smiles} getting URL")
            timer.start()
            get_cid_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
            async with session.get(get_cid_URL, ssl=False) as response:
                get_cid_response = await response.text()
                cid = get_cid_response.strip("\n")
            if cid == '0':
                print(f"Task {smiles}: No chemical in PubChem")
                reactant.in_pubchem = False
                reactant.commercially_available = False
                return reactant
            else:
                reactant.in_pubchem = True
                reactant.pubchem_page = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

            compound_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/XML?heading=Chemical-Vendors"

            async with session.get(compound_url, ssl=False) as response:
                compound_vendors_response = await response.text()

            print(f"  Task {smiles} cid: {cid}")
            print(f"  Task {smiles} {compound_url=}, response = {compound_vendors_response[:70]}")
            timer.stop()

            if "<Message>No data found</Message>" in compound_vendors_response:
                reactant.commercially_available = False
            else:
                reactant.commercially_available = True
    return reactant

async def check_avail_several(smiles_set: set[str]):
    """
    Asynchronously check the availability of several SMILES strings (chemicals) in PubChem
    """
    # Create the queue of work
    work_queue = asyncio.Queue()

    # Add each SMILES to the work queue
    for smiles in smiles_set:
        await work_queue.put(smiles)

    # Run the tasks
    with Timer(text="\nTotal elapsed time: {:.2f}"):
        tasks = [is_commercially_available(smiles, work_queue) for smiles in smiles_set]
        reactants = await asyncio.gather(*tasks)
    
    # n_reactants = len(reactants)
    for index, reactant in enumerate(reactants):
        print(f"Reactant {index}: {reactant}, {reactant.smiles}, in_pubchem: {reactant.in_pubchem}, commercially_available: {reactant.commercially_available}")

if __name__ == "__main__":
    reactants = [
        # In PubChem and is commercially available. Should return True.
        "CCCC",
        # In PubChem and is commercially available. Should return True.
        "CCCC",
        # In PubChem, but not commercially available. Should return False.
        "C3CCCC(C2CCCC(C1CCCCC1)CC2)CCC3",
        # Not in PubChem. Should return False.
        "CCCCCCCCCCCCCCCCCCCCCCc1ccccc1CCCCCCCCCCCC",
    ]

    # Merge list of starting materials into a set, so only query PubChem once for each reactant
    reactants_set = set(reactants)
    asyncio.run(check_avail_several(reactants_set))
