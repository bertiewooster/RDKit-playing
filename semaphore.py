import asyncio

import aiohttp
from codetiming import Timer


class Reactant():
    """Define a reactant's commercial availability."""
        
    def __init__(self, smiles: str):
        """
        Construct a Reactant object to store the commercial availability of a reactant

        :param smiles: SMILES string representing a molecule
        """
        self.smiles = smiles
        # self.mol = Chem.MolFromSmiles(smiles)

        @property
        def in_pubchem(self):
            return self._in_pubchem

        @in_pubchem.setter
        def in_pubchem(self, value: bool):
            """:param value: whether molecule is in PubChem"""
            self._in_pubchem = value

        @property
        def cid(self):
            return self._cid

        @cid.setter
        def cid(self, value: int):
            """:param value: PubChem CID (identifier) for molecule"""
            self._cid = value

        @property
        def commercially_available(self):
            return self._commercially_available

        @commercially_available.setter
        def commercially_available(self, value: bool):
            """:param value: whether molecule is commercially available, per PubChem"""
            self._commercially_available = value

        @property
        def pubchem_page(self):
            return self._pubchem_page

        @pubchem_page.setter
        def pubchem_page(self, value: str):
            """:param value: URL or PubChem page for molecule"""
            self._pubchem_page = value

    def __str__(self):
        str_print = f"Reactant SMILES: {self.smiles}"
        str_print += f", in_pubchem: {self.in_pubchem}"
        if hasattr(self, "cid"):
            str_print += f", cid: {self.cid}"
        if hasattr(self, "commercially_available"):
            str_print += f", commercially_available: {self.commercially_available}"
        if hasattr(self, "pubchem_page"):
            str_print += f", pubchem_page: {self.pubchem_page}"
        return str_print

async def is_commercially_available(smiles):
    """
    Asynchronously check the availability of a queue of SMILES strings (chemicals) in PubChem
    Based on https://realpython.com/python-async-features/#asynchronous-non-blocking-http-calls

    :param reactant_smiles: A SMILES string (representing a molecule)
    :returns: Class Reactant object with information from PubChem
    """
    print(f'Calling API(s) for {smiles}')

    async with aiohttp.ClientSession() as session:

        # Create Reactant object, which will be populated during this function
        reactant = Reactant(smiles)

        timer = Timer(text=f"{{:.2f}}s for {smiles} PubChem API call(s)")

        timer.start()

        # Find the PubChem identifier (CID) for this SMILES string
        get_cid_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
        try:
            async with session.get(get_cid_URL, ssl=False) as response:
                get_cid_response = await response.text()
        except:
            raise ConnectionError
        cid_str = get_cid_response.strip("\n")

        try:
            cid = int(cid_str)
        except ValueError:
            cid = 0

        if cid == 0:
            reactant.in_pubchem = False
            reactant.commercially_available = False
            timer.stop()
            return reactant
        else:
            reactant.cid = cid
            reactant.in_pubchem = True
            reactant.pubchem_page = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

        # Get the compound's PubChem page
        compound_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/XML?heading=Chemical-Vendors"

        async with session.get(compound_url, ssl=False) as response:
            compound_vendors_response = await response.text()

        timer.stop()

        if "<Message>No data found</Message>" in compound_vendors_response:
            reactant.commercially_available = False
        else:
            reactant.commercially_available = True
        print(f"{str(reactant)=}")   

sem = asyncio.Semaphore(2)

async def safe_calls(smiles):
    async with sem:  # semaphore limits num of simultaneous API calls
        return await is_commercially_available(smiles)


async def check_avail_smiles_set(smiles_list):
    smiles_set = set(smiles_list)
    tasks = [asyncio.ensure_future(safe_calls(smiles)) for smiles in smiles_set]
    await asyncio.gather(*tasks)  # await completion of all API calls


def check_avail_smiles_list(smiles_list):
    with Timer(text="-----\n{:.2f}s total elapsed time for PubChem API calls"):
        loop = asyncio.get_event_loop()
        try:
            loop.run_until_complete(check_avail_smiles_set(smiles_list))
        finally:
            loop.run_until_complete(loop.shutdown_asyncgens())
            loop.close()

if __name__ ==  '__main__':
    smiles_list = ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"]
    check_avail_smiles_list(smiles_list)