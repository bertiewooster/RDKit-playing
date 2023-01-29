# Source: https://realpython.com/python-async-features/#asynchronous-non-blocking-http-calls
import asyncio

import aiohttp
from codetiming import Timer


async def is_commercially_available(name, work_queue):
    timer = Timer(text=f"Task {name} elapsed time: {{:.2f}}")
    async with aiohttp.ClientSession() as session:
        while not work_queue.empty():
            smiles = await work_queue.get()
            print(f"Task {name} getting URL")
            timer.start()
            get_cid_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
            async with session.get(get_cid_URL, ssl=False) as response:
                get_cid_response = await response.text()
                cid = get_cid_response.strip("\n")
            if cid == '0':
                print(f"Task {name}: No chemical in PubChem")
                return False

            compound_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/XML?heading=Chemical-Vendors"

            async with session.get(compound_url, ssl=False) as response:
                compound_vendors_response = await response.text()

            print(f"  Task {name} {smiles} cid: {cid}")
            print(f"  Task {name} {smiles} {compound_url=}, response = {compound_vendors_response[:70]}")
            timer.stop()

            if "<Message>No data found</Message>" in compound_vendors_response:
                return False
            else:
                return True


async def main():
    """
    This is the main entry point for the program
    """
    # Create the queue of work
    work_queue = asyncio.Queue()

    # Put some work in the queue
    smiles_list = [
        # In PubChem and is commercially available. Should return True.
        "CCCC",
        # In PubChem, but not commercially available. Should return False.
        "C3CCCC(C2CCCC(C1CCCCC1)CC2)CCC3",
        # Not in PubChem. Should return False.
        "CCCCCCCCCCCCCCCCCCCCCCc1ccccc1CCCCCCCCCCCC",
    ]
    for smiles in smiles_list:
        await work_queue.put(smiles)

    # Run the tasks
    with Timer(text="\nTotal elapsed time: {:.2f}"):
        tasks = [is_commercially_available(smiles, work_queue) for smiles in smiles_list]
        values = await asyncio.gather(*tasks)
    
    for index, value in enumerate(values):
        print(f"{index}: {value}")

if __name__ == "__main__":
    asyncio.run(main())
