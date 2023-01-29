# Source: https://realpython.com/python-async-features/#asynchronous-non-blocking-http-calls
import asyncio

import aiohttp
import requests
from codetiming import Timer


async def task(name, work_queue):
    timer = Timer(text=f"Task {name} elapsed time: {{:.2f}}")
    async with aiohttp.ClientSession() as session:
        while not work_queue.empty():
            smiles = await work_queue.get()
            # print(f"Task {name} getting URL: {url}")
            timer.start()
            get_cid_response = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
            async with session.get(get_cid_response, ssl=False) as response:
                get_cid_response = await response.text()
                cid = get_cid_response.strip("\n")
            print(f"  Task {name} {smiles} cid: {cid}")
            timer.stop()

async def main():
    """
    This is the main entry point for the program
    """
    # Create the queue of work
    work_queue = asyncio.Queue()

    # Put some work in the queue
    for smiles in [
        "CCCC",
        "c1ccccc1",
    ]:
        await work_queue.put(smiles)

    # Run the tasks
    with Timer(text="\nTotal elapsed time: {:.2f}"):
        await asyncio.gather(
            asyncio.create_task(task("One", work_queue)),
            asyncio.create_task(task("Two", work_queue)),
        )

if __name__ == "__main__":
    asyncio.run(main())
