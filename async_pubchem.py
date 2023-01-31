# Source: https://realpython.com/python-async-features/#asynchronous-non-blocking-http-calls
import asyncio

import aiohttp
import rdkit
from codetiming import Timer
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdChemReactions


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

class Reaction():
    """Define a reaction's target, type, reactants, and commercial of reactants."""

    def __init__(self, target, reaction_smarts, name):
        """
        Construct a Reaction object 

        :param input_file: Name of the file to read the Constraint from (string)
        :param output_file: Name of the file to write the valid points to
        :param n_results_str: Number of valid points to find, as a string until determine that it's an integer
        """
        self.target = target
        self.reaction_smarts = reaction_smarts
        self.name = name
        self.target_mol = Chem.MolFromSmiles(self.target)
        self.reactants = dict()

        @property
        def reactants(self):
            return self._reactants

        @reactants.setter
        def reactants(self, value):
            self._reactants = value

        @property
        def reactants_commercially_available(self):
            return self._reactants_commercially_available

    def tally_all_reactants_commercially_available(self):
        for reactant in self.reactants:
            if not self.reactants[reactant].commercially_available:
                self._reactants_commercially_available = False
                return False
        self._reactants_commercially_available = True
        return True

async def is_commercially_available(smiles, work_queue):
    async with aiohttp.ClientSession() as session:
        while not work_queue.empty():
            smiles_again = await work_queue.get()
            print(f"{smiles=} {smiles_again=}")
            # Create Reactant object, which will be populated during this function
            reactant = Reactant(smiles)

            timer = Timer(text=f"{smiles} PubChem API call(s) elapsed time: {{:.2f}}")

            # print(f"Task {smiles} getting URL")
            timer.start()
            get_cid_URL = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"
            async with session.get(get_cid_URL, ssl=False) as response:
                get_cid_response = await response.text()
                cid = get_cid_response.strip("\n")
            if cid == '0':
                # print(f"Task {smiles}: No chemical in PubChem")
                reactant.in_pubchem = False
                reactant.commercially_available = False
                timer.stop()
                return reactant
            else:
                reactant.in_pubchem = True
                reactant.pubchem_page = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

            compound_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/XML?heading=Chemical-Vendors"

            async with session.get(compound_url, ssl=False) as response:
                compound_vendors_response = await response.text()

            # print(f"  Task {smiles} cid: {cid}")
            # print(f"  Task {smiles} {compound_url=}, response = {compound_vendors_response[:70]}")
            timer.stop()

            if "<Message>No data found</Message>" in compound_vendors_response:
                reactant.commercially_available = False
            else:
                reactant.commercially_available = True
    return reactant

async def check_avail_several(smiles_set: set[str]) -> dict[str, object]:
    """
    Asynchronously check the availability of several SMILES strings (chemicals) in PubChem
    """
    # Create the queue of work
    work_queue = asyncio.Queue()

    # Add each SMILES (molecule) to the work queue
    for smiles in smiles_set:
        await work_queue.put(smiles)

    # Determine commercial availability of each reactant
    with Timer(text="\n--Total elapsed time: {:.2f}"):
        tasks = [is_commercially_available(smiles, work_queue) for smiles in smiles_set]
        reactants = await asyncio.gather(*tasks)
    
    # Put reactants in dictionary of SMILES:Reaction object
    smiles_avail = dict()
    for index, reactant in enumerate(reactants):
        # print(f"Reactant {index}: {reactant}, {reactant.smiles}, in_pubchem: {reactant.in_pubchem}, commercially_available: {reactant.commercially_available}")
        smiles_avail[reactant.smiles] = reactant
    
    return smiles_avail

def reverse_reaction(rxn_fwd):
    """rxn_fwd: rdkit.Chem.rdChemReactions.ChemicalReaction"""
    rxn_rev = Chem.ChemicalReaction()
    for i in range(rxn_fwd.GetNumReactantTemplates()):
        rxn_rev.AddProductTemplate(rxn_fwd.GetReactantTemplate(i))
    for i in range(rxn_fwd.GetNumProductTemplates()):
        rxn_rev.AddReactantTemplate(rxn_fwd.GetProductTemplate(i))
    rxn_rev.Initialize()
    return rxn_rev

async def check_reactions(target_reaction_list: list[list[str, str]]):
    """
    target: SMILES
    reaction_smarts: SMARTS
    """

    all_reactants_set = set()

    # List of Reaction objects
    reactions = []
    for target_reaction in target_reaction_list:
        reaction = Reaction(target_reaction[0], target_reaction[1], target_reaction[2])

        # Create forward reaction
        rxn_fwd = Chem.ReactionFromSmarts(reaction.reaction_smarts)

        # Reverse reaction
        rxn_rev = reverse_reaction(rxn_fwd)

        # Run reverse reaction to determine starting materials
        
        reactants = rxn_rev.RunReactants([reaction.target_mol])[0]
        for reactant in reactants:
            reactant_smiles = Chem.MolToSmiles(reactant)
            reaction.reactants[reactant_smiles] = None

        # Add starting materials to set of starting materials
        for reactant in reaction.reactants:
            # reactant_smiles = Chem.MolToSmiles(reactant)
            all_reactants_set.add(reactant)
        
        reactions.append(reaction)

    # Check commercial availability of set of starting materials
    smiles_avail = await check_avail_several(all_reactants_set)

    print("Are starting materials commercially available for reaction -> target:")

    # [[reaction object 0, all reactants available], [reaction object 1, all reactants available],]
    reaction_reactants_avail = [[]]
    for reaction in reactions:
        # Add information to Reaction objects
        for this_reactant in reaction.reactants:
            # Set value for key of reactant SMILES, to reactant object from smiles_avail
            # this_reactant_smiles = Chem.MolToSmiles(this_reactant)
            reaction.reactants[this_reactant] = smiles_avail[this_reactant]
        reaction.tally_all_reactants_commercially_available()

        # Return results
        reaction_reactants_avail.append([reaction, reaction._reactants_commercially_available])
        print(f"{reaction.name} -> {reaction.target}: {reaction._reactants_commercially_available}")

    return reaction_reactants_avail

async def check_rxns():
    bicyclic_target = "OCN1C2CC(C=C2)C1CC1NCCc2ccccc12"
    aniline_target = "c1ccc(NCC2NCCc3ccccc32)cc1"
    cyclobutyl_target = "FC1(F)CC(C2NCCc3ccccc32)C1"

    pictet_spengler_rxn = '[cH1:1]1:[c:2](-[CH2:7]-[CH2:8]-[NH2:9]):[c:3]:[c:4]:[c:5]:[c:6]:1.[#6:11]-[CH1;R0:10]=[OD1]>>[c:1]12:[c:2](-[CH2:7]-[CH2:8]-[NH1:9]-[C:10]-2(-[#6:11])):[c:3]:[c:4]:[c:5]:[c:6]:1'

    # Reaction format: [target (SMILES), reaction_smarts]
    rxn1 = [bicyclic_target, pictet_spengler_rxn, "Pictet Spengler"]
    rxn2 = [aniline_target, pictet_spengler_rxn, "Pictet Spengler"]
    rxn3 = [cyclobutyl_target, pictet_spengler_rxn, "Pictet Spengler"]

    rxns = [rxn1, rxn2, rxn3]
    await check_reactions(rxns)

if __name__ == "__main__":
    asyncio.run(check_rxns())
