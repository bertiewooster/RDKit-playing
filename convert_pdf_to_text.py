import re

import polars as pl
from pypdf import PdfReader

reader = PdfReader("Journal articles/Wiener index ja01193a005.pdf")
number_of_pages = len(reader.pages)
page = reader.pages[1]
text = page.extract_text()

# print(number_of_pages)
# print(text)

# Fix data from tables
# Convert these characters
#   « -> n
#   ^ -> 2
# Delete any line starting with a decimal point or comma

molecules = []
tobss = []

with open("data/wiener_table_III.txt") as f:
    content = f.readlines()

# Show the file contents line by line.
# We added the comma to print single newlines and not double newlines.
# This is because the lines contain the newline character '\n'.
ignore_line_chars = (".", ",")
for line in content:
    if line[0] not in ignore_line_chars:
        # TODO Remove any spaces before "ane "
        end_of_molecule = line.find("ane ")
        no_spaces_in_molecule = line[:end_of_molecule].replace(" ", "")
        molecule_clean = no_spaces_in_molecule.replace("«", "n").replace("^", "2")
        line_clean = molecule_clean + line[end_of_molecule:]

        words = line_clean.split()
        molecule = words[0]
        tobs = words[1]
        molecules.append(molecule)
        tobss.append(tobs)

dataframe = pl.DataFrame({"molecules": molecules,
                          "tobss": tobss})

print(dataframe)

