# Import library
import Bio
from Bio import Entrez
import pandas as pd
import os

# Variables:
lines = []

# Create Dataframe with header and first row with empty values:
header = {'Accession Number': [0], 'Title': [0], 'Organism': [0], 'Sequence length': [0], 'CG %': [0],
          'Protein instability': [0], 'Aromaticity': [0],
          'Isoelectric point': [0]}
dataframe = pd.DataFrame(data=header)

# Read data in:
# Get the Path
path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent)
file = '\encryption.txt'
path = os.getcwd() + file

# Read Data in:
with open(path, 'r') as f:
    lines = f.readlines()
f.close()

# Length of Dataset is 608:
lines_data_set = lines[0:607]
lengthFile = len(lines_data_set)

# Create a list element which we can use for the further filling of the csv:
x = 0
str1 = ""
listNew = []
while x <= lengthFile - 1:
    string = lines_data_set[x]
    string = string.strip('\n')
    listNew += string
    x += 1

# Gather information:
i = 0
while i <= len(listNew):
    Entrez.email = "me21x506@technikum-wien.at"
    handle = Entrez.esearch(db="nuccore", term=listNew[i])
    record = Entrez.read(handle)
    gi_list = record["IdList"]
    i += 1

# Delete unnecessary variables:
del header, file, path, path_parent, f
