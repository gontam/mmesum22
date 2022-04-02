# Import library
from Bio import Entrez
import pandas as pd
import os
import collections

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

# See if same values are in list and delete them:
list_double = [item for item, count in collections.Counter(listNew).items() if count > 1]

# Gather information:
i = 0
id_list = []
while i <= len(list_double):
    if i == len(list_double):
        break
    else:
        list_double[i]
        Entrez.email = "me21x506@technikum-wien.at"
        handle = Entrez.esearch(db="nuccore", term=list_double[i])
        record = Entrez.read(handle)
        id_list.append(record["IdList"])
        i += 1

# Get Information out from id's and save them in csv file:
# Had no further time because a friend from germany was visiting me in Vienna. Therefore I could not make the last step.
#while i <= len(id_list):


# Delete unnecessary variables:
del file, lines, lines_data_set, i, lengthFile, listNew, list_double, path, path_parent, x
