# Import library
from Bio.Seq import Seq
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
print(lines)

# Length of Dataset is 608.


# Delete unnecessary variables:
del header, file, path, path_parent, f
