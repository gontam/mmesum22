# Handed in by Oliver Jovanović
# Import Libraries:
from Bio import Entrez


# Chose gene 'APOL1' from Assignment 1:
# Have a look how many results are there for 'APOL1' at NCBI:
Entrez.email = 'me21x506@technikum-wien.at'
handle = Entrez.egquery(term='APOL1')
record = Entrez.read(handle)
for row in record['eGQueryResult']:
    if row['DbName'] == 'nuccore':
        print(row['Count'])
handle.close()
# 7060 Results!

# Download five different entries:
handle = Entrez.esearch(db='nucleotide', term='APOL1', retmax=5, idtype='acc')
record = Entrez.read(handle)
handle.close()