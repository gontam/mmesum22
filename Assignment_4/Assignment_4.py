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
# 7062 Results!

# Important Information:
# APOL1 can only be found in primal species. It will be differentiated between the primal species:
# Have this gene: - Human, - Gorilla, - Baboon
# Have a pseudo gene: - Orangutan, - Macaque

# Download five different entries from five different primate species:
# Human:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Human', retmax=5, idtype='acc')
record_Human = Entrez.read(handle)

# Gorilla:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Gorilla', retmax=5, idtype='acc')
record_Gorilla = Entrez.read(handle)

# Baboon:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Baboon', retmax=5, idtype='acc')
record_Baboon = Entrez.read(handle)

# Pseudo gene:
# Orangutan:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Baboon', retmax=5, idtype='acc')
record_Orangutan = Entrez.read(handle)

# Macaque:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Macaque', retmax=5, idtype='acc')
record_Macaque = Entrez.read(handle)