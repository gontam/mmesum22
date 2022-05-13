# Handed in by Oliver Jovanović
# Import Libraries:
from Bio import Entrez
from Bio import SeqIO
import pandas as pd


# Chose gene 'APOL1' from Assignment 1:
# Have a look how many results are there for 'APOL1' at NCBI:
Entrez.email = 'me21x506@technikum-wien.at'
handle = Entrez.egquery(term='APOL1')
record_test = Entrez.read(handle)
for row in record_test['eGQueryResult']:
    if row['DbName'] == 'nuccore':
        print(row['Count'])
handle.close()
# 7062 Results!

# Important Information:
# APOL1 can only be found in primal species. It will be differentiated between the primal species:
# Have this gene: - Human, - Gorilla, - Baboon
# Have a pseudo gene: - Orangutan, - Macaque

# Download five different entries from five different primate species:
# Human -- 5:
# Fetch ID:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Human', retmax=5, idtype='acc')
record_Human = Entrez.read(handle)
id_Human = record_Human["IdList"][0]

# Get Other information:
handle = Entrez.efetch(db='nucleotide', id=id_Human, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
for ref in record.annotations['references']:
    title = ref.title
    y = 0
    if y == 0:
        break
title_Human = title
organism_Human = record.annotations['organism']
length_Human = record.__len__()
lst_Human = [id_Human, title_Human, organism_Human, length_Human]

# Gorilla -- 2:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Gorilla', retmax=5, idtype='acc')
record_Gorilla = Entrez.read(handle)
id_Gorilla = record_Gorilla["IdList"][0]

# Get Other information:
handle = Entrez.efetch(db='nucleotide', id=id_Gorilla, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
title_Gorilla = record.annotations['comment'] # No References (Titles).
organism_Gorilla = record.annotations['organism']
length_Gorilla = record.__len__()
lst_Gorilla = [id_Gorilla, title_Gorilla, organism_Gorilla, length_Gorilla]

# Baboon -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Baboon', retmax=5, idtype='acc')
record_Baboon = Entrez.read(handle)
id_Baboon = record_Baboon["IdList"][0]


# Get Other information:
handle = Entrez.efetch(db='nucleotide', id=id_Baboon, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
for ref in record.annotations['references']:
    title = ref.title
    y = 0
    if y == 0:
        break
title_Baboon = title
organism_Baboon = record.annotations['organism']
length_Baboon = record.__len__()
lst_Baboon = [id_Baboon, title_Baboon, organism_Baboon, length_Baboon]

# Pseudo gene:
# Orangutan -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Orangutan', retmax=5, idtype='acc')
record_Orangutan = Entrez.read(handle)
id_Orangutan = record_Orangutan["IdList"][0]

# Get Other information:
handle = Entrez.efetch(db='nucleotide', id=id_Orangutan, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
title_Orangutan = record.annotations['comment']
organism_Orangutan = record.annotations['organism']
length_Orangutan = record.__len__()
lst_Orangutan = [id_Orangutan, title_Orangutan, organism_Orangutan, length_Orangutan]

# Macaque -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Macaque', retmax=5, idtype='acc')
record_Macaque = Entrez.read(handle)
id_Macaque = record_Macaque["IdList"][0]

# Get Other information:
handle = Entrez.efetch(db='nucleotide', id=id_Macaque, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
for ref in record.annotations['references']:
    title = ref.title
    y = 0
    if y == 0:
        break
title_Macaque = title
organism_Macaque = record.annotations['organism']
length_Macaque = record.__len__()
lst_Macaque = [id_Macaque, title_Macaque, organism_Macaque, length_Macaque]

# Create Dataframe with all information:
df = pd.DataFrame(list(zip(lst_Human, lst_Gorilla, lst_Baboon, lst_Orangutan, lst_Macaque)),
                  columns=['Human', 'Gorilla', 'Baboon', 'Orangutan', 'Macaque'])
display(df)

# Delete unnecessary variables:
del handle, record, id_Baboon, lst_Baboon, id_Gorilla, lst_Gorilla, id_Human, lst_Human, id_Macaque, lst_Macaque
del id_Orangutan, lst_Orangutan, length_Orangutan, length_Macaque, length_Human, length_Baboon, length_Gorilla
del organism_Macaque, organism_Baboon, organism_Human, organism_Gorilla, organism_Orangutan
del title, y, title_Macaque, title_Orangutan, title_Baboon, title_Human, title_Gorilla
del record_Baboon, record_test, record_Macaque, record_Orangutan, record_Human, record_Gorilla
del ref, row
