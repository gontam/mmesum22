# Handed in by Oliver Jovanović
# Import Libraries:
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import AlignIO
import pandas as pd
import os
from Bio.Phylo.TreeConstruction import DistanceCalculator

# Chose gene 'APOL1' from Assignment 1:
# Have a look how many results are there for 'APOL1' at NCBI:
Entrez.email = 'me21x506@technikum-wien.at'
handle = Entrez.egquery(term='APOL1')
record_test = Entrez.read(handle)
for row in record_test['eGQueryResult']:
    if row['DbName'] == 'nuccore':
        print("Total number of entries:", row['Count'])
handle.close()

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

# Get further information:
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
GC_Human = GC(record.seq)
Seq_Human = record.seq
lst_Human = [id_Human, title_Human, organism_Human, length_Human, GC_Human]

# Writing sequence file:
Entrez.email = "me21x506@technikum-wien.at"
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=id_Human) as handle:
    seq_record = SeqIO.read(handle, 'fasta')
SeqIO.write(seq_record, 'Human.faa', 'fasta')

# Gorilla -- 2:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Gorilla', retmax=5, idtype='acc')
record_Gorilla = Entrez.read(handle)
id_Gorilla = record_Gorilla["IdList"][0]

# Get further information:
handle = Entrez.efetch(db='nucleotide', id=id_Gorilla, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
title_Gorilla = record.annotations['comment']  # No References (Titles).
organism_Gorilla = record.annotations['organism']
length_Gorilla = record.__len__()
GC_Gorilla = GC(record.seq)
Seq_Gorilla = record.seq
lst_Gorilla = [id_Gorilla, title_Gorilla, organism_Gorilla, length_Gorilla, GC_Gorilla]

# Writing sequence file:
Entrez.email = "me21x506@technikum-wien.at"
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=id_Gorilla) as handle:
    seq_record = SeqIO.read(handle, 'fasta')
SeqIO.write(seq_record, 'Gorilla.faa', 'fasta')

# Baboon -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Baboon', retmax=5, idtype='acc')
record_Baboon = Entrez.read(handle)
id_Baboon = record_Baboon["IdList"][0]


# Get further information:
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
GC_Baboon = GC(record.seq)
Seq_Baboon = record.seq
lst_Baboon = [id_Baboon, title_Baboon, organism_Baboon, length_Baboon, GC_Baboon]

# Writing sequence file:
Entrez.email = "me21x506@technikum-wien.at"
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=id_Baboon) as handle:
    seq_record = SeqIO.read(handle, 'fasta')
SeqIO.write(seq_record, 'Baboon.faa', 'fasta')

# Pseudo gene:
# Orangutan -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Orangutan', retmax=5, idtype='acc')
record_Orangutan = Entrez.read(handle)
id_Orangutan = record_Orangutan["IdList"][0]

# Get further information:
handle = Entrez.efetch(db='nucleotide', id=id_Orangutan, rettype='gb', retmode='text')
record = SeqIO.read(handle, 'genbank')
title_Orangutan = record.annotations['comment']  # No Reference (Titles).
organism_Orangutan = record.annotations['organism']
length_Orangutan = record.__len__()
GC_Orangutan = GC(record.seq)
Seq_Orangutan = record.seq
lst_Orangutan = [id_Orangutan, title_Orangutan, organism_Orangutan, length_Orangutan, GC_Orangutan]

# Writing sequence file:
Entrez.email = "me21x506@technikum-wien.at"
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=id_Orangutan) as handle:
    seq_record = SeqIO.read(handle, 'fasta')
SeqIO.write(seq_record, 'Orangutan.faa', 'fasta')

# Macaque -- 5:
handle = Entrez.esearch(db='nucleotide', term='APOL1 Macaque', retmax=5, idtype='acc')
record_Macaque = Entrez.read(handle)
id_Macaque = record_Macaque["IdList"][0]

# Get further information:
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
GC_Macaque = GC(record.seq)
Seq_Macaque = record.seq
lst_Macaque = [id_Macaque, title_Macaque, organism_Macaque, length_Macaque, GC_Macaque]

# Writing sequence file:
Entrez.email = "me21x506@technikum-wien.at"
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=id_Macaque) as handle:
    seq_record = SeqIO.read(handle, 'fasta')
SeqIO.write(seq_record, 'Macaque.faa', 'fasta')

# Create Dataframe with all information:
lst = ['Accession Number', 'Title', 'Organism', 'Length of Sequence', 'GC Percentage in %']
df = pd.DataFrame(list(zip(lst, lst_Human, lst_Gorilla, lst_Baboon, lst_Orangutan, lst_Macaque)),
                  columns=['Information', 'Human', 'Gorilla', 'Baboon', 'Orangutan', 'Macaque'])
print(df)

# Sequence alignment implementation:
alignment_Human = AlignIO.read(open('Human.faa'), 'fasta')
alignment_Gorilla = AlignIO.read(open('Gorilla.faa'), 'fasta')
alignment_Baboon = AlignIO.read(open('Baboon.faa'), 'fasta')
alignment_Orangutan = AlignIO.read(open('Orangutan.faa'), 'fasta')
alignment_Macaque = AlignIO.read(open('Macaque.faa'), 'fasta')

alignments_all = [alignment_Human, alignment_Gorilla, alignment_Baboon, alignment_Orangutan, alignment_Macaque]

# Rewrite alignments in PHYLIP format:
AlignIO.write(alignments_all, 'alignments.phy', 'phylip')

# Phylogenetic tree implementation:
# I could not implement the Phylogenetic tree. I do not how to correctly format the file.

# Delete unnecessary variables:
del handle, record, id_Baboon, lst_Baboon, id_Gorilla, lst_Gorilla, id_Human, lst_Human, id_Macaque, lst_Macaque
del id_Orangutan, lst_Orangutan, length_Orangutan, length_Macaque, length_Human, length_Baboon, length_Gorilla, lst
del organism_Macaque, organism_Baboon, organism_Human, organism_Gorilla, organism_Orangutan
del title, y, title_Macaque, title_Orangutan, title_Baboon, title_Human, title_Gorilla
del record_Baboon, record_test, record_Macaque, record_Orangutan, record_Human, record_Gorilla
del ref, row, GC_Orangutan, GC_Macaque, GC_Human, GC_Baboon, GC_Gorilla
del Seq_Baboon, Seq_Gorilla, Seq_Human, Seq_Macaque, Seq_Orangutan
