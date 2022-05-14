# Assignment 4 - Biopython/Phylogenetics
# Author: Johanna Schachl
# Date: 01.05.2022
# Last change: 14.05.2022

# Basic information of 5 nucleotide sequences of 5 different species are saved
# Sequence alignment is calculated
# phylogenetic tree is visualized

# gene: IL-8 (LC461682)

# import libraries
from Bio import Entrez, SeqIO, Phylo, AlignIO, Seq
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import csv
import os

Entrez.email = "me21m009@technikum-wien.at"

# variables
geneName = 'IL-8'
geneList = []
organisms = []
numberOrganisms = 0
GCValue = [0, 0, 0, 0, 0]

# create csv-file A4_Schachl_IL-8.csv
with open('A4_Schachl_IL-8.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Accession number", "Title", "Organism", "Sequence length", "GC %"])

    # identify IDs for the first 70
    handle = Entrez.esearch(db="nucleotide", term=geneName, retmax="70")
    record = Entrez.read(handle)
    handle.close()

    # fetch entries
    with open('IL-8.fasta', 'w') as fasta:
        for ID in record["IdList"]:
            # fetch and analyse entries
            handle2 = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", retmode="text")
            record_entry = SeqIO.read(handle2, "genbank")
            geneList.append(record_entry)
            handle2.close()

            # save the first 5 entries of different species
            if (record_entry.annotations["organism"] not in organisms) and (numberOrganisms < 5):
                # save nucleotide sequence
                organisms.append(record_entry.annotations["organism"])
                SeqIO.write(record_entry, fasta, "fasta")
                # calculate GC value
                GCValue[numberOrganisms] = str(GC(record_entry.seq))

                print(numberOrganisms, ': ', organisms[numberOrganisms])

                # write information to csv file
                tableData = [record_entry.annotations['accessions'], record_entry.description,
                             organisms[numberOrganisms], len(record_entry.seq), GCValue[numberOrganisms]]
                writer.writerow(tableData)

                numberOrganisms += 1

# get sequences and find out longest sequence
inputFasta = 'IL-8.fasta'
records = SeqIO.parse(inputFasta, 'fasta')
records = list(records)
maxLen = max(len(record.seq) for record in records)

# pad sequences to equal lengths (based on the longest sequence)
for record in records:
    if len(record.seq) != maxLen:
        sequence = str(record.seq).ljust(maxLen, '.')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxLen for record in records)

# write into temporary output file
outputFasta = '{}_padded.fasta'.format(os.path.splitext(inputFasta)[0])
with open(outputFasta, 'w') as f:
    SeqIO.write(records, f, 'fasta')

# perform alignment
alignment = AlignIO.read(outputFasta, "fasta")
print(alignment)

# phylogenetic tree
constructor = DistanceTreeConstructor(DistanceCalculator('identity'), 'nj')
tree = constructor.build_tree(alignment)
fig = Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length, 3))
