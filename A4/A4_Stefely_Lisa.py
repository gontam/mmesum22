# Assaignment 4: Lisa Stefely

# chosen Gene: GLP1R

#Imports
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.SeqUtils import GC
import csv

#Mail for Entrez
Entrez.email = "me21m007@technikum-wien.at"

gc_value = [0,0,0,0,0]
sequences = []
organisms = []
counter = 0

with open('A4_GLP1R.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Accession number", "Title", "Organism", "Sequence length", "GC %"])

    search_result = Entrez.esearch(db="nucleotide", term="GLP1R[Gene]", retmax="10")
    record = Entrez.read(search_result)

    with open('glp1r.fasta', 'w') as fasta:
        for id in record["IdList"]:

            handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
            glp1r_entry = SeqIO.read(handle, "genbank")
            sequences.append(glp1r_entry)
            handle.close()

            #Stop after 5 different organisms
            if counter>=5:
                break

            # Deleting if already saved organism
            elif (glp1r_entry.annotations["organism"] in organisms):
                continue

            #Saving if new organism
            else:
                organisms.append(glp1r_entry.annotations["organism"])
                # Nucleotide Sequences
                SeqIO.write(glp1r_entry, fasta, "fasta")
                # GC value
                gc_value[counter] = str(GC(glp1r_entry.seq))

            # Write in csv
            writer.writerow([glp1r_entry.id, glp1r_entry.description, glp1r_entry.annotations["organism"], len(glp1r_entry.seq), gc_value[counter]])
            counter += 1

input = './glp1r.fasta'
records = SeqIO.parse(input, 'fasta')
records = list(records)
minlen = min(len(record.seq) for record in records)

# cut in equal lengths (based on smallest sequence)
for rec in records:
        rec.seq = rec.seq[0:minlen]

#write into file
with open(input, 'w') as x:
    SeqIO.write(records, x, 'fasta')

align = AlignIO.read(input, "fasta")
construct = DistanceTreeConstructor(DistanceCalculator('identity'), 'nj')
tree = construct.build_tree(align)
fig = Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length,3))