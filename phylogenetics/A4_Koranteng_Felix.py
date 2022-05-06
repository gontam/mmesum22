from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import csv
import os

Entrez.email = "me21m001@technikum-wien.at"

gcs = [0 for _ in range(5)]

with open('A4_SIRT7_Koranteng_Felix.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Accession number", "Title", "Organism", "Sequence length ", "GC % "])

    db_search = Entrez.esearch(db="nucleotide", term="SIRT7[Gene]", retmax="60")
    record = Entrez.read(db_search)

    sequences = []
    organisms = []
    tracker = 0
    with open('sirt7.faa', 'w') as fasta:
        for gene_id in record["IdList"]:
            
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
            gene_entry = SeqIO.read(handle, "genbank")
            sequences.append(gene_entry)
            handle.close()

            # 5 Organisms
            if len(organisms) == 5:
                break
            
            # Organism already considered -> Unique
            elif(gene_entry.annotations["organism"] in organisms):
                continue

            # Weird Data content -> clean
            elif(gene_entry.annotations["organism"][:5] == "Blomi" or gene_entry.annotations["organism"][:5] == "Tachy" or gene_entry.annotations["organism"][:5] == "Mugil" or gene_entry.annotations["organism"][:5] == "Hippo" or gene_entry.annotations["organism"][:5] == "Girar" or gene_entry.annotations["organism"][:5] == "Scoph" or gene_entry.annotations["organism"][:5] == "Ictal" or gene_entry.annotations["organism"][:5] == "Marmo" or gene_entry.annotations["organism"][:5] == "Micro" or gene_entry.annotations["organism"][:5] == "Meles" or gene_entry.annotations["organism"][:5] == "Ursus"):
                continue

            # New Organism!
            else:
                organisms.append(gene_entry.annotations["organism"])
                
                # Nucleotide Sequences
                SeqIO.write(gene_entry, fasta, "fasta")

                # GC Percentage 
                trans = gene_entry.seq.translate()
                X = ProteinAnalysis(repr(trans))
                gcs[tracker] = X.get_amino_acids_percent()['G']
                
                tracker += 1

            # Storing Basic Info
            writer.writerow([gene_entry.id, gene_entry.description, gene_entry.annotations["organism"], len(gene_entry.seq)])

input_file = '/home/kobbyayolo/Desktop/Semester 2/Bioinformatics/phylogenetics/sirt7.faa'
records = SeqIO.parse(input_file, 'fasta')
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment
output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')


alignment = AlignIO.read(output_file, "fasta")
print(organisms)
print(gcs)
print(alignment)
constructor = DistanceTreeConstructor(DistanceCalculator('identity'), 'nj')
tree = constructor.build_tree(alignment)
print(tree)
fig = Phylo.draw(tree, branch_labels=lambda d: round(d.branch_length, 4))
