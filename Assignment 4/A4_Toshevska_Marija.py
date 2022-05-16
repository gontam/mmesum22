import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import Entrez
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import os
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

file="gene.fna"
Entrez.email = 'be21x016@technikum-wien.at'

handle = Entrez.esearch(db='nucleotide', term='CFTR', retmax=5, idtype='acc')
record = Entrez.read(handle)
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta", retmode="text")
out_handle=open(file,"w")
out_handle.write(handle.read())
out_handle.close()
handle.close()
i=0

#create table in csv file
columns = ["accession number", "title", "organism", "length", "GC percentage"]
matrix= [columns]
for seq_record in SeqIO.parse(file,"fasta"):
    i=i+1
    sequence_id = seq_record.id
    sequence = repr(seq_record.seq)
    sequence_length = len(seq_record)
    parts=seq_record.description.split()
    accession_number=parts[0]
    organism=parts[1]+' '+parts[2]
    print(parts)
    title = parts[4] +' '+parts[5]
    GC_percentage = GC(sequence)
    rows = [accession_number, title, organism, sequence_length, GC_percentage]
    matrix.append(rows)
data=pd.DataFrame(matrix)
data.to_csv('A4_Marija_Tosh.csv', index=False)
records = SeqIO.parse(file, 'fasta')
records = list(records)
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment
output_file = '{}_padded.fasta'.format(os.path.splitext(file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
alignment = AlignIO.read(output_file, "fasta")
#print(alignment)

calculator = DistanceCalculator('identity')
distMatrix = calculator.get_distance(alignment)
#print(distMatrix)

#create a Pyhlogenetic Tree with NJ algorithm
constructor = DistanceTreeConstructor()
NJTree = constructor.nj(distMatrix)
#draw a Phylogenetic Tree
Phylo.draw(NJTree)






