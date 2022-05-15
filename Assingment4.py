import Bio
from Bio.Seq import Seq
import re
import pandas as pd
from Bio.SeqUtils import GC
# Import Libraries:
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
import os
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

file="C:\\Users\\usuario\\Desktop\\FH-TECHNIKUM WIEN\\mmesum22\\data.fna"
Entrez.email = 'be21x016@technikum-wien.at'

handle = Entrez.esearch(db='nucleotide', term='APOL1', retmax=5, idtype='acc')
record = Entrez.read(handle)
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta", retmode="text")
out_handle=open(file,"w")
out_handle.write(handle.read())
out_handle.close()
handle.close()
i=0
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
    title = parts[5] +' '+parts[6]+parts[7]+parts[8]
    GC_percentage = GC(sequence)
    rows = [accession_number, title, organism, sequence_length, GC_percentage]
    matrix.append(rows)
data=pd.DataFrame(matrix)
data.to_excel('results.xlsx', sheet_name='sheet1', index=False)

records = SeqIO.parse(file, 'fasta')
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
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
print(alignment)


calculator = DistanceCalculator('identity')
distMatrix = calculator.get_distance(alignment)
#print(distMatrix)

#construct the pyhlogenetic tree using NJ algorithm
constructor = DistanceTreeConstructor()
NJTree = constructor.nj(distMatrix)

#draw the phylogenetic tree
Phylo.draw(NJTree)



