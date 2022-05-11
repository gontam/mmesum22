#Assignment 4 - Biopython/Phylogenetics
#HEIN Raphael (me21m022)
#Github:


from Bio import SeqIO
from Bio import Entrez
from Bio import Seq
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqUtils import GC #Calculation of the GC Value
import csv
import os

Entrez.email = "me21m022@technikum-wien.at"


g = ["XM_034936850.1","NM_001134016.1","XM_031653524.1","XR_005576740.1","XM_030935770.1"] # Accession numbers

i = 0 # counter var.

# ----------------------------------- Getting Sequence to same length --------------------------------------
# ----------------------------------- https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length --------------------------------------

input_file = "Sequence.txt"
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
output_file = '{}_aligned.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')
alignment = AlignIO.read(output_file, "fasta")

# ----------------------------------- Creating .csv output file --------------------------------------

outputfile = open('A4_Hein.csv', 'w')
writer = csv.writer(outputfile)
writer.writerow(["Accession number", "Title", "Organism", "Sequence length ", "CG % "])

#Loop till total amount of accession numbers have been reached
while i < len(g):

    # Access to GeneBank searching for nucleotide
    handle = Entrez.efetch(db="nucleotide", id=g[i], rettype="gb")

    gene = SeqIO.read(handle, "genbank")                                    #Retreiving DNA
    handle.close()

    protein_sequence = gene.seq.translate()                                 #DNA to Protein
    seq = str(protein_sequence)
    seq=seq.replace('*','')

    #Output of all parameters (limiting calculation to 4 digit)
    writer.writerow([gene.id,gene.description, gene.annotations["organism"], len(gene.seq), str(GC(gene.seq))[0:5]])

    i += 1

print("Output file created --> Look in directory")



# ----------------------------------- Distance Calculation (NJ - algorithm and Plotting  --------------------------------------

distance = DistanceCalculator('identity')
distance_val = distance.get_distance(alignment)
constructor = DistanceTreeConstructor()
NJTree = constructor.nj(distance_val)
Phylo.draw(NJTree)

