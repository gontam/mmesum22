from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import csv


Entrez.email = "johannes.koehrer@gmail.com"
species_list = []
loop_counter = 0

# Create the csv-file and add the headers of the columns
with open("A4_Koehrer_Johannes.csv", 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Accession number, Title, Organism, Sequence, length, CG %"])

# Create the fasta file
with open("sequences_raw.fasta", 'w') as fileObject:
    fileObject.write("")

# Search in Genbank for the gene 'HBB'
handle = Entrez.esearch(db="nucleotide", term="HBB [gene name]", idtype="acc", retmax="15")
record = Entrez.read(handle)

# Collect the available data from the first 15 results
for field in record["IdList"]:
    loop_counter += 1
    print(f"\nLoop  {loop_counter}")
    handle = Entrez.efetch(db="nucleotide", id=field, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Add species only to the list when it appears for the first time
    if record.annotations['organism'] in species_list:
        print(f"Skipped {record.annotations['organism']} - is already in list")
        continue
    else:
        print(f"Added {record.annotations['organism']} to the list")
        species_list.append(record.annotations['organism'])

    # Write relevant data into csv-file
    with open("A4_Koehrer_Johannes.csv", 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)

        csvwriter.writerow([record.id, record.description, record.annotations['organism'],
                            len(record.seq), '%.3f' % GC(record.seq)])

    # Save sequences in fasta-file for alignment
    with open('sequences_raw.fasta', 'a') as fileObject:
        fileObject.write(f">{record.id}\n")
        fileObject.write(str(record.seq) + "\n")

    # Limit the result to 5 different species
    if len(species_list) == 5:
        break


# Parse through the fasta-file and get the length of longest sequence
records = SeqIO.parse('sequences_raw.fasta', 'fasta')
records = list(records)
maxlen = max(len(record.seq) for record in records)

# Pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '.')
        record.seq = Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# Write patted sequence to a new file
with open('sequences_padded.fasta', 'w') as f:
    SeqIO.write(records, f, 'fasta')

# Align the sequences to each other
alignment = AlignIO.read('sequences_padded.fasta', 'fasta')

# Calculate the distance between the sequences
calculator = DistanceCalculator('identity', 'nj')
dm = calculator.get_distance(alignment)

# Generate the Phylogenetic tree visualisation
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
fig = Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length, 3))
