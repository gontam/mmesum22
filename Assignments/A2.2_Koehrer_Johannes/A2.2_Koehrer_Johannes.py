"""
Gathering info from GenBank
Make sure Biopython is available on your machine.
For each gene from the 'Assignment 1 - Research and usage of biological DBs online' gather 5 entries in the GenBank
programmatically (choose most relevant entries and limit your search).
The script must create a table ('A2.2_Chosentopic_Lastname_Name.csv') with the following header:
Accession number, Title, Organism, Sequence length, CG %, Protein instability, Aromaticity, Isoelectric point
Each row must contain the respective data for each of the 5 entries per gene.
Tips:
    Visit the 'Biopython Tutorial and Cookbook'.
    Gene names can be hardcoded.
    Instability, aromaticity, isoelectric point can be calculated after translation.
"""


from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv


def create_file_header():
    # Create a file with the headers
    fields = ['Accession number', 'Title', 'Organism', 'Sequence length', 'CG %',
              'Protein instability', 'Aromaticity', 'Isoelectric point']
    with open("A2.2_SickleCellDisease_Koehrer_Johannes.csv", 'w', newline='') as csv_file:
        csvwriter = csv.writer(csv_file)
        csvwriter.writerow(fields)


def genbank_search():
    # Search Genbank for entries related to the gene list
    gene_list = ['HBB', 'SELP', 'HBA1']
    Entrez.email = "johannes.koehrer@gmail.com"
    collected_data = []

    for gene in gene_list:
        # Search in GenBank for defined term with gene name, limited gene sequence, max 5 entries
        handle = Entrez.esearch(db="nucleotide", term=f"{gene} [gene name] AND 100[SLEN] : 10000[SLEN]",
                                idtype="acc", retmax="5")
        record = Entrez.read(handle)

        for field in record["IdList"]:
            # Save the information from each genbank entry in record
            handle = Entrez.efetch(db="nucleotide", id=field, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            seq = record.seq

            # Translate the DNA sequence for longest protein sequence
            longest_seq = translation(seq)

            # Store obtained data in list
            collected_data.append(store_data(record, longest_seq))

        # Save collected data in file
        save_data(collected_data)
        collected_data.clear()
    print(f"Finished\nSaved whole data in A2.2_SickleCellDisease_Koehrer_Johannes.csv")


def translation(sequence):
    # Translate the DNA sequence for longest protein sequence
    sequence_protein = ""
    longest_sequence = ""
    start_codon = False

    # Check each ORF
    for frame in range(0, 3):
        dna_seq = sequence[frame:]
        # Cut the sequence till only triplets remain
        while len(dna_seq) % 3 != 0:
            dna_seq = dna_seq[:-1]
        # Translate the DNA sequence
        tra = dna_seq.translate()
        for index in range(0, len(tra)):
            # Start when start codon appears or continue
            if tra[index] == 'M' or start_codon is True:
                start_codon = True
                sequence_protein += tra[index]
                # Stop when stop codon appears
                if tra[index] == '*':
                    start_codon = False
                    # Store longest protein sequence
                    if len(longest_sequence) < len(sequence_protein):
                        longest_sequence = sequence_protein
                    sequence_protein = ""
    longest_sequence = longest_sequence[:-1]
    return longest_sequence


def store_data(record, longest_sequence):
    # Store obtained data in a list
    analysed_seq = ProteinAnalysis(longest_sequence)
    rows = [record.id, record.description, record.annotations['organism'], len(record.seq),
            '%.3f' % GC(record.seq), '%.3f' % analysed_seq.instability_index(),
            '%.3f' % analysed_seq.aromaticity(), '%.3f' % analysed_seq.isoelectric_point()]
    return rows


def save_data(data):
    # Save obtained data in the file
    with open("A2.2_SickleCellDisease_Koehrer_Johannes.csv", 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(data)
        print("Save data...")


''' MAIN PROGRAM '''

# Create a file with the headers
create_file_header()

# Search Genbank for entries related to the gene list
genbank_search()
