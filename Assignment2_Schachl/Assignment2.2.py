# MME
# Assignment 2 - Bioinformatics (BI)
# Author: Johanna Schachl, me21m009
# Date: 28.03.2022
# Last Change: 03.04.2022

# 2. Gathering info from GenBank
# script creates table 'A2.2_Lactase_Schachl_Name.csv'

# import Libraries
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


import csv

Entrez.email = "me21m009@technikum-wien.at"

# genes
gene = ['LCT [gene]', 'IL-8[gene]', 'IL-10[gene]', 'TNF-alpha[gene]']

# create csv-file and writer
result = open('A2.2_Lactase_Schachl_Johanna.csv', 'w')
writer = csv.writer(result)

header = ['Accession number', 'Title', 'Organism', 'Sequence length', 'CG %', 'Protein instability', 'Aromaticity', 'Isoelectric point']
writer.writerow(header)


# search in GenBank
for i in gene:
    # identify ID for the first 5 entries
    handle = Entrez.esearch(db="nucleotide", retmax=5, term=i)
    record = Entrez.read(handle)
    handle.close()

    # fetch and analyse entries
    for j in record['IdList']:
        handle = Entrez.efetch(db="nucleotide", id=j, rettype="gb", retmode="text")
        record2_iterator = SeqIO.parse(handle, "genbank")
        record2 = next(record2_iterator)
        handle.close()

        gene_analysis = ProteinAnalysis(repr(record2.seq))

        # save data in cvs-file
        data_gene = [record2.annotations['accessions'], record2.description, record2.annotations['organism'], gene_analysis.length, 'gene_analysis.instability_index() -> Error', gene_analysis.aromaticity(), gene_analysis.isoelectric_point()]
        writer.writerow(data_gene)

# close cvs-file
result.close()