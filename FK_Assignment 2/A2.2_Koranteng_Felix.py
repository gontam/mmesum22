from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData 

import csv

Entrez.email = "me21m001@technikum-wien.at"

#handle = Entrez.einfo(db="nucleotide")
#record = Entrez.read(handle)
#for field in record["DbInfo"]["FieldList"]:
#    print("%(Name)s, %(FullName)s, %(Description)s" % field)
# ACCN, TITL, ORGN, SLEN, 

terms = ["SIRT7[Gene]", "TYR[Gene]", "MITF[Gene]", "PAX3[Gene]"]

with open('A2.2_Melanin_Koranteng_Felix.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Accession number", "Title", "Organism", "Sequence length ", "CG % ", "Protein instability ", "Aromaticity ", "Isoelectric point "])

    for gene in terms:

        db_search = Entrez.esearch(db="nucleotide", term=gene, retmax="5")
        record = Entrez.read(db_search)

        print(record["Count"], " Entries in GenBank for ", gene)
        print("ID of 5 Entries: ")
        print(record["IdList"])

        for gene_id in record["IdList"]:
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
            gene_entry = SeqIO.read(handle, "genbank")
            handle.close()
            X = ProteinAnalysis(str(gene_entry.translate()))

            writer.writerow([gene_entry.id, gene_entry.description, gene_entry.annotations["organism"], len(gene_entry.seq), round(X.get_amino_acids_percent()['G'], 4), 'X.instability_index() -> KeyError: \':\' ' , round(X.aromaticity(), 4) , round(X.isoelectric_point(), 4)])

