#Assignment 2.2 Hein - Used Gene codes:
#KU178287.1
#KJ891296.1
#KR709324.1
#S45136.1
#KJ535072.1

#The user must tell the total amount + which genes
from Bio import SeqIO
from Bio import Entrez #Access to database
from Bio.SeqUtils.ProtParam import ProteinAnalysis #Converting DNA to Protein
from Bio.SeqUtils import GC #Calculation GC Value

import  csv

Entrez.email = "me21m022@technikum-wien.at"

#User selection of total genes
print("How many genes?")

#counter variables
x = input()
i = 0

#Creating new CSV file with heading
outputfile = open('A2.2_GeneticDementia_Hein_Raphael.csv', 'w')
writer = csv.writer(outputfile)
writer.writerow(["Accession number", "Title", "Organism", "Sequence length ", "CG % ",
                 "Protein instability ", "Aromaticity ", "Isoelectric point "])

#Loop till total amount of x has been reached (user selection)
while i < int(x):

    #User selection of genes
    print("ID Please: ")
    id = input()

    # Access to GeneBank searching for nuleotide
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb")

    gene = SeqIO.read(handle, "genbank")                                    #Retreiving DNA
    handle.close()

    #DNA to Protein
    protein_sequence = gene.seq.translate()                                 #DNA to Protein
    seq = str(protein_sequence)
    seq=seq.replace('*','')


    #Output of all parameters (limitng calculation to 4 digit)
    writer.writerow([gene.id,gene.description, gene.annotations["organism"], len(gene.seq), str(GC(gene.seq))[0:5],
                  str(ProteinAnalysis(seq).instability_index())[0:5], str(ProteinAnalysis(seq).aromaticity())[0:5], str(ProteinAnalysis(seq).isoelectric_point())[0:5]])

    i += 1

print("Outputfile created --> Look in directory")