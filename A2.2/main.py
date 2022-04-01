#Assignment 2.2 by Lisa Stefely

import csv

from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#Mail for fetching data
Entrez.email = "me21m007@technikum-wien.at"

#list of selected Nucleotides (5 for GCG, 4 for GLP1R and 5 for GLP2R selected)
list_ids = ["BT006813.1", "AY890069.1", "AY890068.1", "J04040.1", "BC005278.1", "NM_002062.5", "KR138540.1", "BC113493.1", "BC112126.1", "NM_004246.3", "XM_011524077.3", "XM_005256861.3", "XM_017025341.1", "XM_017025340.1"]
new_file = open('A2.2_Glucagon_Stefely_Lisa.csv', 'w')
writer = csv.writer(new_file)

#loop for fetching data and writing it into csv file
for index in range(0,14,1):
    handle = Entrez.efetch(db="nucleotide", id=list_ids[index], rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    protein = record.seq.translate()
    seq = str(protein)
    seq=seq.replace('*','')
    analysed_seq = ProteinAnalysis(seq)
    gc_value = str(GC(record.seq))
    handle.close()
    write_list = ['Acc.-Num.: ' + record.id, ' Title: ' + record.description, ' Organism: ' + record.annotations["source"], ' Seq.-Length: ' + str(len(record.seq)), ' GC-Value: ' + gc_value, ' Instability: '+ str(analysed_seq.instability_index()), ' Aromaticity: ' + str(analysed_seq.aromaticity()), ' Isoelectric Point: ' + str(analysed_seq.isoelectric_point())]
    writer.writerow(write_list)
new_file.close()