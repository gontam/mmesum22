# MME
# Assignment 2 - Bioinformatics (BI)
# Author: Johanna Schachl, me21m009
# Date: 28.03.2022
# Last Change: 02.04.2022

# 1. Translation - from DNA to proteins
# Program translates mRNA sequence provided in 'BD137219.1.fasta' and translates it with the help of 'standard.txt'
# User can choose between one- or three-letter coding.
# Result: coding sequences within the start and stop codons is written to ‘BD137219.1_translated_Schachl.txt’

# load information from standard.txt
import Bio.GenBank

with open('standard.txt') as table_triplets:
    triplets = table_triplets.readlines()

table_triplets.close()

# load information from BD137219.1.fasta and remove first line
gene_code = open('BD137219.1.fasta', 'r')
gene = gene_code.read()
gene_code.close()

# remove description and \n from gene
end_of_description = gene.find("\n")
gene_name = gene[0:end_of_description]
gene = gene[end_of_description+1:]
gene = gene.replace('\n', '')

# variables
i = 0   # counter 1
j = 0   # counter 2
k = 0   # counter 3
m = 0   # counter 4
decoded = ['', '', '']   # decoded gene
ORF = 0  # is open reading frame
counter = 0  # counter for ORFs

# choose between one-letter and three-letter coding
test_str = int(input("Enter '1' for one-letter coding and '2' for three-letter coding:\n"))

# compare strings, find ORFs, and encode data
for i in range(3):
    for j in range(0, len(gene) - 4, 3):
        for k in range(len(triplets)):
            if gene[i + j] == triplets[k][0] and gene[i + j + 1] == triplets[k][1] and gene[i + j + 2] == triplets[k][2]:
                if test_str == 1:
                    if ORF == 0 and triplets[k][4] == 'M':
                        counter += 1
                        ORF = 1
                        decoded[i] += '\nORF ' + str(counter) + ': '
                    if ORF == 1 and triplets[k][4] == '*':
                        ORF = 0
                    decoded[i] += triplets[k][4]
                elif test_str == 2:
                    if ORF == 0 and triplets[k][6:9] == 'Met':
                        counter += 1
                        ORF = 1
                        decoded[i] += '\nORF ' + str(counter) + ': '
                    if ORF == 1 and triplets[k][6:9] == 'Ter':
                        ORF = 0
                    decoded[i] += triplets[k][6:9]

# write results to ‘BD137219.1_translated_Lastname.txt’
result = open('BD137219.1_translated_Schachl.txt', 'w')
result.write('Decoded Aminoacids for ' + gene_name + '\n')

for m in range(3):
    result.write('\nTranslation for 5´3´ Frame ' + str(1+m) + ':\n' + decoded[m] + '\n')

result.close()
