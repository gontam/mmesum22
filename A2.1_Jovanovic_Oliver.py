# Library for table
import pandas as pd
import numpy as np

# Line Vector for storing data in it.
lines = []
dt = []
# Reading in the whole data and store it in the vector "lines".
with open('G:\\My Drive\\Universität\\Master\\Master - 2. Semester (SS22) - '
          'Auslandssemester\\Bioinformatics\\BD137219.1.fasta', 'r') as f:
    lines = f.readlines()
f.close()

# Read in Dictionary to translate sequence to protein.
# with open('G:\\My Drive\\Universität\\Master\\Master - 2. Semester (SS22) - '
# 'Auslandssemester\\Bioinformatics\\standard.txt', 'r') as g:
# dt = g.readlines()
# g.close()
# Getting information how long the fasta file is.
lengthFile = len(lines)
# df = pd.DataFrame(dt)
# print(df)
dt = [['TTT', 'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TAA', 'TGA', 'TTG', 'TCG', 'TAG', 'TGG',
       'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG',
       'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ATG', 'ACG', 'AAG', 'AGG',
       'GTT', 'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG'],
      ['F', 'S', 'Y', 'C', 'F', 'S', 'Y', 'C', 'L', 'S', '*', '*', 'L', 'S', '*', 'W', 'L', 'P', 'H', 'R', 'L', 'P',
       'H', 'R', 'L', 'P', 'Q', 'R', 'L', 'P', 'Q', 'R', 'I', 'T', 'N', 'S', 'I', 'T', 'N', 'S', 'I', 'T', 'K', 'R',
       'M', 'T', 'K', 'R', 'V', 'A', 'D', 'G', 'V', 'A', 'D', 'G', 'V', 'A', 'E', 'G', 'V', 'A', 'E', 'G'],
      ['Phe', 'Ser', 'Tyr', 'Cys', 'Phe', 'Ser', 'Tyr', 'Cys', 'Leu', 'Ser', 'Ter', 'Ter', 'Leu', 'Ser', 'Ter', 'Trp',
       'Leu', 'Pro', 'His', 'Arg', 'Leu', 'Pro', 'His', 'Arg', 'Leu', 'Pro', 'Gln', 'Arg', 'Leu', 'Pro', 'Gln', 'Arg',
       'Ile', 'Thr', 'Asn', 'Ser', 'Ile', 'Thr', 'Asn', 'Ser', 'Ile', 'Thr', 'Lys', 'Arg', 'Met', 'Thr', 'Lys', 'Arg',
       'Val', 'Ala', 'Asp', 'Gly', 'Val', 'Ala', 'Asp', 'Gly', 'Val', 'Ala', 'Glu', 'Gly', 'Val', 'Ala', 'GLu', 'Gly']]
# Not the first line, because there is information about the dataset in it, which is not important to us.
x = 1
str1 = ""
listNew = []
while x <= lengthFile - 1:
    string = lines[x]
    string = string.strip('\n')
    listNew += string
    print('String', x, 'processed.')
    x += 1

# Convert List in String.
str1 = str1.join(listNew)

# Splitting String to analyze genome.
split_strings = []
n = 3
for index in range(0, len(str1), n):
    split_strings.append(str1[index: index + n])
print(split_strings[0])

print('The whole genome is:', str1)
print('The whole list is:', dt)
# Delete unnecessary variables to save space.
del x, string, listNew, lengthFile, f, lines, n, index

# Implementing Logic
v = int(input("W: "))
a = []
j = 0
i = 0
if v == 1:
    while i <= len(split_strings):
        if len(split_strings[i]) == 3:
            if split_strings[i] == dt[0][j]:
                print(split_strings[i], '=', dt[0][j])
                a = dt[1][j]
                with open('encryption.txt', 'w') as f:
                    f.write(dt[1][j])
                i += 1
                j = 0
            elif j == 64:
                break
                f.close()
            else:
                j += 1
        else:
            break
            f.close()

if v == 2:
    while i <= len(split_strings):
        if len(split_strings[i]) == 3:
            if split_strings[i] == dt[0][j]:
                print(split_strings[i], '=', dt[0][j])
                with open('encryption.txt', 'w') as f:
                    f.write(dt[1][j])
                i += 1
                j = 0
            elif j == 64:
                break
                f.close()
            else:
                j += 1
        else:
            break
            f.close()
