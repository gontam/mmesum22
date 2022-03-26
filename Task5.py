import re
file1 = open('BD137219.fasta', 'r')
Lines = file1.readlines()

#remove first line (description)
Lines.pop(0)

StringFromList = ' '.join([str(elem) for elem in Lines])
StringFromList = StringFromList.replace('\n ', '')
ReversString = StringFromList[::-1]

ReversComplement = []
ReversComplement[:0] = ReversString

for i in range(len(ReversComplement)):
    if ReversComplement[i] == 'A':
        ReversComplement[i] = 'T'
    elif ReversComplement[i] == 'T':
        ReversComplement[i] = 'A'
    elif ReversComplement[i] == 'C':
        ReversComplement[i] = 'G'
    elif ReversComplement[i] == 'G':
        ReversComplement[i] = 'C'

print("REVERSE COMPLEMENT: " + ''.join([str(elem) for elem in ReversComplement]) + '\n')
print("ORIGINAL: " + ''.join([str(elem) for elem in StringFromList])+ '\n')
