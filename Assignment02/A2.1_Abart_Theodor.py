def read_fasta(filename):
    f = open(filename, 'r')
    data = f.readlines()
    data = data[1::]  # remove first line
    data = ' '.join([str(line) for line in data])
    return data


def read_standard(filename):
    f = open(filename, 'r')
    data = f.readlines()
    for i in range(data.__len__()):
        data[i] = data[i].split()
    return data


filename_sequence = 'BD137219.1.fasta'
data = read_fasta(filename_sequence)
# print(data)

filename_standard = 'standard.txt'
standard = read_standard(filename_standard)
print(standard)


