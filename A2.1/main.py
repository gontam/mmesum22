# Assignment 2.1 by Lisa Stefely

format_type = input("Enter 'M' or 'Met' as format type:")
# Open files
data_nucleotides = open('BD137219.1.fasta', 'r+')
data_standard = open('standard.txt', 'r')
new_file = open('BD137219.1_translated_Stefely.txt', 'w')

# read data
standard_lines = data_standard.readlines()
data_standard.close()
proteins = data_nucleotides.read()
data_nucleotides.close()

# define start and stop codons
start_codon = 'ATG'
stop_codon_1 = 'TAG'
stop_codon_2 = 'TAA'
stop_codon_3 = 'TGA'

# delete header
index = proteins.find('\n')
proteins = proteins[index:].replace('\n', '')

# define variables for checking if protein
counter = 0 # makes sure, that stop codon is also printed in txt
is_protein =0

# loop for checking triplets and writing into new txt file
for index in range(0,len(proteins)-3,3):# steps of three
    save_codon=proteins[index:index+3]# saving one codon
    if (save_codon == start_codon):# checking if start of protein
        is_protein = 1
    elif ((save_codon == stop_codon_1) or (save_codon == stop_codon_2) or (save_codon == stop_codon_3)):# checking if stop of protein
        counter = 1
    if (is_protein == 1):# if it is protein, write into file
        for j in range(len((standard_lines))):# checking for correct codon in standard.txt
            save_std = standard_lines[j][0:3]
            if save_codon == save_std:
                if format_type == 'M':
                    new_file.write(standard_lines[j][4])
                elif format_type == 'Met':
                    new_file.write(standard_lines[j][6:9])  # data_standard.read(3)
    if (counter == 1):# if stopcodon, make sure that next codon only prints starting with the next start codon
        is_protein = 0
        counter = 0
        new_file.write('\n')# for easier reading of one ORF
new_file.close()