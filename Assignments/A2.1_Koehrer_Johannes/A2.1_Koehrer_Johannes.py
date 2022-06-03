"""
Translation - from DNA to proteins
Important: libraries that perform these operations automatically can not be used.
'standard.txt' contains information regarding which triplets of DNA code for which amino acid.
Using this file,translate the mRNA sequence provided in the 'BD137219.1.fasta'.
The script must provide the user with the possibility to choose whether one-letter ('M') or the three-letter ('Met')
coding will be used.
The result must be written in a file with the title ‘BD137219.1_translated_Lastname.txt’

Tip:
Consider that it is necessary to identify all of the open reading frames (ORFs).
The goal is to find the longest coding sequence within the start and stop codons and provide that one to the user.
Use https://web.expasy.org/translate/ to check your result.
"""


def prompt_user():
    # The user is asked if the translated DNA sequence should be saved as One- or Three-letter code
    user_input = input("Please press\n1 for One-Letter Coding\n3 for Three-Letter Coding\n")
    if user_input == '1':
        print(f"Selected One-Letter Coding")
    elif user_input == '3':
        print(f"Selected Three-Letter Coding")
    return user_input


def get_letter_code(user_input):
    # Based on users input the One- or Three-Letter code will be extracted from the standard file
    letter_code = {}
    with open("standard.txt", 'r') as fileObject:
        file_content = fileObject.readlines()
    # parse through standard.txt and get either the One- or Three-Letter code
    for i in range(0, len(file_content)):
        if user_input == '1':
            letter_code[file_content[i][0:3]] = file_content[i][4]
        else:
            letter_code[file_content[i][0:3]] = file_content[i][6:9]
    return letter_code


def get_sequence_dna():
    # Return the DNA sequence from the file
    with open("BD137219.1.fasta", 'r') as fileObject:
        file_content = fileObject.read()
    # separate header file from the DNA sequence
    sequence_dna = file_content.split('\n', 1)[1]
    # remove all \n within the sequence
    sequence_dna = sequence_dna.replace("\n", "")
    return sequence_dna


def translation(letter_code, sequence_dna):
    # Translates the DNA sequence into the Protein sequence
    sequence_protein = ""
    longest_sequence = ""
    start_codon = False
    # to translate each open reading frame
    for frame in range(0, 3):
        for i in range(frame, len(sequence_dna)-2, 3):
            # shortens the DNA sequence to guarantee triplets
            while len(sequence_dna) % 3 != 0:
                sequence_dna = sequence_dna[:-1]
            base = sequence_dna[i:i + 3]
            # if a start codon appears or already appeared in this run
            if letter_code[base] == 'M' or letter_code[base] == "Met" or start_codon is True:
                start_codon = True
                sequence_protein += letter_code[base]
                # if a stop codon appears, the translation will be stopped
                if letter_code[base] == '*' or letter_code[base] == "Ter":
                    start_codon = False
                    # save the longest translated sequence
                    if len(longest_sequence) < len(sequence_protein):
                        longest_sequence = sequence_protein
                    sequence_protein = ""
        sequence_protein = ""
    if longest_sequence[-1:] == '*':
        return longest_sequence[:-1]
    else:
        return longest_sequence[:-3]


def save_sequence(sequence_protein):
    # Save obtained protein sequence into a new file
    with open("BD137219.1_translated_Koehrer.txt", 'w') as fileObject:
        fileObject.write(sequence_protein)
        print("Protein sequence was saved in 'BD137219.1_translated_Koehrer.txt'")


''' MAIN PROGRAM '''

# Prompt the user if he wants One-letter or Three-letter coding
inputUser = prompt_user()
# Stores either One-letter or Three-letter code
letterCode = get_letter_code(inputUser)
# Store the DNA sequence that will be translated
sequenceDNA = get_sequence_dna()
# Translate the DNA sequence in all ORFs and return the longest protein sequence
longestSequence = translation(letterCode, sequenceDNA)
# Save the translation in a file
save_sequence(longestSequence)
