import numpy as np
import os


def translate_dna_to_protein():
    coding_mode = input("Select coding mode, 1 for 1 letter or 2 for 3 letter coding")
    coding_mode = int(coding_mode)
    dirname = os.path.dirname(__file__)
    results_file = os.path.join(dirname, 'A2_2_BRCA1_Kainz_Jakob.csv')
    # define some basic variables needed later on
    dna = ""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    possible_readings = []
    orf_number = 0
    in_reading = False

    # open all needed files
    prot_dict = np.genfromtxt(os.path.join(dirname, 'standard.txt'), dtype='str')
    dna_file = open(os.path.join(dirname, "BD137219.1.fasta"), "r")
    prot_file = open(os.path.join(dirname, "BD137219.1_translated_Kainz.txt"), "w")
    with open(os.path.join(dirname, "BD137219.1.fasta"), "r") as f:
        next(f)
        for line in f:
            dna += line

    # transform dna to list format
    dna = [dna[h:h + 3] for h in range(0, len(dna), 3)]
    # crosscheck list with dictionary
    for i in dna:
        for j in prot_dict:
            if i in j:
                if i == start_codon and not in_reading:
                    # if it finds a start codon it will create itself another reading frame to fill the data
                    in_reading = True
                    possible_readings.append('')
                if in_reading:
                    # appends all the codons to the current reading frame
                    possible_readings[orf_number] += j[coding_mode]
                    # If a stop codon is found, the number of open reading frames increases and it will start to look
                    # for another start codon
                    if i in stop_codons:
                        orf_number += 1
                        in_reading = False

    # Now we only need to find the longest ORF and save it to our file.
    longest_reading = max(possible_readings, key=len)
    prot_file.write(longest_reading)

    # just to add on, we will print all the ORFs we found:
    print("Saving longest ORF, but FYI those were all ORFs that could be found:")
    print(possible_readings)

    # close the remaining open files
    dna_file.close()
    prot_file.close()
