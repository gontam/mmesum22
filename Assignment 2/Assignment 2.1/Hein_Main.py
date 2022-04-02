
# Searching sequence three letter
def markStopCodons(seq):
    edited_seq = ""
    i = 0
    while i < len(seq):
        codon = str(seq[i:i+3])
        # Replace stop codons
        if codon in ["TAG", "TGA", "TAA"]:
            edited_seq += "*"
        else:
            edited_seq += codon
        # Moving in tripplets
        i += 3
    return edited_seq

# Longest ORF search within the codons
def splitFromStartCodon(longest_orf):
    longest_orf_edited = ""
    i = 0
    is_start = False
    while i < len(longest_orf):
        codon = longest_orf[i:i+3]
        if codon == "ATG":
            is_start = True
        if is_start:
            longest_orf_edited += codon
        i += 3
    return longest_orf_edited

# Output of the longest ORF within the whole sequence
def findLongestOrf(seq, threshold):
    orfs = []
    for strand, nseq in [("+", seq)]:
        for frame in range(3):
            length = 3 * ((len(nseq)-frame) // 3)
            for orf in markStopCodons(nseq[frame:frame+length]).split("*"):
                if len(orf) >= threshold:
                    orfs.append([orf, strand])
    max_orf = max([x[0]+x[1] for x in orfs], key=len)
    longest_orf, strand = max_orf[:-1], max_orf[-1]
    longest_orf_edited = splitFromStartCodon(longest_orf)
    return longest_orf_edited


def main():

    # Protein output variables
    protein_sequence_single = ""
    protein_sequence_tripple = ""

    # User selection if a heading line is included
    data_set = int(input("Source file with heading? (Yes = 1) \n"))
    if data_set == 1:
        file = open('BD137219.1.fasta', 'r')
        lines = file.readlines()
        file.close()

        del lines[0]
        blank_file = open('BD137219.1.fasta', 'w+')

        for line in lines:
            blank_file.write(line)
        blank_file.close()

        print ("First line deleted")

        file = open('BD137219.1.fasta', 'r')
        dna = file.read()

    # No heading - read the data into the variable dna
    else:
        file = open('BD137219.1.fasta', 'r')
        dna = file.read()

    # Deleting spaces
    dna = dna.replace('\n', "")
    dna = dna.replace('\r', "")

    # Codon dictionary
    single = {"TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
               "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
               "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
               "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
               "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
               "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
               "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
               "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
               "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
               "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
               "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
               "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
               "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
               "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
               "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
               "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
    protein = {"F": "Phe", "S": "Ser", "Y": "Tyr", "C": "Cys", "F": "Phe",
              "S": "Ser", "Y": "Tyr", "C": "Cys", "L": "Leu", "S": "Ser",
              "*": "Ter", "*": "Ter", "L": "Leu", "S": "Ser", "*": "Ter",
              "W": "Trp", "L": "Leu", "P": "Pro", "H": "His", "R": "Arg",
              "L": "Leu", "P": "Pro", "H": "His", "R": "Arg", "L": "Leu",
              "P": "Pro", "Q": "Gln", "R": "Arg", "L": "Leu", "P": "Pro",
              "Q": "Gln", "R": "Arg", "I": "Ile", "T": "Thr", "N": "Asn",
              "S": "Ser", "I": "Ile", "T": "Thr", "N": "Asn", "S": "Ser",
              "I": "Ile", "T": "Thr", "K": "Lys", "R": "Arg", "M": "Met",
              "T": "Thr", "K": "Lys", "R": "Arg", "V": "Val", "A": "Ala",
              "D": "Asp", "G": "Gly", "V": "Val", "A": "Ala", "D": "Asp",
              "G": "Gly", "V": "Val", "A": "Ala", "E": "Glu", "G": "Gly",
              "V": "Val", "A": "Ala", "E": "Glu", "G": "Gly",
              }

    # Converting DNA into single letter Protein
    for i in range(0, len(dna) - (3 + len(dna) % 3), 3):
        protein_sequence_single += single[dna[i:i + 3]]

    # Converting DNA into tripple letter Protein
    for i in range(0, len(protein_sequence_single) - (1 + len(protein_sequence_single) % 1), 1):
        protein_sequence_tripple += protein[protein_sequence_single[i:i + 1]]

    # Finding longest ORF
    longestorf = str(findLongestOrf(dna, 180))


    # User selection if the output file should include single or tripple letters
    answer = int(input("Result One-Letter (1)? Every other input = Three-Letter \n"))
    if answer == 1:
        with open('BD137219.1_translated_Hein.txt', 'w') as f:
            f.write('Protein Single: \n' + protein_sequence_single + '\n' + 'Longest ORF: \n' + longestorf)
            f.close()


    else:
        with open('BD137219.1_translated_Hein.txt', 'w') as f:
            f.write('Protein Tripple: \n' + protein_sequence_tripple + '\n' + 'Longest ORF: \n' + longestorf)
            f.close()

main()