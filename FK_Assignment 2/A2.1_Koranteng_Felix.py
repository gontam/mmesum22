from Bio import SeqIO

user_def = int(input("1 -> Single Letter \n2 -> Three-Letter \nSelect Output Option: "))

def get_triplet_reference():
    '''
    Read Standard for Amino Acids from standard.txt
    '''
    code_type = user_def
    if code_type == 1:
        code_type = 4
    else:
        code_type = 6
    triplet_codes = []
    with open('standard.txt') as ref:
        lines = ref.readlines()
        for line in lines:
            if code_type == 4:
                # Single Letter Representation
                triplet_codes.append([line[0:3], line[code_type]])
            else:
                # Triple Letter Representation
                triplet_codes.append([line[0:3], line[code_type:-1]])
    return triplet_codes

def find_rep(code, reference):
    '''
    Get Amino Acid from Codon using Referenced Standard
    '''
    for trial in reference:
        # Check for matching Codon
        if trial[0] == code:
            return trial[1]

def find_longest(frames):
    '''
    Provide Longest identified Frame to User
    '''
    lengths = []
    for frame in frames:
        lengths.append(len(frame))
    for i in range(len(lengths)):
        if lengths[i] == max(lengths):
            idx = i
            print("Longest Identified Frame: ", frames[idx])
            break
    return

def translate():
    '''
    Translate from fasta file -> present all identified reading frames
    '''
    reading_frames = [[],[],[]]
    protein = 'Translated Amino Acids: \n'
    trip_ref = get_triplet_reference()
    #Parse Input from FASTA file
    for record in SeqIO.parse(open('BD137219.1.fasta'), 'fasta'):
        string_DNA = record.seq
        #res = record.seq.translate()
        for j in range(0, 3):
            ORF = False
            protein += 'Off-set =' + str(j) + '\n'
            #Deletion of 0, 1, 2 Nucleotides from Start
            mod_DNA = string_DNA[j:]
            #Loop through in 3s -> codon
            for i in range (int(len(mod_DNA)/3)):
                i *= 3
                codon = mod_DNA[i:i+3]
                amino_acid = find_rep(codon, trip_ref)
                #Check not in Open Reading Frame and Start Codon
                if codon == 'ATG' and ORF == False:
                    protein += '\nNew reading Frame:\n'
                    ORF = True
                    frame = ''
                #Add amino acid to result
                protein += str(amino_acid) + ' '
                if(ORF):
                    frame += amino_acid + ' '
                    #STOP Codons -> Close the Reading frame and 
                    if (codon == 'TAA' or codon == 'TAG' or codon == 'TGA'):
                        protein += '\n'
                        reading_frames[j].append(frame)
                        ORF = False
            protein += '\n\n\n'
    print("Identified Frames: ")
    #print(res)
    for version in reading_frames:
        find_longest(version)

    with open('BD137219.1_translated_Koranteng.txt', 'w') as f:
        f.write(protein)
    return 

translate()
