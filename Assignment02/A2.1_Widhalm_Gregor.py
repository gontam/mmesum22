import os


def readFASTAFile(filename):

    # open and read fasta file
    file1 = open(filename, 'r')
    lines = file1.readlines()

    # remove first line (description)
    lines.pop(0)

    stringFromList = ' '.join([str(elem) for elem in lines])
    stringFromList = stringFromList.replace('\n ', '')

    finalList = list()
    finalList[:0] = stringFromList
    return finalList


def translateToAminoAcid(toTranslateSeparated, oneOrThreeLetterFlag, initialFilename):

    aminoAcidTable = readAminoAcidTable()

    # adapt OneOrThreeLetterFlag to correct column index of amino acid list
    if oneOrThreeLetterFlag == 1:
        oneOrThreeLetterFlag = 1
    elif oneOrThreeLetterFlag == 3:
        oneOrThreeLetterFlag = 2

    translationIndices = []
    translatedResult = ''
    for i in range(len(toTranslateSeparated)):
        # devide into triplets and replace spaces
        nextTriplet = ''.join([str(elem) for elem in toTranslateSeparated[i]]).replace(' ', '')

        for j in range(len(aminoAcidTable)):
            if aminoAcidTable[j][0] == nextTriplet:
                translatedResult = (''.join((translatedResult, aminoAcidTable[j][oneOrThreeLetterFlag]))) + ' '
                translationIndices.append(j)
                break

    # write result to file
    f = open((initialFilename + '_translated_Widhalm.txt'), "w")
    f.write(translatedResult)
    f.close()

    print('Translated result:\n%s' % ('\n'.join(translatedResult[i:i+16] for i in range(0, len(translatedResult), 16))))


def readAminoAcidTable():

    # open and read translation reference file
    tripletAminoAcidStandard = open('standard.txt', 'r')
    tripletAminoAcidStandardLines = tripletAminoAcidStandard.readlines()

    for i in range(len(tripletAminoAcidStandardLines)):
        tripletAminoAcidStandardLines[i] = tripletAminoAcidStandardLines[i].replace('\n', '')
        tripletAminoAcidStandardLines[i] = list(tripletAminoAcidStandardLines[i].split(' '))

    return tripletAminoAcidStandardLines


def findLongestORF(StartCodonsIndices, StopCodonsIndices):
    MaxORF_Overall = 0
    SequenceIndexToTranslate = 0
    StartIndexToTranslate = 0
    StopIndexToTranslate = 0

    lengthsORF = list()
    StartStopPair = list()
    for i in range(0, 3):

        lengthsORF.append(list())
        StartStopPair.append(list())
        CounterStopAlreadyUsed = 0
        # run through indices but only for the shorter list;
        # --> if 10 start codons but only one stop condon --> taking first start codon is sufficient and vice versa
        for j in range(min([len(StartCodonsIndices[i]), len(StopCodonsIndices[i])])):
            # run trough stop codons;
            # if the index of the found stop codon is smaller than the index of the start codon --> skip the stop codon
            for k in range(0, len(StopCodonsIndices[i])):
                if StopCodonsIndices[i][k + CounterStopAlreadyUsed] > StartCodonsIndices[i][j]:
                    lengthsORF[i].append(StopCodonsIndices[i][k + CounterStopAlreadyUsed] - StartCodonsIndices[i][j])
                    StartStopPair[i].append([StartCodonsIndices[i][j], StopCodonsIndices[i][k + CounterStopAlreadyUsed]])
                    CounterStopAlreadyUsed = CounterStopAlreadyUsed + k + 1  # counter to avoid double checking of small stop codon indices
                    break

        MaxORF_length = max(lengthsORF[i])
        MaxORF_start = StartStopPair[i][lengthsORF[i].index(max(lengthsORF[i]))][0]
        MaxORF_stop = StartStopPair[i][lengthsORF[i].index(max(lengthsORF[i]))][1]
        print('maximum ORF length for +%d: %d @ start: %d end: %d' % (i, MaxORF_length, MaxORF_start, MaxORF_stop))
        # check if current max length of ORFs is longer than overall longest ORF
        if MaxORF_length > MaxORF_Overall:
            MaxORF_Overall = MaxORF_length
            SequenceIndexToTranslate = i
            StartIndexToTranslate = MaxORF_start
            StopIndexToTranslate = MaxORF_stop

    return MaxORF_Overall, SequenceIndexToTranslate, StartIndexToTranslate, StopIndexToTranslate


if __name__ == '__main__':
    # ask for 1 or 3 as input until one of those two is entered
    while True:
        print(f"Hello {os.getlogin()}! \nPlease enter whether you prefer one (1) or three (3) letter coding: ")
        OneOrThreeLetterFlag = input()
        try:
            OneOrThreeLetterFlag = int(OneOrThreeLetterFlag)
            if OneOrThreeLetterFlag == 1:
                break
            elif OneOrThreeLetterFlag == 3:
                break
        except ValueError:
            print("Try again ...")

    print("Perfect, let's go with %s letter coding! \n" % ("one" if OneOrThreeLetterFlag == 1 else "three"))
    FastaFileName = 'BD137219.1'
    # call function to read fasta file
    ToTranslate = readFASTAFile((FastaFileName + '.fasta'))

    # definition of start and stop codons as well as lists for start and stop codon indices
    StartCodon = ['A', 'T', 'G']
    StopCodons = [['T', 'A', 'G'], ['T', 'A', 'A'], ['T', 'G', 'A']]
    StartCodonsIndices = [[]]
    StopCodonsIndices = [[]]
    ToTranslateSeparated = [[]]

    # for loop to shift 0, +1 and +2
    for x in range(0, 3):
        StartCodonsCounter = 0
        StopCodonsCounter = 0

        # add list of triplets
        if x == 0:
            ToTranslateSeparated[0] = [ToTranslate[i:i + 3] for i in range(x, len(ToTranslate), 3)]
        else:
            ToTranslateSeparated.append([ToTranslate[i:i + 3] for i in range(x, len(ToTranslate), 3)])

        for i in range(len(ToTranslateSeparated[x])):
            # identify all codon indices
            if ToTranslateSeparated[x][i] == StartCodon:
                if StartCodonsCounter == 0:
                    if x > 0:
                        StartCodonsIndices.append(list())
                    StartCodonsIndices[x] = [i]
                else:
                    StartCodonsIndices[x].append(i)
                StartCodonsCounter += 1

            # identify all stop codon indices
            for j in range(len(StopCodons)):
                if ToTranslateSeparated[x][i] == StopCodons[j]:
                    if StopCodonsCounter == 0:
                        if x > 0:
                            StopCodonsIndices.append(list())
                        StopCodonsIndices[x] = [i]
                    else:
                        StopCodonsIndices[x].append(i)
                    StopCodonsCounter += 1

    # function call to gain information on longest ORF
    MaxORF_Overall, SequenceIndexToTranslate, StartIndexToTranslate, StopIndexToTranslate = findLongestORF(StartCodonsIndices, StopCodonsIndices)

    # translate sequence of longest ORF
    translateToAminoAcid(ToTranslateSeparated[SequenceIndexToTranslate][StartIndexToTranslate:StopIndexToTranslate], OneOrThreeLetterFlag, FastaFileName)

    print("done.")
