#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


def prepair_coding_seq(codons: pd.core.frame.DataFrame, column: str, seq: str) -> list[str]:
    '''
    INPUT:
        1) pandas DataFrame as table with codons and respective aminoacides
        2) column name string (provided after user answers a question) to search in the DataFrame
        3) sequence string from *.fasta file to work with

    PROCEDURE:
        The function searches for the 'ATG'-start codon. After a start codon has been found all triplets after it are checked and respective aminoacides are found in the DataFrame.
        The names of found aminoacides are added to the string.
        If the last checked triplet was one of stop codons, then the current string is appended to a list, internal for-loop breaks and external for-loop searces for the next start codon.
        If current triplet can not be found in the DataFrame, then such protein can not be accomplished and the internal for-loop breaks without appending current string to the list.
        If the last codon checked is not one of the stop codons, then such protein can not be accomplished and the string is not appended to the list.

    OUTPUT:
        A list of strings with the names of aminoacides. Each string is started with Met or M (for ATG) and ended with a termination marker.
    '''
    seq_array = []
    
    for i in range(len(seq)-3):
        seq_string: str = ''
        if seq[i:i+3] == 'ATG':
            for k in range(i, len(seq)-3, 3):
                try:
                    seq_string += codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]]
                except Exception:
                    break

                if codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]] == 'Ter' or codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]] == '*':
                    seq_array.append(seq_string)
                    break

    return seq_array


def find_the_longest_string(string_array: list[str]) -> str:
    '''
    INPUT:
        A list of strings that represent possible found sequences of aminoacides.

    PROCEDURE:
        The function takes items of the list one by one and compares with the saved string.
        If the length of an item is bigger then length of a saved string, then the value of item is saved to the string.

    OUTPUT:
        A string that contains the longest item from the list.
    '''
    the_longest_string = ''

    for item in string_array:
        if len(item)>len(the_longest_string):
            the_longest_string = item

    return the_longest_string


def main() -> None:
    # Read a file with codons
    translation_data = pd.read_csv("standard.txt", sep=' ', header=None, names=["Triplet", "AA1", "AA3"])
    
    # Open a *.fasta file and make one line string from it without a header
    with open("BD137219.1.fasta", "r") as fasta_file:
        fasta_file_ready = fasta_file.read()
        sequence_line = ('\n'.join(fasta_file_ready.split('\n')[1:])).replace('\n', '')

    # Ask user about the output format
    user_choise = input("1 or 3 letters for aminoacide? ")
    # Check if the input of user is correct and exit the program if there is a mistake
    if user_choise != "1" and user_choise != "3":
        print("please enter 1 or 3")
        return None
    
    # Find all possible sequences of aminoacides and save them as a list
    coding_sequences = prepair_coding_seq(translation_data, "AA"+user_choise, sequence_line)

    # Find the longest sequence from the list and write it to the *.txt file
    the_longest_coding_seq = open("BD137219.1_translated.txt", 'w')
    the_longest_coding_seq.write(find_the_longest_string(coding_sequences))
    the_longest_coding_seq.close()


if __name__ == '__main__':
    main()

