#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


def prepair_coding_seq(codons: pd.core.frame.DataFrame, column: str, seq: str) -> list[str]:
    seq_array = []
    
    for i in range(len(seq)-3):
        seq_string: str = ''
        if seq[i:i+3] == 'ATG':
            for k in range(i+3, len(seq)-3, 3):
                try:
                    seq_string += codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]]
                except Exception:
                    seq_string = ''

                if codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]] == 'Ter' or codons[column][codons[codons["Triplet"] == seq[k:k+3]].index[0]] == '*':
                    seq_array.append(seq_string)

    return seq_array


def find_the_longest_string(string_array: list[str]) -> str:
    the_longest_string = ''

    for item in string_array:
        if len(item)>len(the_longest_string):
            the_longest_string = item

    return the_longest_string


def main() -> None:
    translation_data = pd.read_csv("standard.txt", sep=' ', header=None, names=["Triplet", "AA1", "AA3"])

    with open("BD137219.1.fasta", "r") as fasta_file:
        fasta_file_ready = fasta_file.read()
        sequence_line = ('\n'.join(fasta_file_ready.split('\n')[1:])).replace('\n', '')

    user_choise = input("1 or 3 letters for aminoacide? ")
    if user_choise != "1" and user_choise != "3":
        print("please enter 1 or 3")
        return None
    
    coding_sequences = prepair_coding_seq(translation_data, "AA"+user_choise, sequence_line)

    the_longest_coding_seq = open("BD137219.1_translated.txt", 'w')
    the_longest_coding_seq.write(find_the_longest_string(coding_sequences))
    the_longest_coding_seq.close()


if __name__ == '__main__':
    main()

