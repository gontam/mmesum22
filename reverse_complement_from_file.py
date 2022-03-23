#!/usr/bin/env python
# -*- coding: utf-8 -*-

def prepare_sequence(file_name: str) -> str:
    # prepare a string to process by splitting a content of provided file into rows, removing the first row
    # and joining all others to one string    
    opened_file = open(file_name, "r")
    file_to_read = opened_file.read()
    file_to_process = ('\n'.join(file_to_read.split('\n')[1:])).replace("\n", "")
    opened_file.close()

    return file_to_process


def make_reverse_complement(file_to_process: str) -> str:
    # making a reverse complement by reverting the provided string and then creating a new one by appending the value
    # that correspondes to a key from a dictionary
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    file_reverse = file_to_process[::-1]
    reverse_complement = [complement[base] for base in file_reverse]

    return ''.join(reverse_complement)


def main():
    prepared_reverse_complement = make_reverse_complement(prepare_sequence("BD137219.1.fasta"))
    print(prepared_reverse_complement)


if __name__ == '__main__':
    main()

