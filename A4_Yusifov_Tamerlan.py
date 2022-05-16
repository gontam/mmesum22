#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Made by:
#   Ukhnalyov Andrey
#   Yusifov Tamerlan

# Improting libraries
import Bio
from Bio import Entrez, Seq, SeqIO, AlignIO, Phylo
from Bio.SeqIO import FastaIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import csv
import os


def search_in_db(gene: str, entries_number: int, data_base: str) -> tuple[list[list[any]], list[any]]:
    array_to_write = []
    array_of_fetch_results = []
    added_organisms = []
    # Search for gene in provided database
    search = Entrez.esearch(term=gene, retmax=entries_number, db=data_base)
    search_results = Entrez.read(search)
    search.close()

    for element in search_results['IdList']:
        data_fetch = Entrez.efetch(id=element, db=data_base, rettype="gb", retmode="text")
        data_fetch_results = SeqIO.read(data_fetch, "genbank")
        data_fetch.close()

        try:
            str(data_fetch_results.seq)
        except Bio.Seq.UndefinedSequenceError:
            print("Sequence error has happened! Gene ID: %s" % element)
        else:
            if (data_fetch_results.annotations["organism"] not in added_organisms) and len(added_organisms) < 5:
                analysis_results = ProteinAnalysis(str(data_fetch_results.translate()))
                row_to_write = [data_fetch_results.id, data_fetch_results.description,
                                data_fetch_results.annotations["organism"], len(data_fetch_results.seq),
                                analysis_results.get_amino_acids_percent()['G']]

                array_to_write.append(row_to_write)
                array_of_fetch_results.append(data_fetch_results)
                added_organisms.append(data_fetch_results.annotations["organism"])

    return array_to_write, array_of_fetch_results


def make_alignment(alignment_array: list[any]) -> Bio.Align.MultipleSeqAlignment:
    # procedure based on example:
    # https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
    max_len = max(len(element.seq) for element in alignment_array)
    # Make the same length for all sequences
    for element in alignment_array:
        if len(element.seq) != max_len:
            sequence = str(element.seq).ljust(max_len, '.')
            element.seq = Seq.Seq(sequence)
    assert all(len(element.seq) == max_len for element in alignment_array)

    # Write to fasta file for alignment
    output_file = '{}_for_alignment.fasta'.format(os.path.splitext('OCA2_organisms_record.fasta')[0])
    with open(output_file, 'w') as f:
        SeqIO.write(alignment_array, f, "fasta")

    return AlignIO.read(output_file, "fasta")


def main() -> None:
    # Cofigure Email
    Entrez.email = "me21m006@technikum-wien.at"

    # Gene name
    gene_for_search = "OCA2"

    # Making a .csv table
    with open("A4_OCA2.csv", "w") as file:
        csv_table = csv.writer(file)
        # Writing a header
        csv_table.writerow(["Accession number", "Title", "Organism", "Sequence length", "CG %"])

        # Making a list of rows tu put into the table
        rows, list_for_alignment = search_in_db(gene_for_search, 50, "nucleotide")
        # print(type(rows))
        # print(type(list_for_alignment))

        # Putting all rows into the table
        for row in rows:
            csv_table.writerow(row)

        # Make alignment
        alignment_results = make_alignment(list_for_alignment)

        # print(type(alignment_results))

        # Build tree
        tree_builder = DistanceTreeConstructor(DistanceCalculator('identity'), 'nj')
        tree = tree_builder.build_tree(alignment_results)
        print(tree)
        # Draw graph based on tree
        graph = Phylo.draw(tree, branch_labels=lambda branch: round(branch.branch_length, 4))


if __name__ == "__main__":
    main()