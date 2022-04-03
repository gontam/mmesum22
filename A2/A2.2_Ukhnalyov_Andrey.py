#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Importing libraries
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import datetime
import csv

# List of gens to search for
gens = ["OCA2", "PAX3", "EDNRB", "SOX10", "MITF"]


def ammount_of_enteries_from_date(term: str, start_date: str) -> str:
    '''
    INPUT:
        Term to search and start_date of publications

    PROCEDURE:
        The function finds the current date and then makes a search with a build up string-term.

    OUTPUT:
        A functional string that contains the ammount of results found
    '''
    today = datetime.date.today().strftime("%Y/%m/%d")
    search = Entrez.esearch(db="nucleotide", term='('+term+')'+' AND '+'("'+start_date+'"[Publication Date] : "'+today+'"[Publication Date]'+')')
    search_results = Entrez.read(search)
    answer_string = f'The ammount {term} results form {start_date} to {today} is {search_results["Count"]}'

    return answer_string


def search_in_db(topics: list[str], entries_number: int, data_base: str) -> list[list[any]]:
    '''
    INPUT:
        List of topics to search for, ammount of entries neded and a data_base to search in

    PROCEDURE:
        The function makes a search and takes mentioned number of results for each topic.
        Then it creats a list that includes information of interest for each result of each topic and appends it to a main list
    
    OUTPUT:
        List of lists with the information of interest for all found results
    '''
    array_to_write = []

    for topic in topics:
        search = Entrez.esearch(term=topic, retmax=entries_number, db=data_base)
        search_results = Entrez.read(search)
        search.close()

        for element in search_results['IdList']:
            data_fetch = Entrez.efetch(id=element, db=data_base, rettype="gb", retmode="text")
            data_fetch_results = SeqIO.read(data_fetch, "genbank")
            data_fetch.close()

            analysis_results = ProteinAnalysis(str(data_fetch_results.translate()))
            row_to_write = [data_fetch_results.id, data_fetch_results.description, data_fetch_results.annotations["organism"], len(data_fetch_results.seq),
                              analysis_results.get_amino_acids_percent()['G'], 'analysis_results.instability_index() -> returns KeyError', analysis_results.aromaticity(), analysis_results.isoelectric_point()]
            
            array_to_write.append(row_to_write)

    return array_to_write


def main() -> None:
    # Providing identification
    Entrez.email = "me21m018@technikum-wien.at"
    
    # Printing the ammount of ammount of enteries
    print(ammount_of_enteries_from_date("Betacoronavirus", "2019/01/01"))
    
    # Making a .csv table
    with open("A2.2_Heterochromia_Ukhnalyov_Andrey.csv", "w") as file:
        csv_table = csv.writer(file)
        # Writing a header
        csv_table.writerow(["Accession number", "Title", "Organism", "Sequence length", "CG %", "Protein instability", "Aromaticity", "Isoelectric point"])
        
        # Making a list of rows tu put into the table
        rows = search_in_db(gens, 5, "nucleotide")
        
        # Putting all rows into the table
        for row in rows:
            csv_table.writerow(row)


if __name__ == '__main__':
    main()

