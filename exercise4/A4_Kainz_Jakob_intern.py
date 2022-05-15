from Bio import Entrez
from Bio import SeqUtils
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import os
import pandas as pd
from Bio import Seq
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo



# This project was a cooperation of Kainz Jakob and Metzenbauer Gregor.

def fetch_records(email, retmax, term):
    # Identifier for the database
    Entrez.email = email

    # Do to the apparent problem that not all fetched entries can be used, we will fetch more then 5 entries and just
    # stop writting them to the file after we have enough

    handle = Entrez.esearch(db="nucleotide", retmax=retmax, term=term)
    res = Entrez.read(handle)
    idList = res["IdList"]

    # search for all GenBank results for the given IDs
    handle = Entrez.efetch(db='nucleotide', id=idList, rettype='gb', retmode='xml')
    records = Entrez.read(handle)
    unique_organisms = []
    unique_records = []
    for record in records:
        if record["GBSeq_organism"] not in unique_organisms:
            if "GBSeq_sequence" in record.keys() and "GBSeq_primary-accession" in record.keys():
                unique_organisms.append(record["GBSeq_organism"])
                unique_records.append(record)

    # If there are not enough unique entries in the records it needs to fetch again with a larger amount
    if len(unique_organisms) <= 5:
        print('Not enough unique entries found, trying again...')
        unique_records = fetch_records(email, retmax+20, term)

    return unique_records


def write_records(email, term):
    records = fetch_records(email, 20, term)
    print("records collected")
    dirname = os.path.dirname(__file__)
    results_file = os.path.join(dirname, 'Gene_Table.csv')
    gene_sequence_file = os.path.join(dirname, 'Gene_Sequences.csv')

    # definies the parameters to organise the .csv file
    par_accesion = "GBSeq_primary-accession"
    par_title = "GBSeq_locus"
    par_organism = "GBSeq_organism"
    par_length = "GBSeq_length"
    par_sequence = "GBSeq_sequence"

    # Creating the file and it´s header
    result_file = open(results_file, "w")
    result_lines = [
        "Accession number, Title, Organism, Sequence length, CG %"]

    gene_sequence_file = open(gene_sequence_file, "w")
    gene_sequences = []
    written_counter = 0
    for record in records:
        try:
            # Calculating the additional values that are not included in the original record
            gc_percent = SeqUtils.GC(record[par_sequence])

            # Appending all info into a row of the .csv file
            result_lines.append(','.join([
                record[par_accesion],
                record[par_title],
                record[par_organism],
                record[par_length],
                str(gc_percent),
            ]))
            # Generate second file data only containing the necessary information for the tree
            gene_sequences.append(','.join([
                record[par_sequence],
                record[par_organism],
                record[par_accesion],
            ]))
            written_counter += 1
        except:
            continue
        if written_counter == 10:
            break
    result_file.write("\n".join(result_lines))
    gene_sequence_file.write("\n".join(gene_sequences))

def align_sequence():
    # the sequence file contains the gene sequence, its organism name and the accession number
    dirname = os.path.dirname(__file__)
    sequence_file_name = os.path.join(dirname, "./Gene_Sequences.csv")
    sequence = pd.read_csv(sequence_file_name, header=None)
    data = []
    # Creating data containing all the sequences that we need to create the tree

    a = SeqRecord(Seq(sequence[0][0]), id=sequence[1][0]+'_'+sequence[2][0])
    b = SeqRecord(Seq(sequence[0][1]), id=sequence[1][1]+'_'+sequence[2][1])
    c = SeqRecord(Seq(sequence[0][2]), id=sequence[1][2]+'_'+sequence[2][2])
    d = SeqRecord(Seq(sequence[0][3]), id=sequence[1][3]+'_'+sequence[2][3])
    e = SeqRecord(Seq(sequence[0][4]), id=sequence[1][4]+'_'+sequence[2][4])
    records = [a, b, c, d, e]

    # The following lines were for testing purposes which allowed more than 5 entries to be shown on the tree
    if False:
        for index, row in sequence.iterrows():
            data.append(SeqRecord(Seq(row[0]), id=row[1]+'_'+row[2]))
        records = data

    maxlen = max(len(record.seq) for record in records)

    # Fix implemented for ValueError when sequences are not of the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '.')
            record.seq = Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    tmp_file = sequence_file_name
    output_file = '{}_padded.fasta'.format(os.path.splitext(tmp_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')

    # Calculates alignment for the sequences
    align = MultipleSeqAlignment(records)
    print(align)
    # Calculates the distances for the tree
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(align)
    print(dm)

    # Create and draw the tree
    tree_constructor = DistanceTreeConstructor()
    nj_tree = tree_constructor.nj(dm)
    Phylo.draw(nj_tree)

def create_tree(email, gene):
    # Program was tested with a few different genes. 'write_records' method parameter uses search term to create the data
    # records. Check config.py for the genes that we testes
    # The results for the original genes of our team can be found in config.GENE_NAME & config.GENE_NAME_ALTERNATIVE_1.
    write_records(email, gene)
    align_sequence()
