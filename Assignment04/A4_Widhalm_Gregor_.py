import csv
import os
import warnings
import Bio
import matplotlib
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import Entrez
from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator


if __name__ == '__main__':

    # gene = 'TNF' --> used for testing
    gene = 'PER1'

    # define csv filename and header
    csv_filename = 'A4_'+gene+'_Widhalm_Gregor.csv'
    csv_header = ['Accession number', 'Title', 'Organism', 'Sequence length', 'CG %']

    with open(csv_filename, 'w', encoding='UTF8') as f:
        # open csv file and write header
        writer = csv.writer(f)
        writer.writerow(csv_header)

        Entrez.email = "me21m012@technikum-wien.at"     # Always tell NCBI who you are
        # perform gene query
        handle = Entrez.esearch(db="nucleotide", term=(gene + "[Gene]"))
        record = Entrez.read(handle)

        # print how many records could be found
        print(gene + ": " + record["Count"] + " items found.")

        org = []
        records = []
        gcPerc = []
        # run through the list of records until entries of five different organisms were found
        for Id_gene in record["IdList"]:
            handle_id = Entrez.efetch(db="nucleotide", id=Id_gene, rettype="gb", retmode="text")
            record_id = SeqIO.read(handle_id, "genbank")
            handle_id.close()

            # DUE TO FACT THAT SOME ENTRIES DO NOT INCLUDE SEQUENCE!!!
            try:
                str(record_id.seq)
            except Bio.Seq.UndefinedSequenceError:
                print("WARNING: Sequence was not defined for this entry! Gene ID: %s" % Id_gene)
                print("--> excluded")
            else:
                # 5 Organisms
                if len(org) == 5:
                    break
                # find record of organism not included in the organism list yet
                elif record_id.annotations["organism"] not in org:
                    print('YEAH! - NEW ORGANISM: ', record_id.annotations["organism"])
                    org.append(record_id.annotations["organism"])
                    # add current record to the records list
                    records.append(record_id)

                    # catch warning indicating that sequence length is not a multiple of three:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', category=Bio.BiopythonWarning)
                        ProtRes = ProteinAnalysis(str([(record_id.seq.translate())]))

                    # add GC percentage of current sequence to list
                    gcPerc.append(ProtRes.get_amino_acids_percent()['G'])
                    # write information to csv file
                    writer.writerow([record_id.id, record_id.description, record_id.annotations["organism"],
                                     len(record_id.seq), gcPerc[-1]])

    # save all five sequences to one fasta file
    SeqIO.write(records, gene + '_5Organisms.fasta', 'fasta')

    # sequence padding according to:
    # https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
    maximum_len = max(len(record.seq) for record in records)
    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maximum_len:
            sequence = str(record.seq).ljust(maximum_len, '.')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maximum_len for record in records)

    # write to temporary file and do alignment
    output_file = '{}_padded.fasta'.format(os.path.splitext(gene + '_5Organisms.fasta')[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')

    # Tutorial used: https://taylor-lindsay.github.io/phylogenetics/
    # retrieve alignment of padded fasta file
    alignment = AlignIO.read(output_file, "fasta")
    #print(alignment)
    calculator = DistanceCalculator('identity')
    distanceMatrix = calculator.get_distance(alignment)
    #print(distanceMatrix)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    print(tree)
    # plot phylogenetic tree including distances limited to 3 decimal places
    fig = plt.figure(figsize=(13, 5), dpi=200)
    axes = fig.add_subplot(1,1,1)
    Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length, 3), axes=axes)
    fig.savefig("A4_"+gene+"Widhalm_Gregor_Tree")
    print('DONE.')
