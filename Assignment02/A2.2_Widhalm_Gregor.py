import Bio
from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import datetime
from datetime import date


def quantityOfEntries(searchTerm, from_date):
    # query current date and format the datetime for genbank query
    now = datetime.date.today()
    date_now_string = now.strftime("%Y/%m/%d")
    date_string = from_date.strftime("%Y/%m/%d")

    # composing query string which should look like the following:
    # (Betacoronavirus) AND ("2019/01/01"[Publication Date] : "3000"[Publication Date])
    # perform query
    handle_query = (Entrez.esearch(db="nucleotide",
                                   term=('(' + searchTerm + ') AND ("' + date_string + '"[Publication Date] : "' + date_now_string +
                                         '"[Publication Date])'), idtype="acc"))
    records = (Entrez.read(handle_query))
    print("""Number of results of query: (%s) AND ("%s"[Publication Date] : "%s"[Publication Date]): \n%d items"""
          % (searchTerm, date_string, date_now_string, int(records["Count"])))


if __name__ == '__main__':

    csv_filename = 'A2.2_Cortisol_Widhalm_Gregor.csv'
    csv_header = ['Accession number', 'Title', 'Organism', 'Sequence length', 'CG %', 'Protein Instability',
                  'Aromaticity', 'Isoelectric point']

    Entrez.email = "me21m012@technikum-wien.at"     # Always tell NCBI who you are

    quantityOfEntries("Betacoronavirus", date(2019, 1, 1))

    # gene names from assignment a1
    GeneList = ['GLCCI1', 'CRHR1', 'IDO1', 'PER1', 'PER2', 'PER3', 'BMAL1', 'CLOCK', 'CRY1', 'GOT2']

    # https: // www.pythontutorial.net / python - basics / python - write - csv - file /
    with open(csv_filename, 'w', encoding='UTF8') as f:
        writer = csv.writer(f)

        writer.writerow(csv_header)

        for gene in GeneList:

            # http://biopython.org/DIST/docs/tutorial/Tutorial.html#chapter%3Aentrez
            handle = Entrez.esearch(db="nucleotide", term=("Homo sapien[Organism] AND " + gene + "[Gene]"))  # retmax="5")
            record = Entrez.read(handle)

            print(gene + ": " + record["Count"] + " items found.")

            counterPerGene = 1
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
                    # only include five entries per gene
                    if counterPerGene > 5:
                        break

                    protein = record_id.seq.translate()

                    # protein analysis: https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html
                    proteinSeqence = str(protein)
                    # had to replace unknown acids coded with 'X' in order to calc analysis parameters
                    proteinSeqence = proteinSeqence.replace('*', '').replace('X', '')
                    # run the protein analysis
                    resProteinAnalysis = ProteinAnalysis(proteinSeqence)

                    writer.writerow([record_id.id, record_id.description, record_id.annotations["organism"],
                                     len(record_id.seq), resProteinAnalysis.get_amino_acids_percent()['G'],
                                     resProteinAnalysis.instability_index(), resProteinAnalysis.aromaticity(),
                                     resProteinAnalysis.isoelectric_point()])

                    print('current GENE added to list ID: %s, No.: %d' % (Id_gene, counterPerGene))
                    counterPerGene += 1

    f.close()
    print("done.")
