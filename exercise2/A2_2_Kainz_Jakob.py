from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os


def fetch_records():
    # Identifier for the database
    Entrez.email = "me21m014@technikum-wien.at"

    # gets first 5 results for the term "BRCA1" and saves their ID
    handle = Entrez.esearch(db="nucleotide", retmax=5, term="BRCA1")
    res = Entrez.read(handle)
    idList = res["IdList"]

    # search for all GenBank results for the given IDs
    handle = Entrez.efetch(db='nucleotide', id=idList, rettype='gb', retmode='xml')
    records = Entrez.read(handle)
    return records


def write_records():
    records = fetch_records()
    dirname = os.path.dirname(__file__)
    results_file = os.path.join(dirname, 'A2_2_BRCA1_Kainz_Jakob.csv')

    # definies the parameters to organise the .csv file
    par_accesion = "GBSeq_primary-accession"
    par_title = "GBSeq_references"
    par_organism = "GBSeq_organism"
    par_length = "GBSeq_length"
    par_sequence = "GBSeq_sequence"

    # Creating the file and it´s header
    result_file = open(results_file, "w")
    result_lines = [
        "Accession number, Title, Organism, Sequence length, CG %, Protein instability, Aromaticity, Isoelectric point"]

    for record in records:
        # Calculating the additional values that are not included in the original record
        gc_percent = SeqUtils.GC(record[par_sequence])
        protein = Seq(record[par_sequence]).translate(to_stop=True)
        protein_analyse = ProteinAnalysis(str(protein))
        aromaticity = protein_analyse.aromaticity()
        protein_instability = protein_analyse.instability_index()
        isoelectric_point = protein_analyse.isoelectric_point()

        # Appending all info into a row of the .csv file
        result_lines.append(','.join([
            record[par_accesion],
            record[par_title][0]['GBReference_title'],
            record[par_organism],
            record[par_length],
            str(gc_percent),
            str(protein_instability),
            str(aromaticity),
            str(isoelectric_point)
        ]))
    result_file.write("\n".join(result_lines))
