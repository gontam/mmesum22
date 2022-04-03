from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from csv import DictWriter

SEARCH_TERMS = ["HLA-DRB1", "KIR2DS1", "IFN-g"]
SEARCH_COUNTS = [2,2,1]
SEARCH_ORGANISM = "Homo sapiens"
OUTPUT_FILE = "A2.2_MultipleSclerosis_Hauptfeld_Leonhard.csv"
OUTPUT_HEADERS = ['Accession number','Title','Organism','Sequence length','CG %','Protein instability','Aromaticity','Isoelectric point']

Entrez.email = "me21m003@technikum-wien.at"
        
def retrieve_gene(id):
    # Retrieve details from nucleotide DB
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    rec = SeqIO.read(handle, 'genbank')
    handle.close()
    print(f"\t- {rec.annotations['accessions'][0]}")
    # Do an analysis
    analysis = ProteinAnalysis(str(rec.translate(to_stop=True).seq))
    # Return stuff we want in a dict
    return {
        'Accession number': rec.annotations['accessions'][0],
        'Title': rec.description,
        'Organism': rec.annotations['organism'],
        'Sequence length': len(rec),
        'CG %': GC(rec.seq),
        'Protein instability': analysis.instability_index(),
        'Aromaticity': analysis.aromaticity(),
        'Isoelectric point': analysis.isoelectric_point()
    }
def search_gene(term, organism, result_count=2):
    # Look for results in the nucleotide DB
    print(f"{term}")
    handle = Entrez.esearch(db="nucleotide", retmax=result_count, rettype="gb", term=f"{term}[All Fields] AND \"{organism}\"[Organism]")
    search_results = Entrez.read(handle)
    handle.close()
    return search_results['IdList'][:result_count]

if __name__ == "__main__":
    # Gather Gene info
    index = 0
    info = []
    for term in SEARCH_TERMS:
        gene_ids = search_gene(term, SEARCH_ORGANISM, SEARCH_COUNTS[index])
        for gene_id in gene_ids:
            info.append(retrieve_gene(gene_id))
        index += 1
    # Write to csv
    with open(OUTPUT_FILE, 'w') as csvfile:
        writer = DictWriter(csvfile, fieldnames = OUTPUT_HEADERS)
        writer.writeheader()
        writer.writerows(info)

