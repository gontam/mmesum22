from Bio import Entrez
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.SeqUtils import GC
from Bio import Seq
from Bio import SeqRecord
import argparse
from csv import DictWriter
from matplotlib import pyplot as plt

# Default Gene: IL2RA (NC_000010.11)
DEFAULT_EMAIL = "me21m003@technikum-wien.at"
DEFAULT_NAME = "IL2RA"
DEFAULT_SPECIES_COUNT = 5
SEARCH_COUNT = 100
DEFAULT_OUTPUT_CSV = "A4_Hauptfeld_Leonhard.csv"
DEFAULT_OUTPUT_PNG = "A4_Hauptfeld_Leonhard.png"
OUTPUT_HEADERS = ['Accession number','Title','Organism','Sequence length','CG %']

def download_nucleotides(term, organism_count):
    # Look for results in the nucleotide DB
    print(f"\t- Searching for: {term}")
    handle = Entrez.esearch(db="nucleotide", retmax=SEARCH_COUNT, rettype="gb", term=f"{term}[Gene Name]")
    search_results = Entrez.read(handle)
    handle.close()
    # Download 5 nucleotides of unique organisms
    nucleotides = []
    for id in search_results['IdList']:
        # Download nucleotide
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
        nucleotide = SeqIO.read(handle, 'genbank')
        handle.close()
        # Skip if we already have a nucleotide from this species
        if nucleotide.annotations['organism'] in [n.annotations['organism'] for n in nucleotides]:
            print("\t- Skipping duplicate species...")
            continue
        # Verify we have sequence contents (hacky, but no other way found?!)
        try:
            nucleotide.seq.join([])
        except Seq.UndefinedSequenceError:
            print("\t- Skipping nucleotide with no sequence content...")
            continue
        # Add nucleotide
        print(f"\t- Species: {nucleotide.annotations['organism']}, Length: {len(nucleotide)}")
        nucleotides.append(nucleotide)
        # If we've reached target count, end
        if len(nucleotides) >= organism_count: break
    return nucleotides

def analyse_nucleotide(nucleotide):
    return {
        'Accession number': nucleotide.annotations['accessions'][0],
        'Title': nucleotide.description,
        'Organism': nucleotide.annotations['organism'],
        'Sequence length': len(nucleotide),
        'CG %': GC(nucleotide.seq)
    }
def save_nucleotides(nucleotides, filename):
    # Write to csv
    with open(filename, 'w') as csvfile:
        writer = DictWriter(csvfile, fieldnames = OUTPUT_HEADERS)
        writer.writeheader()
        writer.writerows([analyse_nucleotide(n) for n in nucleotides])

def align_nucleotides(nucleotides):
    # Pad all the sequences to a common length, convert to sequence records
    # Some parts taken from https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
    max_len = max(len(n.seq) for n in nucleotides)
    seq_records = []
    print(f"\t- Padding all sequences to common length ({max_len})...")
    for nucleotide in nucleotides:
        if len(nucleotide.seq) != max_len:
            seq = str(nucleotide.seq).ljust(max_len, '.')
            seq_records.append(SeqRecord.SeqRecord(seq, id=nucleotide.annotations['accessions'][0], name=nucleotide.annotations['organism']))
    # Use Bio.Align to align sequences
    return MultipleSeqAlignment(seq_records)
def construct_phylo_tree(alignment):
    calculator = DistanceCalculator()
    # print(calculator.get_distance(alignment))
    # distance_model = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    return tree

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='Acquire nucleotides for a gene from n species and run analysis')
    parser.add_argument('--name', type=str, default=DEFAULT_NAME, help='Name of gene')
    parser.add_argument('--species-count', type=int, default=DEFAULT_SPECIES_COUNT, help='Count of species to analyse')
    parser.add_argument('--email', type=str, default=DEFAULT_EMAIL, help='E-Mail address to identify with GenBank')
    parser.add_argument('--output-csv', type=str, default=DEFAULT_OUTPUT_CSV, help='Default CSV output file name')
    parser.add_argument('--output-png', type=str, default=DEFAULT_OUTPUT_PNG, help='Default PNG output file name')
    args = parser.parse_args()
    # Prepare Entrez
    Entrez.email = args.email
    # Acquire nucleotide sequences
    print("Searching for nucleotides...")
    nucleotides = download_nucleotides(args.name, args.species_count)
    print(f'{len(nucleotides)} nucleotides identified')
    # Save nucleotides
    print("Saving CSV...")
    save_nucleotides(nucleotides, args.output_csv)
    # Align nucleotides
    print("Aligning sequences...")
    alignment = align_nucleotides(nucleotides)
    # Create Tree
    print("Constructing phylogenetic tree...")
    tree = construct_phylo_tree(alignment)
    tree.ladderize()
    print("Saving tree...")
    Phylo.draw(tree,label_func=lambda cl: f'{str(cl)} ({round(cl.branch_length, 2)})',do_show=False)
    plt.savefig(args.output_png)
    print("Done.")