from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import UndefinedSequenceError
import argparse
from csv import DictWriter

# Default Gene: IL2RA (NC_000010.11)
DEFAULT_NAME = "IL2RA"
DEFAULT_SPECIES_COUNT = 1
SEARCH_COUNT = 100
DEFAULT_OUTPUT_CSV = "A4_Hauptfeld_Leonhard.csv"
OUTPUT_HEADERS = ['Accession number','Title','Organism','Sequence length','CG %']

Entrez.email = "me21m003@technikum-wien.at"

def download_nucleotides(term, organism_count):
    # Look for results in the nucleotide DB
    print(f"Searching for: {term}")
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
            print("Skipping duplicate species...")
            continue
        # Verify we have sequence contents (hacky, but no other way found?!)
        try:
            print(nucleotide.seq)
        except UndefinedSequenceError:
            print("Skipping nucleotide with no sequence content...")
            continue
        # Add nucleotide
        print(f"Species: {nucleotide.annotations['organism']}, Length: {len(nucleotide)}")
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

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='Acquire nucleotides for a gene from n species and run analysis')
    parser.add_argument('--name', type=str, default=DEFAULT_NAME, help='Name of gene')
    parser.add_argument('--species-count', type=int, default=DEFAULT_SPECIES_COUNT, help='Count of species to analyse')
    parser.add_argument('--output-csv', type=str, default=DEFAULT_OUTPUT_CSV, help='Default CSV output file name')
    args = parser.parse_args()
    # Acquire nucleotide sequences
    nucleotides = download_nucleotides(args.name, args.species_count)
    print(f'{len(nucleotides)} nucleotides identified')
    print("Saving CSV...")
    # Save nucleotides
    save_nucleotides(nucleotides, args.output_csv)