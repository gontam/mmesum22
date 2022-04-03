import argparse

# Constants
START_CODON = "M"
STOP_CODON = "*"

class TripletError:
    pass
class TranslationError:
    pass
class TranslationTriplet:
    def __init__(self, dna, protein, name):
        self.dna = dna
        self.protein = protein
        self.name = name

def parse_triplets(filename):
    triplets = []
    f = open(filename, "r")
    for l in f:
        line = l.rstrip().split(" ")
        if len(line) != 3: raise TripletError("Invalid triplet file format")
        if len(line[0]) != 3: raise TripletError("Triplets need to be exactly 3 characters long")
        if len(line[1]) != 1: raise TripletError("Translation proteins need to be exactly 1 character long")
        triplets.append(TranslationTriplet(line[0], line[1], line[2]))
    f.close()
    return triplets

class ProteinOpenFrame:
    def __init__(self, frame, start_index=None, stop_index=None):
        self.frame = frame
        self.start_index = start_index
        self.stop_index = stop_index
    def get_sequence(self):
        return self.frame.sequence[self.start_index:self.stop_index]
    def get_length(self):
        return self.stop_index - self.start_index
    def format_sequence(self, format):
        if format == "Met":
            return " ".join(["Met" if p == "M" else "Stop" if p == "*" else p for p in list(self.get_sequence())])
        elif format == "M":
            return self.get_sequence().replace("*", "-")
        else:
            raise ValueError(f"Invalid open frame format: {format}")
class ProteinFrame:
    def __init__(self, offset, reversed=False, sequence=""):
        self.offset = offset
        self.reversed = reversed
        self.sequence = sequence
    def append_sequence(self, sequence):
        self.sequence += sequence
    def read_open_frames(self):
        open_frames = []
        open_frame = None
        protein_index = 0
        for protein in self.sequence:
            if protein == STOP_CODON and open_frame is not None:
                open_frame.stop_index = protein_index + 1  # Include the stop codon
                open_frames.append(open_frame)
                open_frame = None
            elif protein == START_CODON and open_frame is None:
                open_frame = ProteinOpenFrame(self, protein_index)
            protein_index += 1
        return open_frames

class FastaError(Exception):
    pass
class FastaSequence:
    def __init__(self, description, sequence=""):
        self.description = description
        self.sequence = sequence
    def append_sequence(self, sequence):
        self.sequence += sequence
    def translate(self, triplets, reverse=False):
        # Copy stored sequence, reversed or not
        sequence = self.sequence[::-1] if reverse else self.sequence
        # Reorganize triplets into dictionary for easier access
        triplets = { t.dna: t for t in triplets }
        frames = { i: ProteinFrame(i, reverse) for i in range(0,3) }
        for i in range(0, len(sequence), 3):
            # Get the next three characters in the three frames (offsets)
            for offset in range(0,3):
                dna = sequence[i+offset:i+offset+3]
                # If a problem comes along, you just skip it (skip it good)
                if len(dna) != 3: continue
                # Try to find the translation sequence
                if dna not in triplets: raise TranslationError(f"Sequence {dna} not found for translation")
                triplet = triplets[dna]
                # Translate
                frames[offset].append_sequence(triplet.protein)
        return list(frames.values())

def parse_fasta(filename):
    sequences = []
    sequence = None
    f = open(filename, "r")
    for l in f:
        line = l.rstrip()
        if line[0] == ">":
            if sequence is not None: sequences.append(sequence)
            sequence = FastaSequence(description=line[1:])
        elif line[0] == ";": continue
        else:
            if sequence is None: raise FastaError("Headerless sequence")
            sequence.append_sequence(line)
    f.close()
    # Append any remaining last sequence
    if sequence is not None: sequences.append(sequence)
    return sequences

# Parse Arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse a FASTA formatted sequence and write the largest open reading frame to a file')
    parser.add_argument('--input', type=str, default="BD137219.1.fasta", help='File in FASTA format to read')
    parser.add_argument('--output', type=str, default="BD137219.1_translated_Hauptfeld.txt", help='File to write the open reading frame to')
    parser.add_argument('--format', type=str, default="Met", choices=["M", "Met"], help='Format to output the open reading frame in')
    parser.add_argument('--triplets', type=str, default="standard.txt", help='File with translation pairs')
    parser.add_argument('--no-reverse', action='store_true', help="Don't parse the sequence in reverse")
    parser.add_argument('--sequence', type=int, default=0, help='Index of sequence from FASTA file to parse')
    args = parser.parse_args()
    # Parse protein translations
    triplets = parse_triplets(args.triplets)
    print(f"{len(triplets)} triplets found")
    # Parse sequence file
    fasta_sequences = parse_fasta(args.input)
    print(f"{len(fasta_sequences)} sequences found:")
    for s in fasta_sequences:
        print(s.description)
    # Translate, both forward and reverse
    protein_frames = fasta_sequences[args.sequence].translate(triplets, reverse=False)
    if not args.no_reverse:
        protein_frames.extend(fasta_sequences[args.sequence].translate(triplets, reverse=True))
    print(f"{len(protein_frames)} frames translated")
    # Identify open frames
    protein_open_frames = []
    for f in protein_frames:
        protein_open_frames.extend(f.read_open_frames())
    print(f"{len(protein_open_frames)} open frames identified")
    # Identify longest frame
    longest_open_frame = max(protein_open_frames, key=lambda orf: orf.get_length())
    print(f"Longest frame length is {longest_open_frame.get_length()}")
    # Dump to file in requested format
    f = open(args.output, "w")
    f.write(longest_open_frame.format_sequence(args.format))
    f.close()
    print(f"Dumped result in {args.format} format to {args.output}")