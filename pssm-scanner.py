# Imports
import math


# Parse Gene Sequences .fna
def parse_fna(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    sequence_lines = [line.strip() for line in lines if not line.startswith('>')]
    clean_seq = "".join(sequence_lines).upper()

    return clean_seq


input_sequences = parse_fna("gene.fna")


# Parse Binding Sites .fasta
def parse_fasta(filepath):
    sequences = []
    prev_header = ''
    current_sequence = ''

    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if prev_header:
                    sequences.append((prev_header, current_sequence))

                prev_header = line[1:]
                current_sequence = ''

            else:
                current_sequence += line

        if prev_header:
            sequences.append((prev_header, current_sequence))

    return sequences


sites = parse_fasta("MA0011.2.fasta")


# Input Validation & Normalisation
def normalise(sequences):
    normalised_sequences = []
    current_sequence = ''

    for header, sequence in sequences:
        for base in sequence:
            if not base.isspace():
                current_sequence += base.upper()
        normalised_sequences.append((header, current_sequence))
        current_sequence = ''

    return normalised_sequences


sequences = normalise(sites)


def validate(sequences):
    valid_sequences = []
    valid_bases = "ACGT"

    for header, sequence in sequences:
        valid = True
        for base in sequence:
            if base not in valid_bases:
                valid = False
        if valid:
            valid_sequences.append((header, sequence))

    return valid_sequences


sequences = validate(sequences)


# Constructing Position Frequency Matrix (PFM)
def create_pfm(sequences):
    length = len(sequences[0][1])
    pfm = {base: [0] * length for base in 'ACGT'}

    for _, seq in sequences:
        for i, base in enumerate(seq):
            pfm[base][i] += 1

    return pfm


pfm = create_pfm(sequences)


# Constructing Position Probability Matrix (PPM)
def create_ppm(pfm, sequences, k=1.0, bg=0.25):
    # k --> pseudo-count (Laplace smoothing) to prevent log(0) domain errors during PSSM conversion
    # bg --> background frequency of each nucleotide in the target genome
    num_seqs = len(sequences)
    ppm = {}
    for base, counts in pfm.items():
        # Add (bg * k) to numerator and k to denominator
        ppm[base] = [(count + bg * k) / (num_seqs + k) for count in counts]
    return ppm


ppm = create_ppm(pfm, sequences)


# Constructing Position-Specific Scoring Matrix (PSSM)
def create_pssm(ppm):
    """
    Converts probability of each base into log scores
    """
    bg_dict = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    pssm = {}

    for base, probs in ppm.items():
        bg = bg_dict[base]
        # calculate log2(prob / bg) for each position in the list
        pssm[base] = [math.log2(p / bg) for p in probs]

    return pssm


pssm = create_pssm(ppm)
