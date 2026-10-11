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


# Create IUPAC
def ppm_to_iupac(ppm, threshold=0.2):
    """
    Takes a PPM dictionary  and generates an IUPAC degenerate consensus sequence string
    using a frequency threshold
    """

    # IUPAC lookup dictionary
    iupac_mapping = {
        ('A',): 'A',
        ('C',): 'C',
        ('G',): 'G',
        ('T',): 'T',
        ('A', 'G'): 'R',  # Purine
        ('C', 'T'): 'Y',  # Pyrimidine
        ('C', 'G'): 'S',  # Strong
        ('A', 'T'): 'W',  # Weak
        ('G', 'T'): 'K',  # Keto
        ('A', 'C'): 'M',  # Amino
        ('A', 'C', 'G', 'T'): 'N'  # Any
    }

    # Get motif length from the first row of PPM
    motif_len = len(next(iter(ppm.values())))
    consensus_chars = []

    for i in range(motif_len):
        valid_bases = []

        # Check each nucleotide's probability at column position i
        for base in ['A', 'C', 'G', 'T']:
            if ppm[base][i] >= threshold:
                valid_bases.append(base)

        # Sort the bases alphabetically to ensure they match the lookup dict order
        valid_bases.sort()
        base_tuple = tuple(valid_bases)

        # Assign IUPAC letter
        consensus_char = iupac_mapping.get(base_tuple, 'N')
        consensus_chars.append(consensus_char)

    return "".join(consensus_chars)


iupac_seq = ppm_to_iupac(ppm)


def scan_iupac_sequence(target_seq, iupac_consensus):
    """
    Scans a target sequence using a sliding window and an IUPAC consensus string
    """
    consensus = iupac_consensus.upper()
    motif_len = len(consensus)
    results = []

    iupac_rules = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
        'R': {'A', 'G'},   # Purine
        'Y': {'C', 'T'},   # Pyrimidine
        'S': {'G', 'C'},   # Strong
        'W': {'A', 'T'},   # Weak
        'K': {'G', 'T'},   # Keto
        'M': {'A', 'C'},   # Amino
        'N': {'A', 'C', 'G', 'T'}  # Any
    }

    # Slide the window across the target sequence
    for i in range(len(target_seq) - motif_len + 1):
        window = target_seq[i:i + motif_len]
        valid_window = True

        # Check each base pair position against the corresponding IUPAC character
        for j, base in enumerate(window):
            consensus_char = consensus[j]
            allowed_bases = iupac_rules.get(consensus_char, {'A', 'C', 'G', 'T'})

            # Handle unknown characters
            if base not in allowed_bases:
                valid_window = False
                break

        if valid_window:
            results.append((i, window))

    return results  # Outputs --> (index, window).


motifs = scan_iupac_sequence(input_sequences, iupac_seq)
