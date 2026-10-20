# Imports
import math
import numpy as np
import time
import argparse

# Input
parser = argparse.ArgumentParser(description="PSSM DNA motif detection")

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Path to the genome FNA file"
)

args = parser.parse_args()

start_time = time.time()


# Parse Gene Sequences .fna
def parse_fna(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    sequence_lines = [line.strip() for line in lines if not line.startswith('>')]
    clean_seq = "".join(sequence_lines).upper()

    return clean_seq


input_sequences = parse_fna(args.input)


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
def create_ppm(pfm, sequences, k=1.0):
    # k --> pseudo-count (Laplace smoothing) to prevent log(0) domain errors during PSSM conversion
    # bg --> background frequency of each nucleotide in the target genome
    bg_freq = {'A': 0.285, 'T': 0.285, 'G': 0.215, 'C': 0.215}
    num_seqs = len(sequences)
    ppm = {}
    for base, counts in pfm.items():
        # Add (bg * k) to numerator and k to denominator
        bg = bg_freq[base]
        ppm[base] = [(count + bg * k) / (num_seqs + k) for count in counts]
    return ppm


ppm = create_ppm(pfm, sequences)


# Constructing Position-Specific Scoring Matrix (PSSM)
def create_pssm(ppm):
    """
    Converts probability of each base into log scores
    """
    bg_dict = {'A': 0.285, 'T': 0.285, 'G': 0.215, 'C': 0.215}
    pssm = {}

    for base, probs in ppm.items():
        bg = bg_dict[base]
        # calculate log2(prob / bg) for each position in the list
        pssm[base] = [math.log2(p / bg) for p in probs]

    return pssm


pssm = create_pssm(ppm)


# Motif Detection
def calc_threshold(pssm, target_p_value=0.001, precision=100):
    """
    Calculates threshold using dynamic programming to determine
    the probability of random DNA sequences achieving each possible score.
    Selects a threshold based on a specified maximum false-positive probability (0.001).
    """

    bg = {'A': 0.285, 'T': 0.285, 'G': 0.215, 'C': 0.215}
    bases = 'ACGT'

    dist = np.array([1.0])  # initialise distribution with 100% probability at score 0
    lowest = 0  # lowest possible cumulative PSSM score

    for i in range(len(pssm['A'])):
        score = {b: int(round(pssm[b][i] * precision)) for b in bases}
        min_score, max_score = min(score.values()), max(score.values())

        # initialise temporary array sized for the updated score range
        new = np.zeros(len(dist) + max_score - min_score)

        for b in bases:
            # calculate array index shift for the current nucleotide score
            idx_shift = score[b] - min_score
            new[idx_shift:idx_shift + len(dist)] += dist * bg[b]

        # update probability distribution & lowest possible cumulative score
        dist, lowest = new, lowest + min_score

    # calculates cumulative probability of getting each score or higher
    p_values = np.cumsum(dist[::-1])[::-1]  # P(score >= lowest + i)
    idx = np.flatnonzero(p_values <= target_p_value)  # find score indices with match probabilities at/below the target

    if idx.size == 0:
        # unreachable target: default to the max score
        return (lowest + len(dist) - 1) / precision

    # otherwise converts selected array index into the PSSM threshold score
    return (lowest + idx[0]) / precision


score = calc_threshold(pssm)


def scan_sequence(target_seq, pssm, threshold):
    """
    Scans a target sequence using a PSSM and returns matching windows and scores
    """
    motif_len = len(next(iter(pssm.values())))  # Get motif length
    results = []

    # Initialise sliding window across gene sequence length
    for i in range(len(target_seq) - motif_len + 1):
        window = target_seq[i:i + motif_len]
        score = 0.0
        valid_window = True

        # Score each base pair position inside the current window
        for j, base in enumerate(window):
            if base in pssm:  # Check if nucleotide exists in PSSM dictionary
                score += pssm[base][j]  # Add log score for this base at this position
            else:
                valid_window = False  # Handle unknown characters
                break

        if valid_window:
            if threshold is None or score >= threshold:
                results.append((i, window, score))

    return results  # Outputs --> (index, window, score)


motifs = scan_sequence(input_sequences, pssm, score)
end_time = time.time()

print("Motifs:", motifs)
print(f"Total time:     {end_time - start_time:.3f}s")
