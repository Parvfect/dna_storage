
from crc_encoding import \
    convert_dna_to_binary_string, convert_binary_string_to_dna
from itertools import cycle
import random
from utils import read_synthesized_strands_from_file
from tqdm import tqdm
from collections import Counter

def kmer_repeats(dna, k=3, window=20, max_count=2):
    bad_regions = 0

    for start in range(0, len(dna) - window + 1):
        w = dna[start:start+window]

        # extract k-mers in this window
        kmers = [w[i:i+k] for i in range(len(w) - k + 1)]
        counts = Counter(kmers)

        # find over-represented kmers inside this window
        local_bad = [(km, c) for km, c in counts.items() if c > max_count]
        if local_bad:
            bad_regions += 1

    return bad_regions == 0


def check_obeys_dna_constraints(
        strand, max_hp_length=3, max_gc_content=0.7):
    """
    Checks if homopolymers are less than 3 and 
    GC content is less than 45 %
    """

    run = 1
    gc = 0
    prev_nt = ""
    violations = 0
    for nt in strand:
        if nt == prev_nt:
            run += 1
        else:
            run = 1

        if nt == 'G' or nt == 'C':
            gc += 1
        
        if run > max_hp_length:
            return False
        prev_nt = nt
        
    
    gc_content = gc / len(strand)

    if gc_content > max_gc_content:
        return False
    
    if not (kmer_repeats(strand, k=4) or kmer_repeats(strand, k=3)):
        return False

    print(f"GC content is {gc_content}")
    
    return True


def generate_random_xor_seed(seed_length):
    return "".join([
        str(random.randint(0, 1)) for i in range(seed_length)])


def add_xor_mask(bitstream, seed=None):

    if seed is None:
        seed = generate_random_xor_seed(128)

    return ''.join(
        '1' if b != mask else '0' for b, mask in zip(
            bitstream, cycle(seed)))


def get_valid_xor_seed_and_strand(strand, xor_seed_length=20):
    
    xor_seed_dna = ""
    final_strand = ""
    for i in range(4):
        start = i * 300
        end = (i + 1) * 300
        if end > len(strand):
            end = len(strand)
        substrand = strand[start: end]
        bin_strand = convert_dna_to_binary_string(
            substrand)
        while True:
            xor_seed = generate_random_xor_seed(xor_seed_length)
            xored_bin_string = add_xor_mask(bin_strand, xor_seed)
            new_strand = convert_binary_string_to_dna(xored_bin_string)
            if check_obeys_dna_constraints(new_strand):
                final_strand += new_strand
                xor_seed_dna += convert_binary_string_to_dna(xor_seed)
                break

    return final_strand, xor_seed_dna


def get_bitstream_from_strand(strand, xor_seed):

    fixed_bitstream = ""
    step_length = len(xor_seed) // 4
    
    for i in range(4):
        start = i * 300
        end = (i + 1) * 300
        if end > len(strand):
            end = len(strand)
        substrand = strand[start: end]
        xor_seed_sub = xor_seed[
            step_length * i : step_length * (i + 1)]
        
        bin_strand = convert_dna_to_binary_string(
            substrand)
        bin_xor_subseed = convert_dna_to_binary_string(xor_seed_sub)

        t = add_xor_mask(bin_strand, bin_xor_subseed)
        fixed_bitstream += add_xor_mask(bin_strand, bin_xor_subseed)

    return fixed_bitstream


