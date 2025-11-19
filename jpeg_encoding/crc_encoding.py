
import random
import zlib

base_mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
int_mapping = {'00': 'A', '01': 'C', "10": 'G', "11": 'T'}

def convert_dna_to_binary_string(strand):
    return "".join([base_mapping[i] for i in strand])

def convert_binary_string_to_dna(binary_str):
    return "".join(int_mapping[binary_str[i:i+2]] for i in range(
        0, len(binary_str), 2))

def convert_dna_to_byte_string(strand):
    binary_str = convert_dna_to_binary_string(strand)

    padding = (8 - len(binary_str) % 8) % 8
    binary_str = binary_str + "0" * padding  # pad with zeros at the end
    
    byte_string = int(binary_str, 2).to_bytes(
        len(binary_str) // 8, byteorder="big")

    return byte_string

def generate_random_strand(strand_length: int):
    return "".join([random.choice(
        ['A', 'C', 'T', 'G']) for i in range(strand_length)])

def convert_integer_to_dna(num):
    binary_str = bin(num)[2:].zfill(32)
    if len(binary_str) % 2 != 0:
        binary_str = "0" + binary_str
    return convert_binary_string_to_dna(binary_str)

def get_crc_strand(strand):
    crc = zlib.crc32(convert_dna_to_byte_string(strand))
    dna_crc = convert_integer_to_dna(crc)
    final_strand = strand + dna_crc
    return final_strand