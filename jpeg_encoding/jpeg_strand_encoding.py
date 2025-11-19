
from jpeg_encoding.crc_encoding import get_crc_strand
from jpeg_encoding.xor import get_valid_xor_seed_and_strand, get_bitstream_from_strand, check_obeys_dna_constraints
import numpy as np
from textwrap import wrap
import math
import random
import os
from itertools import cycle
from datetime import datetime
import json
from datetime import datetime


bits_to_base = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
alphabet = ["A", "C", "G", "T"]

def add_xor_mask(bitstream, seed='01001'):
    return ''.join(
        '1' if b != mask else '0' for b, mask in zip(
            bitstream, cycle(seed)))

def index_to_bases(idx: int, length: int = 8) -> str:
    bits = format(idx, f'0{2*length}b')  # 2 bits per base
    return ''.join(bits_to_base[bits[i:i+2]] for i in range(0, len(bits), 2))

def map_bits_ternary(bits):
    # convert bitstream to integers 0,1,2
    ints = [int(bits[i:i+2], 2) % 3 for i in range(0, len(bits), 2)]

    dna = []
    prev = None

    for x in ints:
        # available bases excluding previous base
        allowed = [b for b in alphabet if b != prev]
        base = allowed[x % len(allowed)]
        dna.append(base)
        prev = base

    return ''.join(dna)


def encode_strands(
        input_file_jpeg, id_length, strand_length,fix_strand_ids=True,
        crc_encoding=False, timestamp=None):
    """
    Encodes strands by adding XOR mask of 0101 and basic bit to dna conversion
    Splits to length such that its equally balanced
    """

    fasta_file = "bird_strands.fasta"
    bitstream_strands = "bird_bit_strands.fasta"
    cfg_file = "cfg.json"
    preceding_path = 'jpeg_encoding/data/'
    if os.path.isdir(preceding_path):
        new_path = os.path.join(preceding_path, timestamp)
        os.makedirs(new_path)

        fasta_file = os.path.join(new_path, fasta_file)
        cfg_file = os.path.join(new_path, cfg_file)
        bitstream_strand_file = os.path.join(new_path, bitstream_strands)

    with open(input_file_jpeg, "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)

    bitstream = ''.join(format(b, '08b') for b in data)
    ptr = 0
    dna = ''.join(bits_to_base[bitstream[i:i+2]] for i in range(
        0, len(bitstream), 2))

    payload_length = strand_length
    total_bases = len(dna)
    n_strands = math.ceil(total_bases / strand_length)

    # Adjust if too long for the payload section
    bases_per_strand = min(n_strands, payload_length)
    chunks = wrap(dna, strand_length)        

    with open(bitstream_strand_file, "w") as f:
        for i, strand in enumerate(chunks):
            f.write(f">strand_{i}\n{strand}\n")


    new_strands = []
    xor_seeds = []
    for strand in chunks:
        xor_strand, xor_seed = get_valid_xor_seed_and_strand(strand)
        new_strands.append(xor_strand)
        xor_seeds.append(xor_seed)
    
    strand_ids = xor_seeds

    # Add primers (truncate or pad as needed)
    strands = [
        f"{strand_ids[i]}{strand}" for i, strand in enumerate(
            new_strands)]
    
    padding = np.zeros(n_strands)
    padded_strands = []
    strand_length = len(strands[0])

    for ind, strand in enumerate(strands):
        padding_length = 0
        if len(strand) < strand_length: # Padding
            padding_length = strand_length - len(strand)
            strand += "A" * (padding_length)

        padded_strands.append(strand)
        padding[ind] = padding_length
    
    strands = padded_strands

    if crc_encoding:
        strands = [
            get_crc_strand(strand) for strand in strands]
    

    with open(fasta_file, "w") as f:
        for i, strand in enumerate(strands):
            f.write(f">strand_{strand_ids[i]}\n{strand}\n")

    cfg_dict = {
        'date': str(datetime.today()),
        'image_encoded': input_file_jpeg,
        'fasta_encoded': fasta_file,
        'XOR seed': xor_seed,
        'strand_ids': strand_ids,
        'padding': list(padding)
    }

     # ==== WRITE config file ====
    with open(cfg_file, "w") as f:
        json.dump(cfg_dict, f)

    print(f"Encoded {len(strands)} strands written to {fasta_file}")
    return strands, strand_ids, total_bases, n_strands, padding


def save_partially_decoded_jpeg(
        recovered_strands, n_strands,
        id_length, strand_ids, filename="decoded_partial.jpg",
        maintain_order=True, crc_encoding=True,
        padding=[0, 0, 0, 0, 0, 2], timestamp=None):

    # Remove primers and join payloads safely
    recovered_payloads = []
    strand_order = []
    strand_ids = []
    for ind, s in enumerate(recovered_strands):
        if len(s) > id_length:
            recovered_strand = s[id_length:]
            recovered_id = s[:id_length]

            if crc_encoding:
                recovered_strand = recovered_strand[:-16]
            
            recovered_payloads.append(recovered_strand)
            strand_ids.append(recovered_id)

            if recovered_id in strand_ids:
                strand_order.append(strand_ids.index(recovered_id))
        else:
            print("Warning: strand too short, skipping:", s[:20], "…")

    # Sort recovered_strands by strand_order - reject if not strictly correct
    
    if sorted(strand_order) != list(np.arange(n_strands))[:len(recovered_payloads)]:
        print('Strands are in incorrect order! Decoding incomplete')
        return False
    
    if not recovered_payloads:
        raise ValueError("No valid strand payloads found")
    
    if len(recovered_payloads) > 1:
        indices = np.argsort(strand_order)
        recovered_payloads = np.array(recovered_payloads)[
            indices]
        strand_ids = np.array(strand_ids)[indices]

    if padding is not None:
        recovered_payloads_ = []

        for ind, i in enumerate(recovered_payloads):
            strand = i
            if padding[ind] > 0:
                strand = i[:-int(padding[ind])]
            recovered_payloads_.append(strand)

        recovered_payloads = recovered_payloads_

    recovered_payload = ''.join(recovered_payloads)
    #recovered_payload = recovered_payload[:total_bases]  # trim to max length if needed

    base_to_bits = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}

    filtered_payload = ''.join(
        b for b in recovered_payload if b in base_to_bits)
    
    recovered_bits = ''.join(
        get_bitstream_from_strand(p, id) for p, id in zip(recovered_payloads, strand_ids)
    )

    bit_len = len(recovered_bits) - (
        len(recovered_bits) % 8)
    if bit_len == 0:
        raise ValueError("Recovered bits too short for even one byte")

    recovered_bits = recovered_bits[:bit_len]

    byte_data = int(recovered_bits, 2).to_bytes(bit_len // 8, byteorder='big')

    decoded_file = filename
    print(decoded_file)

    if os.path.isdir('data'):
        preceding_path = 'jpeg_encoding/data/'
        new_path = os.path.join(preceding_path, timestamp)
        decoded_file = os.path.join(new_path, decoded_file)

    with open(decoded_file, "wb") as f:
        f.write(byte_data)

    return True  # Would like a more robust test here


if __name__ == '__main__':

    id_length = 20
    dir = os.path.join("data", "bird")

    decoded_file = "bird_decoded.jpg"
    cfg_file = "cfg.json"
    xor_seed = '01001'
    input_file_jpeg = r"C:\Users\Parv\Doc\RA\Projects\dna_storage\jpeg_encoding\data\bird.jpg"

    strand_length = 1093          # total bases per strand (including primer)
    primer_prefix = "ACGTACGTACGT" # 20nt primer (easy to read/recognize)
    crc_encoding=True

    timestamp = str(datetime.now()).replace('-', '_').replace(':','_')
    # Not encoding anymore
    
    strands, strand_ids, total_bases, n_strands, padding = encode_strands(
        input_file_jpeg=input_file_jpeg, id_length=id_length,
        strand_length=strand_length, crc_encoding=crc_encoding, timestamp=timestamp)
    
    for i in range(n_strands):
        save_partially_decoded_jpeg(
            recovered_strands=strands[:i+1],
            n_strands=n_strands, id_length=id_length,
            strand_ids=strand_ids, filename=f"decoded_partial_{i}.jpg",
            crc_encoding=crc_encoding, padding=padding, timestamp=timestamp)