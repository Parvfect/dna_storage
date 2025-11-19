
from Bio import SeqIO, SeqRecord
from tqdm import tqdm
import random
from typing import List


def reverse_complement(dna: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def parse_biopython(input_fastq):
    for record in SeqIO.parse(input_fastq, 'fastq'):
        yield record

def get_fastq_records(fastq_filepath):
    records = []
    for i, record in enumerate(tqdm(parse_biopython(fastq_filepath))):
        records.append(record)
    return records

def read_synthesized_strands_from_file(file_path, ids=False):

    sequences = []
    ids = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('>'):
            ids.append(line[1:].strip())

        if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            sequences.append(line.strip())

    if ids:
        return sequences, ids
    
    return sequences

def sample_reads(reads, ids, sampling_rate=None, n_samples=1000):
    """
    Sample reads and ids """

    if sampling_rate:
        n_samples = int(len(reads) * sampling_rate)

    sample_indices = [random.randint(0, len(reads)) for i in range(n_samples)]
    sampled_reads = [reads[i] for i in sample_indices]
    sampled_ids = [ids[i] for i in sample_indices]

    return sampled_reads, sampled_ids

def get_badread_strand_id(record):
    return record.description.split()[1].split(',')[0]

def create_fasta_file(ids: List[str], strands: List[str], output_filepath: str):
    with open(output_filepath, 'w') as f:
        for i, strand in enumerate(strands):
            f.write(f">{ids[i]}\n")
            f.write(strand + '\n\n')

    print(f"File saved as {output_filepath}")

def create_random_strand(strand_length: int) -> str:
    return "".join(
        [random.choice(['A', 'C', 'T', 'G']) for i in range(
            strand_length)])

## Removing strand pool primers

def remove_adapters_zero_error(strand_pool, strand_ids):
    
    front_adapter = 'CATGAAGTGGAGTAGCGGCG'
    reverse_adapter = 'GCGACCGATGCTCACTGATC'

    rev_front = reverse_complement(reverse_adapter)
    rev_back = reverse_complement(front_adapter)

    cleaned_strand_pool = []
    cleaned_strand_ids = []

    for ind, i in enumerate(strand_pool):
        
        if front_adapter in i[:20]:
            if reverse_adapter in i[-21:]:
                cleaned_strand = i[20: -21]
                
                if len(cleaned_strand) > 200 and len(cleaned_strand) < 230:
                    cleaned_strand_pool.append(cleaned_strand)
                    cleaned_strand_ids.append(strand_ids[ind])
        if rev_front in i[:20]:
            if rev_back in i[-21:]:
                cleaned_strand = reverse_complement(i[20:-21])
                if len(cleaned_strand) > 200 and len(cleaned_strand) < 230:
                    cleaned_strand_pool.append(cleaned_strand)
                    cleaned_strand_ids.append(strand_ids[ind])

    return cleaned_strand_pool, cleaned_strand_ids