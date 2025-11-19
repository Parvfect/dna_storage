
from utils import read_synthesized_strands_from_file
import Levenshtein
from clustering import Clustering
from utils import get_fastq_records
from strand_reconstruction import make_prediction
from tqdm import tqdm
from utils import reverse_complement
from crc_encoding import get_crc_strand
from jpeg_strand_encoding import save_partially_decoded_jpeg
import os

def validate_crc(strand, info_length=1113):
    return get_crc_strand(strand[:info_length]) == strand

original_strands = read_synthesized_strands_from_file('data/final_run/bird_strands.fasta')[0]

records = get_fastq_records(r'data/final_run/birding.fastq')
strands = [str(i.seq) for i in records]
ids = [i.id for i in records]
strand_length = 1129

strands_length_filtered = [i for i in strands if len(i) > strand_length - 5 and len(i) < strand_length + 10]


clustering_obj = Clustering(strand_pool=strands_length_filtered, reference_length=1128, n_reference_strands=6, distance_threshold=100)

clustering_obj.run_pipeline(fix_orientation=True)


def generate_candidates_crc_validated(
        clustered_seqs, n_clusters=15, n_attempts=5,
        strand_length=200, ma_sample_size=10):
    
    validated_strands = []
    for ind, i in tqdm(enumerate(clustered_seqs[:n_clusters]), total=n_clusters):  # Iterate through clusters
        # Can be RC remember!
        for k in range(n_attempts): # Repeat n_attempts time
            candidate = make_prediction(i, sample_size=ma_sample_size)  # Make candidate prediction
            rev = reverse_complement(candidate)  # Obtain RC vector
            if validate_crc(candidate):  # Validate CRC code for forward and rc prediction
                validated_strands.append(candidate)
                break
            elif validate_crc(rev):
                validated_strands.append(rev)
                break
            else:
                continue
    
    validated_strands = list(set(validated_strands))
    print(f"{len(validated_strands)} valid strands found")

    return validated_strands

validated_strands = generate_candidates_crc_validated(clustering_obj.clustered_seqs)

no_crc_original_strands = [i[20:-16] for i in original_strands]
t = [i[20:-16] for i in validated_strands]

print(len(set(t).intersection(
    no_crc_original_strands)))

ids = ["GGGATTTAGTGACATAATCG", "GGAGACGAGCCACGTGTCAT", "AAACGAGAAGGCCTTGCCCA", "TGAGGCAAGCTCCACCACGG", "CTCCCTCCCCGCAGCATATT", "CAGAGCTGAGTAGAGATCCA"]

id_length = 20

# Build mapping from ID → strand (assuming one strand per ID)
id_to_strand = {s[:id_length]: s for s in validated_strands}

savepath = 'data/final_run'

for k in range(len(ids)):
    # recover strands whose IDs appear in ids[0:k+1], in that exact order
    recovered_strands = [
        id_to_strand[i] for i in ids[:k+1] if i in id_to_strand
    ]

    save_partially_decoded_jpeg(
        recovered_strands=recovered_strands,
        n_strands=6,
        id_length=id_length,
        strand_ids=ids,
        filename=os.path.join(savepath, f"decoded_partial_{k}.jpg"),
        crc_encoding=True
    )
