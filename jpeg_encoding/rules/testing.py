

from utils import read_synthesized_strands_from_file
from jpeg_encoding.rules.FastDNARules import FastDNARules

x = FastDNARules()

def test_strand(strand):
    return x.apply_all_rules(strand)