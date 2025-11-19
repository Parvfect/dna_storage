
from clustering.pool_preprocessing import remove_adapters_from_strands
from clustering.heirarchal_clustering import cluster_strands
from clustering.strand_reconstruction import get_clustered_seqs, get_candidates, get_candidate_orientation
from clustering.evaluation import evaluate_candidates
from typing import List, Tuple, Dict
from utils import reverse_complement
import numpy as np
import matplotlib.pyplot as plt


class Clustering:

    def __init__(
            self, strand_pool: List[str], reference_length: int, original_strands: List[str] = [],
            strand_pool_ids: List[str] = [], original_strand_ids: List[str] = None, distance_threshold: int = 40,
            front_adapter: str = "AATGTACTTCGTTCAGTTACGTATTGCT", reverse_adapter=None, n_reference_strands: int = 0):

        self.strand_pool = strand_pool
        self.n_strands_pool = len(strand_pool)
        self.reference_length = reference_length
        self.original_strands = original_strands
        if n_reference_strands == 0:
            self.n_reference_strands = len(original_strands)
        else:
            self.n_reference_strands = n_reference_strands
        self.strand_pool_ids = strand_pool_ids
        self.original_strands_ids = original_strand_ids
        self.distance_threshold = distance_threshold
        self.front_adapter = front_adapter
        self.reverse_adapter = reverse_adapter

    def filter_by_length(
            self, max_length:int = 50, min_length: int = 5, ids: bool = False) -> List[str]:

        filtered_indices = [ind for ind in range(self.n_strands_pool) if len(
            self.strand_pool[ind]) < self.reference_length + max_length and len(
                self.strand_pool[ind]) > self.reference_length - min_length]
                
        strand_pool = [self.strand_pool[ind] for ind in filtered_indices]
        self.strand_pool = strand_pool

        print(
            f"{(self.n_strands_pool - len(self.strand_pool))/ self.n_strands_pool} strands filtered out")
        self.n_strands_pool = len(self.strand_pool)
        
        if ids:
            strand_pool_ids = [self.strand_pool_ids[ind] for ind in filtered_indices]
            self.strand_pool_ids = strand_pool_ids

            return self.strand_pool, self.strand_pool_ids
        
        return self.strand_pool
    
    def remove_adapters(
            self, front_adapter: str = "AATGTACTTCGTTCAGTTACGTATTGCT",
            reverse_adapter: str = None,
            overwrite: bool = True) -> Tuple[List[str], List[str]]:

        strand_pool, strand_pool_ids = remove_adapters_from_strands(
            strands=self.strand_pool, original_strand_length=self.reference_length,
            ids=self.strand_pool_ids, starting_adapter=front_adapter)
        
        strand_pool, strand_pool_ids = remove_adapters_from_strands(
            strands=strand_pool, original_strand_length=self.reference_length,
            ids=strand_pool_ids, starting_adapter=reverse_adapter)
        
        if overwrite:
            self.strand_pool, self.strand_pool_ids = strand_pool, strand_pool_ids
            self.n_strands_pool = len(self.strand_pool)
            
        return strand_pool, strand_pool_ids

    def cluster_strand_pool(self, distance_threshold:int = None, strand_pool: List[str] = None):

        if not distance_threshold:
            distance_threshold = self.distance_threshold

        if strand_pool:
            cluster_dict = cluster_strands(
                strand_pool=strand_pool, distance_threshold=distance_threshold)
        else:
            cluster_dict = cluster_strands(
                strand_pool=self.strand_pool, distance_threshold=distance_threshold)
            strand_pool = self.strand_pool            

        self.clusters = cluster_dict['clusters']
        self.reversed_markers = cluster_dict['reversed_markers']
        self.cluster_heads = cluster_dict['cluster_heads']

        assert [len(self.clusters[i]) > len(self.clusters[i + 1]) for i in range(len(self.clusters) - 1)]
        print("Clusters are sorted")
        
        self.clustered_seqs = get_clustered_seqs(
            clusters=self.clusters, reversed_markers=self.reversed_markers, strand_pool=strand_pool)
        print("Orientation fixed in the strand pool")

        return self.clustered_seqs
    
    def generate_candidates(
            self, n_candidates: int, n_samples: int = 15, clustered_seqs: List[List[str]] = None,
            fix_orientation=False) -> List[str]:

        if clustered_seqs:  
            clustered_seqs = clustered_seqs[:n_candidates]
            self.candidates = get_candidates(clustered_seqs=clustered_seqs, n_samples=n_samples)
        else:
            clustered_seqs = self.clustered_seqs[:n_candidates]
            self.candidates = get_candidates(clustered_seqs=clustered_seqs, n_samples=n_samples)

        if fix_orientation:
            print("Fixing candidate orientations")
            self.candidates = self.fix_candidate_orientations()
            
        return self.candidates
    
    def fix_candidate_orientations(self, candidates: List[str] = None) -> List[str]:
        if not candidates:
            candidates = self.candidates

        n_candidates = len(candidates)
        reversed_markers = get_candidate_orientation(
            original_strands=self.original_strands, candidates=candidates)
        n_reversed = sum(reversed_markers)
        print(f"{n_reversed/n_candidates} candidates are reversed")
        candidates = [reverse_complement(
            candidates[ind]) if reversed_markers[ind] else candidates[ind] for ind in range(
                len(candidates))]
        return candidates
    
    def evaluate_candidates(
            self, candidates: List[str] = None, hist: bool = False, metric: str = "identity") -> Dict[str, np.ndarray]:
        
        if candidates:
            self.evaluation_dict = evaluate_candidates(
                original_strands=self.original_strands,
                candidates=candidates,
                metric=metric
            )
        else:
            self.evaluation_dict = evaluate_candidates(
                original_strands=self.original_strands,
                candidates=self.candidates,
                metric=metric
            )

        self.reference_recoveries = self.evaluation_dict['reference_recoveries']
        self.reference_strand_indices = self.evaluation_dict['reference_strand_indices']
        self.recovery_rates = self.evaluation_dict['recovery_rates']

        if hist:
            plt.hist(self.reference_recoveries)
            plt.xlabel("Recovery rate")
            plt.ylabel("Number of reference strands")
            plt.title(f"Recovery of reference strands for candidates using {metric}")
            plt.show()

        return self.evaluation_dict
    
    def fsm(self, candidates: List[str] = None, hard = False) -> bool:
        
        if hard:
            assert len(candidates) == self.n_reference_strands
        
        found = 0
        for i in self.original_strands:
            if i in candidates:
                found += 1
                continue
        
        if found == self.n_reference_strands:
            print("Found all")
            return True

        if hard:
            return False

        else:
            print(f"Found {found}")
            return False
        
    def run_pipeline(
            self, fsm: bool = False, fix_orientation: bool = False,
            remove_adapters: bool = False):

        #print("Filtering strands by length")
        #self.filter_by_length(ids=True)

        if remove_adapters:
            print("Removing adapters")
            self.remove_adapters(front_adapter=self.front_adapter, reverse_adapter=self.reverse_adapter)

        print("Clustering strands")
        self.cluster_strand_pool()

        print(f"Generating {self.n_reference_strands + 100} candidates")
        self.generate_candidates(n_candidates=self.n_reference_strands + 100, fix_orientation=fix_orientation)

        if self.original_strands:
            print("Evaluating candidates")
            self.evaluate_candidates(candidates=self.candidates, hist=True)

        if fsm:
            print(
                "Check back to see if this feature is implemented!")