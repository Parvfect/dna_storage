

import numpy as np
from collections import Counter
from typing import List


def handle_long_run_dna(traces, pointers, current_base, q):
    """
    Handles a long run of a DNA base during alignment.
    """
    observed_lengths = []

    for trace, idx in zip(traces, pointers):
        run_length = 0
        while idx + run_length < len(trace) and trace[idx + run_length] == current_base:
            run_length += 1
        observed_lengths.append(run_length)
    
    observed_lengths.sort()
    median_length = observed_lengths[len(observed_lengths) // 2]
    estimated_true_length = int(round(median_length / (1 - q)))

    updated_pointers = [idx + length for idx, length in zip(pointers, observed_lengths)]
    estimated_run = current_base * estimated_true_length

    return estimated_run, updated_pointers


def majority_base_at_pointers(traces, pointers):
    """
    Returns the majority base among all current pointer positions.
    """
    current_bases = [
        trace[idx] for trace, idx in zip(traces, pointers) if idx < len(trace)
    ]
    if not current_bases:
        return None, 0
    counter = Counter(current_bases)
    majority_base, count = counter.most_common(1)[0]
    return majority_base, count


def reconstruct_dna_strand(traces: List[str], q: float, n: int):
    """
    Reconstructs a DNA strand from multiple deletion-corrupted traces.
    
    Args:
        traces - List of deletion profiles
        q - Deletion rate
        n - Original reference length
    """
    reconstructed = []
    pointers = [0] * len(traces)
    sqrt_n = int(np.sqrt(n))

    while len(reconstructed) < n:
        majority_base, count = majority_base_at_pointers(traces, pointers)
        if majority_base is None:
            break  # All pointers have reached the end

        # Check if it's a long run
        run_lengths = []
        for trace, idx in zip(traces, pointers):
            run_length = 0
            while idx + run_length < len(trace) and trace[idx + run_length] == majority_base:
                run_length += 1
            run_lengths.append(run_length)

        long_run_detected = sum(length >= sqrt_n for length in run_lengths) >= len(traces) // 2

        if long_run_detected:
            # Use long run handler
            estimated_run, pointers = handle_long_run_dna(traces, pointers, majority_base, q)
            reconstructed.append(estimated_run)
        else:
            # Use majority vote for one base
            reconstructed.append(majority_base)
            for i in range(len(traces)):
                if pointers[i] < len(traces[i]) and traces[i][pointers[i]] == majority_base:
                    pointers[i] += 1

    return ''.join(reconstructed[:n])

def get_top_n_percent_reads(synth_strands, n=90):

    lens = [len(i) for i in synth_strands]
    percentile = np.percentile(lens, n)

    return [i for i in synth_strands if len(i) >= percentile]

def bidirectional_alignment(traces, q, n, filter_percentile=90):
    """
    Reconstructs a DNA strand from multiple deletion-corrupted traces.
    Uses forward and reverse prediction to obtain final prediction. 
    Filters large deletion rates out.
    
    Args:
        traces - List of deletion profiles
        q - Deletion rate
        n - Original reference length
    """

    traces = get_top_n_percent_reads(traces, filter_percentile)

    forward_strand = reconstruct_dna_strand(traces, q, n)
    reverse_strand = reconstruct_dna_strand([i[::-1] for i in traces], q, n)[::-1]

    return forward_strand[:int(n / 2)] + reverse_strand[int(n / 2):]
    #return forward_strand