# DNA Storage Toolkit

A lightweight collection of algorithms and utilities for experimenting with DNA-based data storage pipelines.  
Includes modules for clustering, multiple alignment, JPEG-based encoding, and a simple synthesis noise model.

---

## Features

### Clustering
Tools for grouping and denoising DNA sequence reads.
- Edit-distance hierarchical clustering  
- Multiple sequence alignment (MSA) for consensus recovery

### JPEG Encoding
A pure scan-based JPEG → DNA encoder, designed to be deterministic and testable for storage and reconstruction workflows.

### Naive Synthesis Model
A simple coin-flip per-base synthesis noise simulator:
- Random base substitutions  
- Useful for testing robustness of clustering, decoding, and consensus operations

---

## Repository Structure
dna_storage/
├── clustering/ # Edit-distance clustering + MSA tools
├── jpeg_encoding/ # Scan-based JPEG→DNA encoder
├── synthesis/ # Naive coin-flip synthesis model
└── utils/ # General helpers and IO


## Getting Started

```bash
git clone <your-repo-url>
cd dna_storage
pip install -r requirements.txt
```

## Notes
This repository is research-oriented and intended for experimentation.
Components are modular and easy to extend for custom DNA storage pipelines.
Email for issues: parvagrw02@gmail.com
