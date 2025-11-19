
# JPEG Encoding Module

This module provides a simple, scan-based pipeline for converting JPEG bitstreams into DNA strands.  
It is designed for experimentation with minimal assumptions and lightweight constraints handling.

---

## Overview

The encoder operates directly on the **raw JPEG bitstream**, without performing any JPEG decoding or transformation.  
The bitstream is split into fixed-size segments and converted into DNA strands using a basic binary→nucleotide mapping.

To satisfy biochemical constraints (e.g., limiting homopolymer length or GC content), each strand is XORed with a **pseudo-random seed**, ensuring improved sequence diversity and structure.

---

## Pipeline

1. **Read JPEG file**  
   The raw binary data is loaded as a bitstream.

2. **Split into segments**  
   The bitstream is chunked into fixed-length blocks corresponding to the target DNA strand size.

3. **Seed-based XOR**  
   Each block is XORed with a per-strand seed:
   - helps break long runs of identical bits  
   - improves GC balance  
   - reduces homopolymers  
   - preserves decodability (seed stored or encoded)

4. **Binary → DNA encoding**  
   Uses a simple mapping (e.g., 2-bit → 1 nucleotide) to produce the final DNA strands.

---

## Features

- Pure scan-based JPEG→DNA conversion  
- Deterministic and streaming-friendly  
- XOR seed constraint smoothing  
- Modular encoder components  
- Easy to integrate with clustering or synthesis modules

---
