# guide_location

A pipeline for mapping CRISPR guide RNA sequences to genomic locations and visualizing their positions across viral gene features.

## Overview

Given a set of guide RNA sequences (e.g., pan-sarbecovirus PspCas13b guides) and a GenBank accession number, this pipeline:

1. Downloads the reference genome (FASTA + GenBank annotation) from NCBI
2. Extracts subgenomic feature sequences (S, N, E, M, ORFs, etc.)
3. Builds a BLAST database from the reference
4. BLASTs guide sequences against the database to determine genomic coordinates
5. Outputs results in BED format and generates an interactive HTML visualization of guide positions across each gene

## Scripts

- `scripts/prepare_guide_input.py` — Converts an Excel guide file (with `ID` and `Spacer` columns) to a tab-separated CSV suitable for pipeline input
- `scripts/get_guide_loaction.py` — Main pipeline script: fetches genome data, builds BLAST DB, runs BLAST, and produces the visualization

## Usage

```bash
# Prepare guide input from Excel
python scripts/prepare_guide_input.py guides.xlsx

# Run the full pipeline
python scripts/get_guide_loaction.py \
  --accession MT007544.1 \
  --output_dir data/raw/MT007544.1 \
  --guide_input data/pan-sarbecov_pspcas13b_guides.csv
```

## Dependencies

- Python 3
- Biopython
- pandas, numpy
- plotly
- BLAST+ (`makeblastdb`, `blastn`)
- NCBI Entrez (via Biopython)
