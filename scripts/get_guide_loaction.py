"""
This script retrieves the genomic location of a guide RNA sequence by performing a BLAST search against a reference database. 
It takes a guide sequence in FASTA format, queries it against a BLAST database and visualize the result.
Author: Rebecca Bengtsson
Date: 2026-07-15
"""

import argparse
import subprocess
import os, sys
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from pathlib import Path
from pyexpat import features
from tracemalloc import start

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))
import src.entrez as entrez
import src.utils as utils
import src.visualization as visualization



def parse_args():
    parser = argparse.ArgumentParser(description="Extract subgenomic sequences from given genome")
    parser.add_argument("-a", "--accession", help="GenBank accession number", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save extracted sequences", type=Path, required=True)
    parser.add_argument("-i", "--guide_input", help="Input file containing guide sequences", type=Path, required=True)
    return parser.parse_args()



def fetch_files(accession, gbk_file_path, fasta_file_path):
    record = entrez.fetch_gbk(accession)
    SeqIO.write(record, gbk_file_path, "genbank")
    record = entrez.fetch_fasta(accession)
    SeqIO.write(record, fasta_file_path, "fasta")



def _build_reference_db(accession, output_dir, gbk_file_path, fasta_file_path):
    """
    Build a BLAST database from the reference genome.
    """
    
    # #### Step 1 - Download reference genome from NCBI
    fetch_files(accession, gbk_file_path, fasta_file_path)
    print(f"Fetching GenBank and FASTA files for accession {accession} completed.")

    feature_locations = entrez.get_feature_locations(gbk_file_path)

    #### Step 2 - Extract feature coordinates and save to file
    feature_list = [feature for feature in feature_locations.keys()]
    for feature in feature_list:
        start = feature_locations[feature]['start']
        end = feature_locations[feature]['end']
        entrez.fetch_gene_sequence(feature, start, end, fasta_file_path, output_dir)
    
    #### Step 3 - Build reference BLAST database
    db_path = utils.make_blast_db(fasta_file_path)
    print(f"BLAST database created at: {db_path}")

    return db_path



def _blast_query(guide_input, db_path, output_dir):
    """
    Perform a BLAST search of the guide sequence against the reference database.
    """
    
    #### Read input file results
    df = pd.read_csv(guide_input, sep="\t", header=0)
    df['target sequence'] = df['guide sequence'].apply(utils.reverse_complement)

    records = [
        SeqRecord(Seq(seq), id=str(seq_id), description="")
        for seq_id, seq in zip(df['guide_id'], df["target sequence"])
    ]

    query_guides_file = output_dir / "query_guides.fasta"
    
    #### Write query guides to FASTA file
    SeqIO.write(records, query_guides_file, "fasta")

    #### Run BLAST search
    blastn_output_file = output_dir / "blastn_output.txt"
    utils.blastn_query(query_guides_file, db_path, blastn_output_file)

    return blastn_output_file



def _generate_bed_file(blastn_output_file, output_dir):

    blastn_df =utils.format_blastn_results(blastn_output_file)
    
    # Write BED file
    df_bed = blastn_df[['sseqid', 'sstart', 'send', 'qseqid']]
    bed_file_path = output_dir / "blastn_output.bed"
    df_bed.to_csv(bed_file_path, sep="\t", header=False, index=False)
    print(f"BED file created at: {bed_file_path}")

    return bed_file_path


def main():

    args = parse_args()
    accession = args.accession
    guide_input = args.guide_input

    output_dir = args.output_dir / accession
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    gbk_file_path = output_dir / f"{accession}.gbk"
    fasta_file_path = output_dir / f"{accession}.fasta"

    db_path = _build_reference_db(accession, output_dir, gbk_file_path, fasta_file_path)
    blastn_output_file = _blast_query(guide_input, db_path, output_dir)

    blastn_output_file = output_dir / "blastn_output.txt"
    bed_file_path = _generate_bed_file(blastn_output_file, output_dir)

    # bed_file_path = output_dir / "blastn_output.bed"
    visualization.visualize_guide_positions(output_dir, bed_file_path, outfile_name="PansarbeCov_guide_positions.html")

    
    



if __name__ == "__main__":
    main()