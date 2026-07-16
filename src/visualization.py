"""
This script provides functions for visualizing guide RNA positions using Plotly. 
"""

import os, sys
from zipfile import Path
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from . import entrez
from pathlib import Path
import numpy as np
import pandas as pd


def get_gene_length(gene_coords_dict):
    for gene, coordinates in gene_coords_dict.items():
        start, end = coordinates["start"], coordinates["end"]
        gene_lengths = int(end) - int(start)
        gene_coords_dict[gene] = {"start": start, "end": end, "length": gene_lengths}
    return gene_coords_dict


def visualize_guide_positions(output_dir, bed_file_path, 
                              outfile_name="guide_positions.html"):
    
    # Get gene coordinates from Entrez
    acession = Path(bed_file_path).parts[-2]

    gbk_dir_path = os.path.dirname(bed_file_path)
    gbk_file_path = os.path.join(gbk_dir_path, f"{acession}.gbk")
    gene_coords = entrez.get_feature_locations(gbk_file_path)

    gene_coords = get_gene_length(gene_coords)


    # Create a subplot for each unique gene
    df_bed = pd.read_csv(bed_file_path, sep="\t", header=None, names=["gene", "start", "end", "guide_id"])
    genes = df_bed['gene'].unique()
    n = len(genes)

    def assign_tracks(df, start_col='start', end_col='end'):
        """
        Assign a track (y-level) to each interval so overlapping guides
        are stacked, non-overlapping guides share a track.
        """
        df = df.sort_values(start_col).reset_index(drop=True)
        track_ends = []  # tracks[i] = end position of last interval placed on that track
        
        tracks = []
        for _, row in df.iterrows():
            placed = False
            for t, end in enumerate(track_ends):
                if row[start_col] >= end:  # no overlap with last item on this track
                    track_ends[t] = row[end_col]
                    tracks.append(t)
                    placed = True
                    break
            if not placed:
                track_ends.append(row[end_col])
                tracks.append(len(track_ends) - 1)

        df['track'] = tracks
        return df

    
    # Create subplots for each gene
    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[f"{gene}" for gene in genes],
        shared_xaxes=False,
        vertical_spacing=0.1
    )


    for row_idx, gene in enumerate(genes, start=1):
        gene_df = df_bed[df_bed['gene'] == gene]
        gene_df = assign_tracks(gene_df)
        gene_length = gene_coords[gene]['length']
        
        for _, guide in gene_df.iterrows():
            fig.add_trace(
                go.Scatter(
                    x=[guide['start'], guide['end']],
                    y=[guide['track'], guide['track']],
                    mode='lines',
                    line=dict(color='blue', width=4),
                    marker=dict(size=8),
                    name=f"Guide {guide['guide_id']}",
                    hoverinfo='text',
                    hovertext=f"Guide ID: {guide['guide_id']}<br>Start: {guide['start']}<br>End: {guide['end']}"
                ),
                row=row_idx, col=1
            )
    
    fig.update_layout(
        title_text="Pan-SarbeCov Guide RNA Positions",
        height=300 * n,  # Adjust height based on number of genes
        showlegend=False,
        xaxis_title="Gene Position")
    
    # Save the figure as an HTML file
    output_file_path = os.path.join(output_dir, outfile_name)
    fig.write_html(output_file_path)

        

