#!/usr/bin/env python3
"""quantify_inversions.py
Quantify inversion sizes (bp) and gene count from synteny block ranges.
Outputs detailed list, global stats, and chromosome-level statistics.
"""
import os
import sys
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Quantify inversion sizes and gene counts from block ranges.")
    parser.add_argument("--ranges", required=True, help="Path to block_ranges.tsv")
    parser.add_argument("--orient", help="Path to block_orientation.tsv (optional, will try to auto-detect)")
    parser.add_argument("--out-dir", required=True, help="Output directory to save results")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Load block ranges
    if not os.path.exists(args.ranges):
        print(f"Error: Ranges file not found at {args.ranges}")
        sys.exit(1)
        
    df = pd.read_csv(args.ranges, sep='\t')
    
    # Ensure Size1 and Size2 columns exist
    if 'Size1' not in df.columns:
        df['Size1'] = df['End1'] - df['Start1']
    if 'Size2' not in df.columns:
        df['Size2'] = df['End2'] - df['Start2']

    # 2. Try to load block orientation for gene counts
    orient_path = args.orient
    if not orient_path:
        # Auto-detect block_orientation.tsv in same directory
        dir_name = os.path.dirname(args.ranges)
        possible_orient = os.path.join(dir_name, "block_orientation.tsv")
        if os.path.exists(possible_orient):
            orient_path = possible_orient

    orient_df = None
    if orient_path and os.path.exists(orient_path):
        print(f"Loading block orientation from {orient_path}")
        orient_df = pd.read_csv(orient_path, sep='\t')
        if 'NumGenes' in orient_df.columns:
            df = df.merge(orient_df[['Block_ID', 'NumGenes']], on='Block_ID', how='left')
    else:
        print("Warning: block_orientation.tsv not found or loaded. Gene counts will not be available in reports.")

    # 3. Filter for inverted blocks
    # Support both 'inverted' and 'minus' as orientation values
    inv = df[df['Orientation'].str.lower().isin(['inverted', 'minus'])].copy()
    print(f"Detected {len(inv)} inverted blocks out of {len(df)} total blocks.")

    # 4. Generate outputs
    
    # Output A: inversions_detailed.tsv
    out_detailed = os.path.join(args.out_dir, "inversions_detailed.tsv")
    detailed_cols = ['Block_ID', 'Chr1', 'Start1', 'End1', 'Size1', 'Chr2', 'Start2', 'End2', 'Size2']
    if 'NumGenes' in df.columns:
        detailed_cols.append('NumGenes')
    detailed_cols.append('Orientation')
    
    inv_detailed = inv[detailed_cols].sort_values(by=['Chr1', 'Start1'])
    inv_detailed.to_csv(out_detailed, sep='\t', index=False)
    print(f"Detailed inversions list written to {out_detailed}")

    # Output B: inversions_global_stats.tsv
    out_global = os.path.join(args.out_dir, "inversions_global_stats.tsv")
    
    if len(inv) > 0:
        stats = {
            'Metric': ['Count', 'Total_bp', 'Mean_bp', 'Median_bp', 'Min_bp', 'Max_bp'],
            'Genome1_Size': [
                len(inv),
                inv['Size1'].sum(),
                inv['Size1'].mean(),
                inv['Size1'].median(),
                inv['Size1'].min(),
                inv['Size1'].max()
            ],
            'Genome2_Size': [
                len(inv),
                inv['Size2'].sum(),
                inv['Size2'].mean(),
                inv['Size2'].median(),
                inv['Size2'].min(),
                inv['Size2'].max()
            ]
        }
        if 'NumGenes' in df.columns:
            stats['Gene_Count'] = [
                len(inv),
                inv['NumGenes'].sum(),
                inv['NumGenes'].mean(),
                inv['NumGenes'].median(),
                inv['NumGenes'].min(),
                inv['NumGenes'].max()
            ]
    else:
        stats = {
            'Metric': ['Count', 'Total_bp', 'Mean_bp', 'Median_bp', 'Min_bp', 'Max_bp'],
            'Genome1_Size': [0, 0, 0, 0, 0, 0],
            'Genome2_Size': [0, 0, 0, 0, 0, 0]
        }
        if 'NumGenes' in df.columns:
            stats['Gene_Count'] = [0, 0, 0, 0, 0, 0]
            
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(out_global, sep='\t', index=False)
    print(f"Global statistics written to {out_global}")

    # Output C: inversions_by_chromosome.tsv (Grouped by Chr1)
    out_chr = os.path.join(args.out_dir, "inversions_by_chromosome.tsv")
    if len(inv) > 0:
        chr_grouped = inv.groupby('Chr1').agg(
            Num_Inversions=('Block_ID', 'count'),
            Total_Size1_bp=('Size1', 'sum'),
            Mean_Size1_bp=('Size1', 'mean'),
            Median_Size1_bp=('Size1', 'median'),
            Total_Size2_bp=('Size2', 'sum'),
            Mean_Size2_bp=('Size2', 'mean')
        ).reset_index().sort_values(by='Num_Inversions', ascending=False)
        
        # If gene count is present, add it
        if 'NumGenes' in df.columns:
            chr_genes = inv.groupby('Chr1').agg(
                Total_Genes=('NumGenes', 'sum'),
                Mean_Genes=('NumGenes', 'mean')
            ).reset_index()
            chr_grouped = chr_grouped.merge(chr_genes, on='Chr1', how='left')
    else:
        chr_grouped = pd.DataFrame(columns=['Chr1', 'Num_Inversions', 'Total_Size1_bp', 'Mean_Size1_bp', 'Total_Size2_bp', 'Mean_Size2_bp'])
        
    chr_grouped.to_csv(out_chr, sep='\t', index=False)
    print(f"Chromosome-level statistics written to {out_chr}")

if __name__ == '__main__':
    main()
