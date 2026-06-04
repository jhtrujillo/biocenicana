#!/usr/bin/env python3
"""intersect_sv_blocks.py
Simple intersection of SV regions with block orientation data.
"""
import os, sys, argparse, pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Intersect SV regions with block orientation data")
    parser.add_argument("--sv", default="genomica_comparativa/r570/tables/sv_regions.bed", help="Path to SV regions BED/TSV")
    parser.add_argument("--orient", default="data/block_orientation.tsv", help="Path to block orientation TSV")
    parser.add_argument("--output", default="genomica_comparativa/r570/tables/sv_block_overlap.tsv", help="Path to output TSV")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Load SVs
    sv_df = pd.read_csv(args.sv, sep='\t', header=0)
    # Load block orientation (contains Chr1, Chr2 fields)
    block_df = pd.read_csv(args.orient, sep='\t', header=0)

    # Simplistic join on chromosome matching (Chr1 or Chr2)
    merged = []
    for _, sv in sv_df.iterrows():
        matches = block_df[(block_df['Chr1'] == sv['CHROM']) | (block_df['Chr2'] == sv['CHROM'])]
        for _, blk in matches.iterrows():
            merged.append({
                'Block_ID': blk['Block_ID'],
                'SVTYPE': sv['SVTYPE'],
                'Chrom': sv['CHROM'],
                'SV_Start': sv['START'],
                'SV_End': sv['END'],
                'Orientation': blk['Orientation'],
                'NumGenes': blk['NumGenes'],
                'OrphanFlag': 'NA'  # placeholder – real implementation would check PAV list
            })
    out_df = pd.DataFrame(merged)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"SV‑block overlap written to {args.output}")

if __name__ == '__main__':
    main()

