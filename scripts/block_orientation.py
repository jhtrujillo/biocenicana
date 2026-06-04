#!/usr/bin/env python3
"""block_orientation.py
Compute block orientation (direct vs inverted) from reporte_comparativo.tsv.
"""
import os, sys, argparse, pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Compute block orientation (direct vs inverted)")
    parser.add_argument("--report", default="genomica_comparativa/r570/reporte_comparativo.tsv", help="Path to comparative report")
    parser.add_argument("--output", default="data/block_orientation.tsv", help="Path to output TSV")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # Load report (tab‑separated, first line is header)
    cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1','Gene2_ID','Chr2','Start2','End2','Strand2']
    df = pd.read_csv(args.report, sep='\t', usecols=range(len(cols)), names=cols, header=0)

    # Group by block
    summary = []
    for block, sub in df.groupby('Block_ID'):
        # majority orientation
        same = (sub['Strand1'] == sub['Strand2']).sum()
        orientation = 'direct' if same >= len(sub)/2 else 'inverted'
        # representative chromosomes (most common)
        chr1 = sub['Chr1'].mode()[0]
        chr2 = sub['Chr2'].mode()[0]
        summary.append([block, orientation, len(sub), chr1, chr2])

    out_df = pd.DataFrame(summary, columns=['Block_ID','Orientation','NumGenes','Chr1','Chr2'])
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Block orientation written to {args.output}")

if __name__ == '__main__':
    main()

