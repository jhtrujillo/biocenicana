#!/usr/bin/env python3
"""compute_block_ranges.py
Generate a BED‑like table with the genomic span of every collinear block.
"""
import os, sys, argparse, pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Compute block ranges")
    parser.add_argument("--report", default="genomica_comparativa/r570/reporte_comparativo.tsv", help="Path to comparative report")
    parser.add_argument("--output", default="data/block_ranges.tsv", help="Path to output TSV")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
            'Gene2_ID','Chr2','Start2','End2','Strand2']
    rep = pd.read_csv(args.report, sep='\t', header=0, usecols=range(len(cols)), names=cols)

    # Helper to decide orientation (majority of strand concordance)
    def orientation(group):
        same = (group['Strand1'] == group['Strand2']).sum()
        return 'direct' if same >= len(group)/2 else 'inverted'

    records = []
    for block_id, sub in rep.groupby('Block_ID'):
        chr1 = sub['Chr1'].mode()[0]
        chr2 = sub['Chr2'].mode()[0]
        start1 = sub['Start1'].min()
        end1   = sub['End1'].max()
        start2 = sub['Start2'].min()
        end2   = sub['End2'].max()
        size1  = end1 - start1
        size2  = end2 - start2
        orient = orientation(sub)
        records.append([block_id, chr1, start1, end1, size1, chr2, start2, end2, size2, orient])

    out_df = pd.DataFrame(records,
        columns=['Block_ID','Chr1','Start1','End1','Size1','Chr2','Start2','End2','Size2','Orientation'])
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f'Wrote {len(out_df)} block ranges to {args.output}')

if __name__ == '__main__':
    main()

