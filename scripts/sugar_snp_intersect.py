#!/usr/bin/env python3
"""sugar_snp_intersect.py
Intersect VCF SNPs with sugar-gene genomic ranges.
Output: genomica_comparativa/r570/tables/sugar_snp_overlap.tsv
Columns: Gene_ID, Chr, Start, End, Block_ID, Status, Num_SNPs, SNP_Positions
"""
import os, csv, collections
import pandas as pd

import os, csv, collections, argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Intersect VCF SNPs with sugar-gene genomic ranges")
    parser.add_argument("--vcf", default="benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf", help="VCF file path")
    parser.add_argument("--sugar-ids", default="data/sugar_gene_ids.txt", help="Sugar gene IDs file path")
    parser.add_argument("--report", default="genomica_comparativa/r570/reporte_comparativo.tsv", help="Comparative report TSV")
    parser.add_argument("--output", default="genomica_comparativa/r570/tables/sugar_snp_overlap.tsv", help="Output TSV path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)


    # ── Load sugar gene IDs ────────────────────────────────────────────────────
    with open(args.sugar_ids) as f:
        sugar_ids = set(l.strip() for l in f if l.strip())
    print(f"Loaded {len(sugar_ids)} sugar gene IDs")

    # ── Load report ────────────────────────────────────────────────────────────
    cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
            'Gene2_ID','Chr2','Start2','End2','Strand2']
    rep = pd.read_csv(args.report, sep='\t', header=0,
                      usecols=range(len(cols)), names=cols, low_memory=False)
    # Filter to sugar genes present in report
    sugar_rep = rep[rep['Gene1_ID'].isin(sugar_ids)].copy()
    print(f"Sugar genes found in report: {len(sugar_rep)}")

    # ── Normalise chromosome names (chr01 → 1, chr10 → 10) ────────────────────
    def norm_chr(name):
        s = str(name).lower().replace('chr', '').lstrip('0')
        return s if s else '0'

    sugar_rep['Chr1_norm'] = sugar_rep['Chr1'].apply(norm_chr)

    # ── Build SNP index: {norm_chrom: sorted list of positions} ───────────────
    snp_index = collections.defaultdict(list)
    with open(args.vcf) as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            parts = line.split('\t', 8)
            if len(parts) < 2:
                continue
            chrom = norm_chr(parts[0])
            pos   = int(parts[1])
            snp_index[chrom].append(pos)

    # Sort each list for binary search
    for chrom in snp_index:
        snp_index[chrom].sort()
    print(f"Loaded SNPs for {len(snp_index)} chromosomes")

    # ── Intersect ──────────────────────────────────────────────────────────────
    import bisect
    rows = []
    for _, row in sugar_rep.iterrows():
        chrom   = row['Chr1_norm']
        start   = int(row['Start1'])
        end     = int(row['End1'])
        snps    = snp_index.get(chrom, [])
        lo      = bisect.bisect_left(snps, start)
        hi      = bisect.bisect_right(snps, end)
        hits    = snps[lo:hi]
        rows.append({
            'Gene_ID':       row['Gene1_ID'],
            'Chr':           row['Chr1'],
            'Start':         start,
            'End':           end,
            'Block_ID':      row['Block_ID'],
            'Status':        row['Status'],
            'Num_SNPs':      len(hits),
            'SNP_Positions': ';'.join(map(str, hits[:20]))  # cap at 20
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"Written {len(out_df)} sugar genes to {args.output}")
    print(f"  → {(out_df['Num_SNPs'] > 0).sum()} genes have ≥1 overlapping SNP")

if __name__ == '__main__':
    main()

