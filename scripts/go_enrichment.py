#!/usr/bin/env python3
"""go_enrichment.py
GO/KEGG enrichment via g:Profiler REST API for structurally affected genes
(genes in inverted blocks + orphan genes).
"""
import os, sys, argparse, requests, pandas as pd

API_URL  = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'

def write_placeholder(path, note, n_genes):
    pd.DataFrame([{
        'source': 'NOTE',
        'term_id': 'N/A',
        'term_name': note,
        'p_value': 'N/A',
        'gene_count': n_genes,
        'intersection_size': 0,
        'genes': ''
    }]).to_csv(path, sep='\t', index=False)
    print(f"Placeholder written to {path}")

def main():
    parser = argparse.ArgumentParser(description="GO/KEGG enrichment via g:Profiler")
    parser.add_argument("--orient", default="data/block_orientation.tsv", help="Block orientation TSV")
    parser.add_argument("--report", default="genomica_comparativa/r570/reporte_comparativo.tsv", help="Comparative report TSV")
    parser.add_argument("--output", default="genomica_comparativa/r570/tables/go_enrichment.tsv", help="Output TSV path")
    parser.add_argument("--organism", default="sbicolor", help="g:Profiler organism ID")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────────
    orient_df = pd.read_csv(args.orient, sep='\t', low_memory=False)
    cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
            'Gene2_ID','Chr2','Start2','End2','Strand2']
    rep = pd.read_csv(args.report, sep='\t', header=0,
                      usecols=range(len(cols)), names=cols, low_memory=False)

    # Inverted block gene IDs
    inv_blocks = set(orient_df[orient_df['Orientation'] == 'inverted']['Block_ID'])
    inv_genes  = set(rep[rep['Block_ID'].isin(inv_blocks)]['Gene1_ID'].dropna())

    # Orphan gene IDs (cap to 500 for API)
    orphan_mask  = rep['Status'].isin(['Orphan_G1', 'Orphan_G2'])
    orphan_genes = set(rep[orphan_mask]['Gene1_ID'].dropna())

    # Combined — prioritise inverted genes, then orphans, cap at 2000
    combined = list(inv_genes)[:1000] + list(orphan_genes - inv_genes)[:1000]
    combined = combined[:2000]
    print(f"Query genes: {len(combined)} ({len(inv_genes)} inverted, {len(orphan_genes)} orphan)")

    # ── g:Profiler REST call ───────────────────────────────────────────────────
    payload = {
        'organism':       args.organism,
        'query':          combined,
        'sources':        ['GO:BP', 'GO:MF', 'KEGG'],
        'user_threshold': 0.05,
        'no_evidences':   False,
        'all_results':    False,
        'no_iea':         False,
        'domain_scope':   'annotated',
        'significance_threshold_method': 'g_SCS'
    }

    try:
        resp = requests.post(API_URL, json=payload, timeout=90)
        resp.raise_for_status()
        result   = resp.json()
        enriched = result.get('result', [])
        if enriched:
            rows = []
            for entry in enriched:
                rows.append({
                    'source':            entry.get('source', ''),
                    'term_id':           entry.get('native', ''),
                    'term_name':         entry.get('name', ''),
                    'p_value':           entry.get('p_value', ''),
                    'gene_count':        entry.get('query_size', ''),
                    'intersection_size': entry.get('intersection_size', ''),
                    'genes':             ','.join(entry.get('intersections', []))
                })
            out_df = pd.DataFrame(rows)
            out_df.to_csv(args.output, sep='\t', index=False)
            print(f"GO enrichment: {len(out_df)} terms written to {args.output}")
        else:
            write_placeholder(args.output,
                f'g:Profiler returned 0 enriched terms for {args.organism} query',
                len(combined))
    except Exception as e:
        print(f"g:Profiler API error: {e}")
        write_placeholder(args.output,
            f'g:Profiler API unavailable: {e}',
            len(combined))

if __name__ == '__main__':
    main()

