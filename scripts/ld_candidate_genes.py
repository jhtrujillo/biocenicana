#!/usr/bin/env python3
"""
Phase 4: LD analysis and haplotype definition around candidate genes.
For each candidate gene, calculates r² between all markers within ±window_kb
and identifies haplotype blocks and tag SNPs.
"""

import argparse
import csv
import math
from collections import defaultdict

CANDIDATE_GENES = [
    ("CC01t039050.1", "chr01", 55618891, 55625629, "Sucrose synthase 2 (SuSy2)",       "Purificadora"),
    ("CC01t042770.1", "chr01", 60506673, 60519537, "Sucrose synthase 4",               "Presentacion"),
    ("CC01t044570.1", "chr01", 62833534, 62834770, "Alkaline/neutral invertase",       "Purificadora"),
    ("CC03t157420.1", "chr03", 24241458, 24246070, "Neutral/alkaline invertase mitoc", "Presentacion"),
    ("CC03t185480.1", "chr03", 62738478, 62743743, "Sucrose-phosphate synthase (SPS)", "Neutral"),
    ("CC03t185500.1", "chr03", 62761410, 62766866, "putative SPS",                     "Presentacion"),
    ("CC04t207510.1", "chr04", 25129201, 25134796, "Neutral/alkaline invertase clp",   "Purificadora"),
    ("CC04t208410.1", "chr04", 26667005, 26670298, "Cytosolic invertase 1",            "Purificadora"),
    ("CC04t230890.1", "chr04", 55774345, 55777956, "Sucrose synthase (SuSy)",          "Neutral"),
    ("CC04t196460.1", "chr04",  6117451,  6122149, "putative Alkaline/neutral inv",    "Presentacion"),
    ("CC04t199390.1", "chr04", 10079686, 10086951, "Sucrose-phosphate synthase",       "Presentacion"),
    ("CC04t217770.1", "chr04", 39342372, 39345580, "Cytosolic invertase 1",            "Presentacion"),
    ("CC05t236990.1", "chr05",  8919861,  8922720, "Alkaline/neutral invertase",       "Purificadora"),
    ("CC09t347620.1", "chr09", 31997030, 32014463, "putative SPS4",                    "Neutral"),
    ("CC10t084020.1", "chr10", 27527598, 27532806, "Sucrose-phosphate synthase",       "Neutral"),
    ("CC10t092720.1", "chr10", 40579110, 40588278, "Sucrose synthase",                 "Presentacion"),
    ("CC10t092800.1", "chr10", 40651907, 40661068, "Sucrose synthase (SuSy)",          "Positiva"),
    ("CC06t277320.1", "chr06",  None,     None,    "Beta-fructofuranosidasa Ka/Ks=29", "Positiva"),
    ("CC03t152280.1", "chr03",  None,     None,    "ERD6-like transporter Ka/Ks=26",   "Positiva"),
    ("CC03t152310.1", "chr03",  None,     None,    "ERD6-like transporter Ka/Ks=20",   "Positiva"),
    ("CC04t193040.1", "chr04",  None,     None,    "Alkaline/neutral invertase",        "Positiva"),
]

def normalize_chr(c):
    return c.lower().lstrip("chr").lstrip("0") or "0"

def pearson_r2(xs, ys):
    n = len(xs)
    if n < 5: return float('nan')
    mx = sum(xs)/n; my = sum(ys)/n
    num = sum((xs[i]-mx)*(ys[i]-my) for i in range(n))
    dx  = sum((xs[i]-mx)**2 for i in range(n))
    dy  = sum((ys[i]-my)**2 for i in range(n))
    if dx == 0 or dy == 0: return float('nan')
    r = num / math.sqrt(dx*dy)
    return r*r

def load_region_markers(vcf_path, chr_norm, start, end):
    """Load markers within [start-window, end+window] for one chromosome."""
    markers = {}  # pos -> dosages list
    samples = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#'):
                samples = line.strip().split('\t')[9:]
                continue
            cols = line.strip().split('\t')
            if len(cols) < 10: continue
            if normalize_chr(cols[0]) != chr_norm: continue
            pos = int(cols[1])
            if pos < start or pos > end: continue

            fmt = cols[8].split(':')
            acn_idx = next((i for i,f in enumerate(fmt) if f=='ACN'), -1)
            ds = []
            for g in cols[9:]:
                if g.startswith('.'):
                    ds.append(float('nan'))
                    continue
                fields = g.split(':')
                if acn_idx != -1 and len(fields) > acn_idx:
                    try: ds.append(float(fields[acn_idx].split(',')[-1]))
                    except: ds.append(float('nan'))
                else:
                    ds.append(float('nan'))
            markers[pos] = ds
    return markers, samples

def compute_ld_matrix(markers_dict):
    """Compute pairwise r² for all marker pairs."""
    positions = sorted(markers_dict.keys())
    n = len(positions)
    r2_matrix = {}
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = positions[i], positions[j]
            xi = [d for d in markers_dict[pi] if not math.isnan(d)]
            xj = [d for d in markers_dict[pj] if not math.isnan(d)]
            # Use only shared non-missing
            pairs = [(markers_dict[pi][k], markers_dict[pj][k])
                     for k in range(min(len(markers_dict[pi]), len(markers_dict[pj])))
                     if not math.isnan(markers_dict[pi][k]) and not math.isnan(markers_dict[pj][k])]
            if len(pairs) < 5:
                r2_matrix[(pi, pj)] = float('nan')
                continue
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            r2_matrix[(pi, pj)] = pearson_r2(xs, ys)
    return positions, r2_matrix

def find_tag_snps(positions, r2_matrix, r2_threshold=0.8):
    """
    Greedy tag SNP selection: select minimum set of markers such that
    every marker is captured by at least one tag with r² >= threshold.
    """
    n = len(positions)
    captured = set()
    tags = []

    # Build adjacency: which markers does each marker capture?
    captures = {p: {p} for p in positions}
    for i in range(n):
        for j in range(i+1, n):
            pi, pj = positions[i], positions[j]
            r2 = r2_matrix.get((pi, pj), float('nan'))
            if not math.isnan(r2) and r2 >= r2_threshold:
                captures[pi].add(pj)
                captures[pj].add(pi)

    remaining = set(positions)
    while remaining:
        # Pick marker that captures most remaining uncaptured
        best = max(remaining, key=lambda p: len(captures[p] & remaining))
        tags.append(best)
        captured.update(captures[best])
        remaining -= captures[best]

    return tags

def run(vcf_path, out_dir, window_kb=500):
    import os; os.makedirs(out_dir, exist_ok=True)
    window_bp = window_kb * 1000

    summary_rows = []
    tag_snp_rows = []

    for gene_id, chr_name, start, end, func, cat in CANDIDATE_GENES:
        if start is None:
            print(f"[LD] {gene_id}: sin coordenadas — omitido")
            continue

        chr_norm = normalize_chr(chr_name)
        region_start = max(0, start - window_bp)
        region_end   = end + window_bp
        gene_mid = (start + end) // 2

        print(f"[LD] {gene_id} ({chr_name}:{start:,}-{end:,}) ventana ±{window_kb}kb...")
        markers, samples = load_region_markers(vcf_path, chr_norm,
                                               region_start, region_end)

        if len(markers) < 2:
            print(f"  Sin marcadores suficientes en la región ({len(markers)})")
            summary_rows.append({
                'Gene': gene_id, 'Funcion': func, 'Categoria': cat,
                'Chr': chr_name, 'Gene_start': start, 'Gene_end': end,
                'N_markers_region': len(markers), 'N_markers_en_gen': 0,
                'N_tag_snps': 0, 'Max_r2': 'NA', 'Mean_r2': 'NA',
                'Tag_SNPs': 'NA', 'Estado': 'Sin datos'
            })
            continue

        positions, r2_matrix = compute_ld_matrix(markers)

        # Markers inside the gene body
        in_gene = [p for p in positions if start <= p <= end]

        # r² stats
        r2_vals = [v for v in r2_matrix.values() if not math.isnan(v)]
        max_r2  = round(max(r2_vals), 3) if r2_vals else float('nan')
        mean_r2 = round(sum(r2_vals)/len(r2_vals), 3) if r2_vals else float('nan')

        # Tag SNPs (r² >= 0.8)
        tags = find_tag_snps(positions, r2_matrix, r2_threshold=0.8)
        tags_str = ";".join(f"{chr_name}:{p}" for p in sorted(tags))

        # High LD pairs with gene markers
        high_ld_with_gene = []
        for p_out in positions:
            if p_out in in_gene: continue
            for p_in in in_gene:
                key = (min(p_out, p_in), max(p_out, p_in))
                r2  = r2_matrix.get(key, float('nan'))
                if not math.isnan(r2) and r2 >= 0.5:
                    high_ld_with_gene.append((p_out, p_in, round(r2, 3)))

        high_ld_with_gene.sort(key=lambda x: -x[2])

        print(f"  {len(positions)} marcadores en región | {len(in_gene)} en gen | "
              f"{len(tags)} tag SNPs | max r²={max_r2}")

        summary_rows.append({
            'Gene': gene_id, 'Funcion': func, 'Categoria': cat,
            'Chr': chr_name, 'Gene_start': start, 'Gene_end': end,
            'N_markers_region': len(positions),
            'N_markers_en_gen': len(in_gene),
            'N_tag_snps': len(tags),
            'Max_r2': max_r2, 'Mean_r2': mean_r2,
            'Tag_SNPs': tags_str,
            'Estado': 'OK'
        })

        # Tag SNP detail rows
        for tag_pos in sorted(tags):
            dist_to_gene = min(abs(tag_pos - start), abs(tag_pos - end))
            in_gen = "Sí" if start <= tag_pos <= end else "No"
            # Find best r² of this tag with a gene-body marker
            best_r2_with_gene = 'NA'
            if in_gene:
                r2s = []
                for p_in in in_gene:
                    key = (min(tag_pos, p_in), max(tag_pos, p_in))
                    r2  = r2_matrix.get(key, float('nan'))
                    if not math.isnan(r2): r2s.append(r2)
                if r2s:
                    best_r2_with_gene = round(max(r2s), 3)
            tag_snp_rows.append({
                'Gene': gene_id, 'Funcion': func, 'Categoria': cat,
                'Chr': chr_name, 'Tag_pos': tag_pos,
                'En_gen': in_gen,
                'Dist_gen_bp': dist_to_gene if in_gen == "No" else 0,
                'r2_con_gen': best_r2_with_gene,
                'Marcador_ID': f"{normalize_chr(chr_name)}_{tag_pos}"
            })

    # Write outputs
    print("\n[LD] Guardando resultados...")
    fields_sum = ['Gene','Funcion','Categoria','Chr','Gene_start','Gene_end',
                  'N_markers_region','N_markers_en_gen','N_tag_snps',
                  'Max_r2','Mean_r2','Tag_SNPs','Estado']
    with open(f"{out_dir}/ld_summary.tsv", 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields_sum, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(summary_rows)
    print(f"  → {out_dir}/ld_summary.tsv")

    fields_tag = ['Gene','Funcion','Categoria','Chr','Tag_pos','En_gen',
                  'Dist_gen_bp','r2_con_gen','Marcador_ID']
    with open(f"{out_dir}/tag_snps.tsv", 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields_tag, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(tag_snp_rows)
    print(f"  → {out_dir}/tag_snps.tsv")

    # Print summary table
    print("\n" + "="*95)
    print(f"{'Gen':<20} {'Chr':<7} {'Mks_región':>10} {'Mks_gen':>8} {'Tags':>5} {'max_r²':>7} {'Categoría'}")
    print("="*95)
    for r in summary_rows:
        print(f"{r['Gene']:<20} {str(r['Chr']):<7} {str(r['N_markers_region']):>10} "
              f"{str(r['N_markers_en_gen']):>8} {str(r['N_tag_snps']):>5} "
              f"{str(r['Max_r2']):>7}  {r['Categoria']}")
    print("="*95)
    print(f"\nTotal tag SNPs identificados: {len(tag_snp_rows)}")
    print(f"Candidatos para marcadores KASP: genes con marcadores dentro del gen (En_gen=Sí)")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Phase 4: LD y tag SNPs alrededor de genes candidatos.")
    p.add_argument("--vcf",     required=True, help="VCF de la población")
    p.add_argument("--out-dir", default="genomica_comparativa/mapa_genetico/fase4_ld")
    p.add_argument("--window-kb", type=int, default=500,
                   help="Ventana de búsqueda en kb a cada lado del gen (default: 500)")
    args = p.parse_args()
    run(args.vcf, args.out_dir, args.window_kb)
