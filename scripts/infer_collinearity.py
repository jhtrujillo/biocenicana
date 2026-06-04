#!/usr/bin/env python3
import sys
import os
import re
from collections import defaultdict

def parse_gff(gff_path):
    """
    Parses a GFF3 file and returns:
    1. A dictionary mapping gene_id -> (chrom, start, end, index)
    2. A dictionary mapping chrom -> list of (gene_id, start, end) sorted by position
    """
    print(f"Parsing GFF: {gff_path}...")
    gene_map = {}
    chrom_genes = defaultdict(list)
    
    if not os.path.exists(gff_path):
        print(f"Error: GFF file not found: {gff_path}")
        sys.exit(1)
        
    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in ('gene', 'mRNA', 'transcript'):
                continue
                
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            # Extract ID
            gene_id = ""
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match:
                gene_id = id_match.group(1).strip()
            else:
                continue
                
            # If it's a gene feature, we keep it as a fallback, but mRNA is preferred.
            # R570 GFF has gene ID like SoffiXsponR570.10Ag000100.v2.1 and mRNA ID like SoffiXsponR570.10Ag000100.1.
            # S. spontaneum has mRNA features.
            # We store the ID exactly.
            if gene_id not in gene_map:
                gene_map[gene_id] = (chrom, start, end)
                chrom_genes[chrom].append((gene_id, start, end))
                
    # Sort and index
    indexed_gene_map = {}
    for chrom, genes in chrom_genes.items():
        genes.sort(key=lambda x: x[1]) # Sort by start coordinate
        for idx, (gene_id, start, end) in enumerate(genes):
            indexed_gene_map[gene_id] = (chrom, start, end, idx)
            
    print(f"Loaded {len(indexed_gene_map)} genes/transcripts from {gff_path}.")
    return indexed_gene_map

def parse_collinearity(collinearity_path):
    """
    Parses a collinearity file and returns a list of syntenic pairs: (gene1, gene2, evalue)
    """
    print(f"Parsing collinearity: {collinearity_path}...")
    pairs = []
    if not os.path.exists(collinearity_path):
        print(f"Error: Collinearity file not found: {collinearity_path}")
        sys.exit(1)
        
    with open(collinearity_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    data = parts[1].strip().split()
                    if len(data) >= 3:
                        g1 = data[0]
                        g2 = data[1]
                        try:
                            ev = float(data[2])
                        except ValueError:
                            ev = 0.0
                        pairs.append((g1, g2, ev))
    print(f"Loaded {len(pairs)} syntenic pairs from {collinearity_path}.")
    return pairs

def clean_1940_id(x):
    # CC-01-1940 IDs are identical in both files, but let's clean them to match them if there's minor isoform difference (e.g. .1 vs .2)
    # Actually, they are CC01t...1 in both, so let's keep them exact.
    return x

def main():
    if len(sys.argv) < 6:
        print("Usage: python3 infer_collinearity.py <gff_r570> <gff_spont> <coll_1940_vs_r570> <coll_1940_vs_spont> <output_collinearity>")
        sys.exit(1)
        
    gff_r570_path = sys.argv[1]
    gff_spont_path = sys.argv[2]
    coll_r570_path = sys.argv[3]
    coll_spont_path = sys.argv[4]
    out_path = sys.argv[5]
    
    # 1. Parse GFFs
    r570_genes = parse_gff(gff_r570_path)
    spont_genes = parse_gff(gff_spont_path)
    
    # 2. Parse collinearity files
    pairs_r570 = parse_collinearity(coll_r570_path)
    pairs_spont = parse_collinearity(coll_spont_path)
    
    # 3. Build mapping: gene_1940 -> list of r570 genes
    map_1940_to_r570 = defaultdict(list)
    for g_1940, g_r570, ev in pairs_r570:
        is_1940_g1 = g_1940.startswith("CC01")
        g_1940_actual = g_1940 if is_1940_g1 else g_r570
        g_r570_actual = g_r570 if is_1940_g1 else g_1940
        map_1940_to_r570[clean_1940_id(g_1940_actual)].append((g_r570_actual, ev))
        
    # Build mapping: gene_1940 -> list of spont genes
    map_1940_to_spont = defaultdict(list)
    for g_1940, g_spont, ev in pairs_spont:
        is_1940_g1 = g_1940.startswith("CC01")
        g_1940_actual = g_1940 if is_1940_g1 else g_spont
        g_spont_actual = g_spont if is_1940_g1 else g_1940
        map_1940_to_spont[clean_1940_id(g_1940_actual)].append((g_spont_actual, ev))
        
    # 4. Infer R570 vs Spontaneum pairs
    inferred_pairs = []
    seen = set()
    for g_1940, r570_list in map_1940_to_r570.items():
        if g_1940 in map_1940_to_spont:
            for g_r570, ev_r570 in r570_list:
                for g_spont, ev_spont in map_1940_to_spont[g_1940]:
                    pair_key = (g_r570, g_spont)
                    if pair_key not in seen:
                        seen.add(pair_key)
                        inferred_pairs.append((g_r570, g_spont, max(ev_r570, ev_spont)))
                        
    print(f"Inferred {len(inferred_pairs)} raw syntenic pairs.")
    
    # 5. Group into collinear blocks
    valid_pairs = []
    for g_r570, g_spont, ev in inferred_pairs:
        if g_r570 in r570_genes and g_spont in spont_genes:
            valid_pairs.append((g_r570, g_spont, ev))
            
    print(f"Valid pairs with coordinates in GFFs: {len(valid_pairs)}")
    
    # Group pairs by (chrom_r570, chrom_spont)
    grouped_pairs = defaultdict(list)
    for g_r570, g_spont, ev in valid_pairs:
        chr_r570, _, _, idx_r570 = r570_genes[g_r570]
        chr_spont, _, _, idx_spont = spont_genes[g_spont]
        grouped_pairs[(chr_r570, chr_spont)].append((g_r570, idx_r570, g_spont, idx_spont, ev))
        
    # For each chromosome combination, cluster into collinear blocks
    max_gap = 25
    min_block_size = 5
    blocks = []
    
    for (chr_r, chr_s), pairs in grouped_pairs.items():
        pairs.sort(key=lambda x: x[1])
        for orientation in ('plus', 'minus'):
            current_block = []
            for p in pairs:
                g_r, idx_r, g_s, idx_s, ev = p
                if not current_block:
                    current_block.append(p)
                    continue
                    
                _, last_idx_r, _, last_idx_s, _ = current_block[-1]
                gap_r = idx_r - last_idx_r
                gap_s = (idx_s - last_idx_s) if orientation == 'plus' else (last_idx_s - idx_s)
                
                if 0 < gap_r <= max_gap and 0 < gap_s <= max_gap:
                    current_block.append(p)
                elif gap_r > max_gap or gap_s > max_gap or gap_s <= 0:
                    if len(current_block) >= min_block_size:
                        blocks.append((chr_r, chr_s, orientation, current_block))
                    current_block = [p]
                    
            if len(current_block) >= min_block_size:
                blocks.append((chr_r, chr_s, orientation, current_block))
                
    # 6. Write output in MCScanX collinearity format
    print(f"Clustered into {len(blocks)} collinearity blocks.")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("############### Parameters ###############\n")
        f.write("# MATCH_SCORE: 50\n")
        f.write("# MATCH_SIZE: 5\n")
        f.write("# GAP_PENALTY: -1\n")
        f.write("# OVERLAP_WINDOW: 5\n")
        f.write("# E_VALUE: 1e-05\n")
        f.write("# MAX GAPS: 25\n")
        f.write("############### Statistics ###############\n")
        total_collinear_genes = sum(len(b[3]) * 2 for b in blocks)
        f.write(f"# Number of collinear genes: {total_collinear_genes}\n")
        f.write("##########################################\n")
        
        for block_idx, (chr_r, chr_s, orient, block_pairs) in enumerate(blocks):
            score = len(block_pairs) * 50.0
            evalue = min(p[4] for p in block_pairs)
            f.write(f"## Alignment {block_idx}: score={score:.1f} e_value={evalue} N={len(block_pairs)} {chr_r}&{chr_s} {orient}\n")
            for pair_idx, (g_r, idx_r, g_s, idx_s, ev) in enumerate(block_pairs):
                f.write(f"  {block_idx}-{pair_idx}:\t{g_r}\t{g_s}\t{ev}\n")
                
    print(f"Successfully generated collinearity file: {out_path}")

if __name__ == "__main__":
    main()
