#!/usr/bin/env python3
# scripts/check_sucrose_genes.py
# Análisis evolutivo de genes asociados al metabolismo y transporte de sacarosa

import os
import csv
import re

import argparse

parser = argparse.ArgumentParser(description="Analyze sucrose metabolism/transport genes")
parser.add_argument("--gff1", default="benchmarks/genomas/1940/CC-01-1940.gff3")
parser.add_argument("--gff2", default="benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3")
parser.add_argument("--report", default="genomica_comparativa/r570/reporte_comparativo.tsv")
parser.add_argument("--kaks", default="genomica_comparativa/r570/kaks_1940_vs_r570.tsv")
parser.add_argument("--name1", default="CC 1940")
parser.add_argument("--name2", default="R570")
args = parser.parse_args()

gff1_path = args.gff1
gff2_path = args.gff2
report_path = args.report
kaks_path = args.kaks
name1 = args.name1
name2 = args.name2

def load_gff_notes(gff_path):
    notes = {}
    if not os.path.exists(gff_path):
        return notes
    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ('gene', 'mRNA'):
                continue
            
            attributes = parts[8]
            id_match = re.search(r'ID=([^;]+)', attributes)
            if not id_match:
                continue
            gene_id = id_match.group(1).strip()
            
            note_match = re.search(r'Note=([^;]+)', attributes)
            if note_match:
                note_val = note_match.group(1).strip()
                notes[gene_id] = note_val
    return notes

def norm_g1(x):
    x = x.strip().replace("gene:", "").replace("transcript:", "").replace("mrna:", "")
    x = re.sub(r"\.\d+$", "", x)
    x = re.sub(r"-mRNA-\d+$", "", x)
    x = re.sub(r"-[1-9]T$", "", x)
    x = re.sub(r"-[1-9]P$", "", x)
    x = re.sub(r"CC(\d+)t(\d+)", r"CC\1g\2", x)
    return x

def norm_g2(x):
    return norm_g1(x)

print("Cargando notas funcionales de GFFs...")
notes1 = load_gff_notes(gff1_path)
notes2 = load_gff_notes(gff2_path)
print(f"{name1}: {len(notes1)} anotaciones cargadas.")
print(f"{name2}: {len(notes2)} anotaciones cargadas.")

# Load kaks
kaks_data = {}
if os.path.exists(kaks_path):
    print(f"Cargando Ka/Ks de: {kaks_path}")
    with open(kaks_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader) # skip header
        for row in reader:
            if not row or len(row) < 5:
                continue
            g1_norm = norm_g1(row[0])
            g2_norm = norm_g2(row[1])
            try:
                ka = float(row[2])
                ks = float(row[3])
                ratio = float(row[4])
                kaks_data[(g1_norm, g2_norm)] = (ka, ks, ratio)
            except ValueError:
                continue

print(f"Parejas con Ka/Ks cargadas: {len(kaks_data)}")

keywords = ['sucrose', 'sps', 'sut', 'susy', 'invertase', 'sweet', 'sugar transporter', 'fructosyltransferase', 'galactosyltransferase']
sugar_genes_mapped = []

if os.path.exists(report_path):
    print(f"Analizando reporte comparativo: {report_path}")
    
    # First pass: load lines and count sizes/sugar counts per block
    rows_data = []
    block_total = {}
    block_sugar = {}
    block_strands = {}
    
    with open(report_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            status = row['Status']
            g1_id = row['Gene1_ID']
            g2_id = row['Gene2_ID']
            
            desc1 = notes1.get(g1_id, notes1.get(g1_id.split('.')[0], ""))
            desc2 = notes2.get(g2_id, notes2.get(g2_id.split('.')[0], ""))
            
            is_sugar = False
            for kw in keywords:
                if kw in desc1.lower() or kw in desc2.lower() or kw in g1_id.lower() or kw in g2_id.lower():
                    is_sugar = True
                    break
            
            block_id = row['Block_ID']
            if status == 'Syntenic':
                block_total[block_id] = block_total.get(block_id, 0) + 1
                if is_sugar:
                    block_sugar[block_id] = block_sugar.get(block_id, 0) + 1
                
                # track strand matches for block orientation
                if block_id not in block_strands:
                    block_strands[block_id] = {'same': 0, 'diff': 0}
                if row.get('Strand1') == row.get('Strand2'):
                    block_strands[block_id]['same'] += 1
                else:
                    block_strands[block_id]['diff'] += 1
            
            rows_data.append((row, desc1, desc2, is_sugar))

    # Second pass: build mapped lists
    for row, desc1, desc2, is_sugar in rows_data:
        if is_sugar:
            status = row['Status']
            g1_id = row['Gene1_ID']
            g2_id = row['Gene2_ID']
            g1_norm = norm_g1(g1_id)
            g2_norm = norm_g2(g2_id)
            ka, ks, ratio = kaks_data.get((g1_norm, g2_norm), (None, None, None))
            
            block_id = row['Block_ID']
            
            # Refine status based on orientation
            if status == 'Syntenic':
                if block_id in block_strands:
                    if block_strands[block_id]['same'] >= block_strands[block_id]['diff']:
                        status = 'Syntenic (Direct)'
                    else:
                        status = 'Syntenic (Inverted)'

            tot = block_total.get(block_id, 0) if 'Syntenic' in status else 0
            sug = block_sugar.get(block_id, 0) if 'Syntenic' in status else 0
            
            sugar_genes_mapped.append({
                'Block_ID': block_id,
                'Block_Total': tot,
                'Block_Sugar': sug,
                'Status': status,
                'G1_ID': g1_id,
                'G1_Chr': row['Chr1'],
                'G1_Desc': desc1,
                'G2_ID': g2_id,
                'G2_Chr': row['Chr2'],
                'G2_Desc': desc2,
                'Ka': ka,
                'Ks': ks,
                'Ratio': ratio
            })

print(f"\nTotal de genes de metabolismo/transporte de azúcar identificados: {len(sugar_genes_mapped)}")

# Group by status
print("\nDistribución de estado sinténico en genes de azúcar:")
states = {}
for g in sugar_genes_mapped:
    states[g['Status']] = states.get(g['Status'], 0) + 1
for k, v in states.items():
    print(f"  {k}: {v}")

# Find genes under positive selection (ratio > 1)
positive_sel = [g for g in sugar_genes_mapped if g['Ratio'] is not None and g['Ratio'] > 1.0]
positive_sel.sort(key=lambda x: x['Ratio'], reverse=True)
print(f"\nGenes de azúcar bajo SELECCIÓN POSITIVA (Ka/Ks > 1.0): {len(positive_sel)}")
for g in positive_sel[:15]:
    print(f"  G1: {g['G1_ID']} ({g['G1_Chr']}) <-> G2: {g['G2_ID']} ({g['G2_Chr']}) | Bloque: {g['Block_ID']} (Total={g['Block_Total']} genes, azúcar={g['Block_Sugar']}) | Ka/Ks = {g['Ratio']:.4f} (Ka={g['Ka']:.4f}, Ks={g['Ks']:.4f}) | Desc: {g['G1_Desc'] or g['G2_Desc']}")

# Find genes under purifying selection (ratio < 0.1)
purifying_sel = [g for g in sugar_genes_mapped if g['Ratio'] is not None and g['Ratio'] < 0.1]
purifying_sel.sort(key=lambda x: x['Ratio'])
print(f"\nGenes de azúcar bajo SELECCIÓN PURIFICADORA (Ka/Ks < 0.1): {len(purifying_sel)}")
for g in purifying_sel[:15]:
    print(f"  G1: {g['G1_ID']} ({g['G1_Chr']}) <-> G2: {g['G2_ID']} ({g['G2_Chr']}) | Bloque: {g['Block_ID']} (Total={g['Block_Total']} genes, azúcar={g['Block_Sugar']}) | Ka/Ks = {g['Ratio']:.4f} (Ka={g['Ka']:.4f}, Ks={g['Ks']:.4f}) | Desc: {g['G1_Desc'] or g['G2_Desc']}")

# Find orphans (PAVs) in name1 related to sucrose/sugar
orphans_g1 = [g for g in sugar_genes_mapped if g['Status'] == 'Orphan_G1']
print(f"\nGenes de azúcar HUÉRFANOS (PAVs) en {name1} (G1): {len(orphans_g1)}")
for g in orphans_g1[:10]:
    print(f"  G1 ID: {g['G1_ID']} | Chr: {g['G1_Chr']} | Desc: {g['G1_Desc']}")

# Find orphans (PAVs) in name2 related to sucrose/sugar
orphans_g2 = [g for g in sugar_genes_mapped if g['Status'] == 'Orphan_G2']
print(f"\nGenes de azúcar HUÉRFANOS (PAVs) en {name2} (G2): {len(orphans_g2)}")
for g in orphans_g2[:10]:
    print(f"  G2 ID: {g['G2_ID']} | Chr: {g['G2_Chr']} | Desc: {g['G2_Desc']}")
