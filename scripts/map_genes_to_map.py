#!/usr/bin/env python3
"""
Phase 2: Position candidate genes on the genetic linkage map.

Reads gene IDs from a file, extracts physical coordinates from the GFF3,
then finds the nearest flanking markers in the genetic map and interpolates
the cM position.

Handles:
  - Gene-level IDs (CC00g352050) → looks up mRNA children in GFF3
  - Transcript-level IDs (CC01t039050.1) → looks up directly in GFF3
  - Chromosome name normalization (chr01 ↔ 1)
"""

import argparse
import csv
import sys
from collections import defaultdict


def normalize_chr(name: str) -> str:
    """Normalize chromosome names: 'chr01' -> '1', '1' -> '1'."""
    n = name.lower().lstrip("chr").lstrip("0") or "0"
    return n


def load_gene_ids(genes_file: str) -> tuple:
    """
    Load gene IDs from a file (one per line, ignores # comments).
    Supports optional second TSV column for category (e.g. Positiva/Neutral/Purificadora).
    Returns: (list of ids, dict id->categoria)
    """
    ids = []
    categories = {}
    with open(genes_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            gene_id = parts[0]
            ids.append(gene_id)
            if len(parts) >= 2:
                categories[gene_id] = parts[1].strip()
    return ids, categories


def load_gff3_coords(gff_path: str, gene_ids: list) -> dict:
    """
    Extract physical coordinates from GFF3 for a list of gene IDs.
    Supports both gene-level (CCxxgxxxxxx) and transcript-level (CCxxtxxxxxx.x) IDs.
    Returns dict: gene_id -> (chr, start, end, strand, note)
    """
    target_genes  = set(gene_ids)
    # Also build transcript → gene mapping for gene-level IDs
    gene_to_mrna  = defaultdict(list)   # gene_id -> [transcript_id, ...]
    coords        = {}                  # gene_id or transcript_id -> (chr, start, end, strand, note)

    print(f"[Phase 2] Reading GFF3: {gff_path}")
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            feat_type = cols[2]
            if feat_type not in ("gene", "mRNA"):
                continue

            chr_   = cols[0]
            start  = int(cols[3])
            end    = int(cols[4])
            strand = cols[6]
            attrs  = {a.split("=")[0]: a.split("=")[1]
                      for a in cols[8].split(";") if "=" in a}

            feat_id = attrs.get("ID", "")
            parent  = attrs.get("Parent", "")
            note    = attrs.get("Note", "")

            if feat_type == "gene" and feat_id in target_genes:
                coords[feat_id] = (chr_, start, end, strand, note)

            if feat_type == "mRNA":
                # Store transcript coords
                if feat_id in target_genes:
                    coords[feat_id] = (chr_, start, end, strand, note)
                # Link transcript to parent gene
                if parent in target_genes:
                    gene_to_mrna[parent].append(feat_id)
                    # Use transcript coords for the gene if not already set
                    if parent not in coords:
                        coords[parent] = (chr_, start, end, strand, note)

    # For gene-level IDs with no direct match, try first transcript
    for gene_id in target_genes:
        if gene_id not in coords and gene_id in gene_to_mrna:
            first_t = gene_to_mrna[gene_id][0]
            if first_t in coords:
                coords[gene_id] = coords[first_t]

    found = sum(1 for g in target_genes if g in coords)
    print(f"[Phase 2] {found}/{len(target_genes)} genes found in GFF3")
    return coords


def load_map(map_path: str) -> dict:
    """Returns markers_by_chr: dict[chr_normalized -> sorted list of (pos, cM, LG, marker_id)]"""
    markers_by_chr = defaultdict(list)
    with open(map_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                pos = int(row["Pos_Phys"])
                cm  = float(row["Position(cM)"])
            except (ValueError, KeyError):
                continue
            chr_norm = normalize_chr(row["Chr_Phys"])
            markers_by_chr[chr_norm].append((pos, cm, row["LinkageGroup"], row["Marker"]))
    for chr_norm in markers_by_chr:
        markers_by_chr[chr_norm].sort(key=lambda x: x[0])
    return markers_by_chr


def interpolate_cm(gene_pos: int, markers: list):
    """Find flanking markers and interpolate cM position."""
    if not markers:
        return None, None, None, None, None

    left = None
    right = None
    for m in markers:
        if m[0] <= gene_pos:
            left = m
        elif right is None:
            right = m
            break

    if left is None and right is None:
        return None, None, None, None, None

    if left is not None and right is not None:
        span_bp = right[0] - left[0]
        span_cm = right[1] - left[1]
        frac = (gene_pos - left[0]) / span_bp if span_bp > 0 else 0
        cm_interp = left[1] + frac * span_cm
        dist = min(abs(gene_pos - left[0]), abs(gene_pos - right[0]))
        lg = left[2] if abs(gene_pos - left[0]) <= abs(gene_pos - right[0]) else right[2]
        return round(cm_interp, 2), lg, left, right, dist

    if left is not None:
        return round(left[1], 2), left[2], left, None, gene_pos - left[0]

    return round(right[1], 2), right[2], None, right, right[0] - gene_pos


def run(map_path: str, gff_path: str, genes_file: str, output_path: str):
    print(f"[Phase 2] Loading gene list: {genes_file}")
    gene_ids, categories = load_gene_ids(genes_file)
    print(f"[Phase 2] {len(gene_ids)} gene IDs loaded ({len(categories)} with category)")

    coords = load_gff3_coords(gff_path, gene_ids)

    print(f"[Phase 2] Loading genetic map: {map_path}")
    markers_by_chr = load_map(map_path)
    total_markers = sum(len(v) for v in markers_by_chr.values())
    print(f"[Phase 2] {total_markers} markers across {len(markers_by_chr)} chromosomes")

    results = []
    no_coords   = 0
    no_markers  = 0
    positioned  = 0

    for gene_id in gene_ids:
        cat = categories.get(gene_id, "")
        if gene_id not in coords:
            no_coords += 1
            results.append({
                "Gene": gene_id, "Chr_Phys": "N/A", "Start": "N/A", "End": "N/A",
                "Gene_Mid_Phys": "N/A", "LG_Asignado": "N/A", "Pos_cM": "N/A",
                "Marcador_Izq": "N/A", "Pos_Izq_bp": "N/A", "Pos_Izq_cM": "N/A",
                "Marcador_Der": "N/A", "Pos_Der_bp": "N/A", "Pos_Der_cM": "N/A",
                "Dist_Marcador_bp": "N/A", "Funcion": "N/A",
                "Categoria": cat, "Estado": "Sin coordenadas GFF3"
            })
            continue

        chr_name, start, end, strand, note = coords[gene_id]
        gene_mid  = (start + end) // 2
        chr_norm  = normalize_chr(chr_name)
        chr_markers = markers_by_chr.get(chr_norm, [])

        if not chr_markers:
            no_markers += 1
            results.append({
                "Gene": gene_id, "Chr_Phys": chr_name, "Start": start, "End": end,
                "Gene_Mid_Phys": gene_mid, "LG_Asignado": "N/A", "Pos_cM": "N/A",
                "Marcador_Izq": "N/A", "Pos_Izq_bp": "N/A", "Pos_Izq_cM": "N/A",
                "Marcador_Der": "N/A", "Pos_Der_bp": "N/A", "Pos_Der_cM": "N/A",
                "Dist_Marcador_bp": "N/A", "Funcion": note, "Estado": "Sin marcadores en chr"
            })
            continue

        cm, lg, left, right, dist = interpolate_cm(gene_mid, chr_markers)
        positioned += 1

        results.append({
            "Gene":           gene_id,
            "Chr_Phys":       chr_name,
            "Start":          start,
            "End":            end,
            "Gene_Mid_Phys": gene_mid,
            "LG_Asignado":   lg if lg else "N/A",
            "Pos_cM":         cm if cm is not None else "N/A",
            "Marcador_Izq":  left[3]  if left  else "N/A",
            "Pos_Izq_bp":    left[0]  if left  else "N/A",
            "Pos_Izq_cM":    left[1]  if left  else "N/A",
            "Marcador_Der":  right[3] if right else "N/A",
            "Pos_Der_bp":    right[0] if right else "N/A",
            "Pos_Der_cM":    right[1] if right else "N/A",
            "Dist_Marcador_bp": dist if dist is not None else "N/A",
            "Funcion":        note,
            "Categoria":      cat,
            "Estado":         "Interpolado"  if (left and right) else
                              "Extrapolado"
        })

    # Print summary
    print(f"\n[Phase 2] Posicionados: {positioned}/{len(gene_ids)} genes")
    print(f"[Phase 2] Sin coordenadas GFF3: {no_coords}")
    print(f"[Phase 2] Sin marcadores en cromosoma: {no_markers}")

    # LG distribution
    lg_counts = defaultdict(int)
    for r in results:
        if r["LG_Asignado"] != "N/A":
            lg_counts[r["LG_Asignado"]] += 1
    if lg_counts:
        print(f"\n[Phase 2] Distribución por LG (top 15):")
        for lg, n in sorted(lg_counts.items(), key=lambda x: -x[1])[:15]:
            print(f"  {lg}: {n} genes")

    # Write TSV
    fields = ["Gene","Chr_Phys","Start","End","Gene_Mid_Phys","LG_Asignado",
              "Pos_cM","Marcador_Izq","Pos_Izq_bp","Pos_Izq_cM",
              "Marcador_Der","Pos_Der_bp","Pos_Der_cM","Dist_Marcador_bp",
              "Funcion","Categoria","Estado"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    print(f"\n[Phase 2] Results written to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Position genes from a list onto the genetic linkage map."
    )
    parser.add_argument("--map",    required=True, help="Genetic map .map file")
    parser.add_argument("--gff",    required=True, help="GFF3 file with gene coordinates")
    parser.add_argument("--genes",  required=True, help="File with gene IDs (one per line)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()
    run(args.map, args.gff, args.genes, args.output)
