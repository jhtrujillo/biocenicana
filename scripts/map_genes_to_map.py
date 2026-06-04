#!/usr/bin/env python3
"""
Phase 2: Position candidate genes on the genetic linkage map.

For each candidate gene (known physical coordinates from GFF3), finds the
nearest flanking markers in the genetic map and interpolates the cM position.
Handles chromosome name normalization (chr01 ↔ 1).
"""

import argparse
import csv
import sys
from collections import defaultdict


# ── Candidate genes (physical coordinates from CC 01-1940 GFF3) ─────────────
CANDIDATE_GENES = [
    # Confirmed in both presentation and informe
    ("CC01t039050.1", "chr01", 55618891, 55625629, "Sucrose synthase 2 (SuSy2)",          "Purificadora"),
    ("CC01t042770.1", "chr01", 60506673, 60519537, "Sucrose synthase 4",                  "Presentacion"),
    ("CC01t044570.1", "chr01", 62833534, 62834770, "Alkaline/neutral invertase",          "Purificadora"),
    ("CC03t157420.1", "chr03", 24241458, 24246070, "Neutral/alkaline invertase, mitoc.",  "Presentacion"),
    ("CC03t185480.1", "chr03", 62738478, 62743743, "Sucrose-phosphate synthase (SPS)",    "Neutral"),
    ("CC03t185500.1", "chr03", 62761410, 62766866, "putative SPS",                        "Presentacion"),
    ("CC04t207510.1", "chr04", 25129201, 25134796, "Neutral/alkaline invertase 3, clp",   "Purificadora"),
    ("CC04t208410.1", "chr04", 26667005, 26670298, "Cytosolic invertase 1",               "Purificadora"),
    ("CC04t230890.1", "chr04", 55774345, 55777956, "Sucrose synthase (SuSy)",             "Neutral"),
    ("CC04t196460.1", "chr04",  6117451,  6122149, "putative Alkaline/neutral invertase", "Presentacion"),
    ("CC04t199390.1", "chr04", 10079686, 10086951, "Sucrose-phosphate synthase",          "Presentacion"),
    ("CC04t217770.1", "chr04", 39342372, 39345580, "Cytosolic invertase 1",               "Presentacion"),
    ("CC05t236990.1", "chr05",  8919861,  8922720, "Alkaline/neutral invertase",          "Purificadora"),
    ("CC09t347620.1", "chr09", 31997030, 32014463, "putative SPS4",                       "Neutral"),
    ("CC10t084020.1", "chr10", 27527598, 27532806, "Sucrose-phosphate synthase",          "Neutral"),
    ("CC10t092720.1", "chr10", 40579110, 40588278, "Sucrose synthase",                    "Presentacion"),
    ("CC10t092800.1", "chr10", 40651907, 40661068, "Sucrose synthase (SuSy)",             "Positiva"),
    # New from Ka/Ks informe (not in presentation)
    ("CC06t277320.1", "chr06",      None,     None, "Beta-fructofuranosidase (Invertasa)", "Positiva Ka/Ks=29.2"),
    ("CC03t152280.1", "chr03",      None,     None, "Sugar transporter ERD6-like 5",       "Positiva Ka/Ks=25.9"),
    ("CC03t152310.1", "chr03",      None,     None, "Sugar transporter ERD6-like 5",       "Positiva Ka/Ks=19.6"),
    ("CC04t193040.1", "chr04",      None,     None, "Alkaline/neutral invertase",          "Positiva Ka/Ks=20.6"),
]


def normalize_chr(name: str) -> str:
    """Normalize chromosome names: 'chr01' -> '1', 'chr1' -> '1', '1' -> '1'."""
    n = name.lower().lstrip("chr").lstrip("0") or "0"
    return n


def load_map(map_path: str):
    """
    Returns:
        markers_by_chr: dict[chr_normalized -> sorted list of (pos_phys, cM, LG, marker_id)]
    """
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

    # Sort each chromosome list by physical position
    for chr_norm in markers_by_chr:
        markers_by_chr[chr_norm].sort(key=lambda x: x[0])

    return markers_by_chr


def interpolate_cm(gene_pos: int, markers: list):
    """
    Given a sorted list of (pos_phys, cM, LG, marker_id) on the same chromosome,
    find the two flanking markers and interpolate the cM position.

    Returns:
        (cm_interp, lg, left_marker, right_marker, dist_nearest_bp)
    """
    if not markers:
        return None, None, None, None, None

    # Find flanking markers
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

    # Both flanks available → linear interpolation
    if left is not None and right is not None:
        span_bp = right[0] - left[0]
        span_cm = right[1] - left[1]
        if span_bp > 0:
            frac = (gene_pos - left[0]) / span_bp
            cm_interp = left[1] + frac * span_cm
        else:
            cm_interp = left[1]
        dist_nearest = min(abs(gene_pos - left[0]), abs(gene_pos - right[0]))
        # Use LG of the nearest marker
        lg = left[2] if abs(gene_pos - left[0]) <= abs(gene_pos - right[0]) else right[2]
        return round(cm_interp, 2), lg, left, right, dist_nearest

    # Only left flank (gene is beyond last marker)
    if left is not None:
        dist = gene_pos - left[0]
        return round(left[1], 2), left[2], left, None, dist

    # Only right flank (gene is before first marker)
    dist = right[0] - gene_pos
    return round(right[1], 2), right[2], None, right, dist


def run(map_path: str, output_path: str, gff_path: str = None):
    print("[Phase 2] Loading genetic map...")
    markers_by_chr = load_map(map_path)
    total_markers = sum(len(v) for v in markers_by_chr.values())
    print(f"[Phase 2] {total_markers} markers across {len(markers_by_chr)} chromosomes loaded.")

    # Optionally update gene positions from GFF3
    gff_positions = {}
    if gff_path:
        print(f"[Phase 2] Reading gene positions from GFF3: {gff_path}")
        try:
            with open(gff_path) as f:
                for line in f:
                    if line.startswith("#") or "\tmRNA\t" not in line:
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue
                    attrs = {a.split("=")[0]: a.split("=")[1]
                             for a in cols[8].split(";") if "=" in a}
                    gene_id = attrs.get("ID", "")
                    if gene_id:
                        gff_positions[gene_id] = (cols[0], int(cols[3]), int(cols[4]))
        except FileNotFoundError:
            print(f"[Phase 2] GFF3 not found, using built-in coordinates.")

    results = []
    for gene_id, chr_name, start, end, func, category in CANDIDATE_GENES:
        # Override positions if GFF3 provided
        if gene_id in gff_positions:
            chr_name, start, end = gff_positions[gene_id]

        if start is None:
            results.append({
                "Gene": gene_id, "Funcion": func, "Categoria_KaKs": category,
                "Chr_Phys": chr_name, "Start_Phys": "N/A", "End_Phys": "N/A",
                "Gene_Mid_Phys": "N/A", "LG_Asignado": "N/A", "Pos_cM": "N/A",
                "Marcador_Izq": "N/A", "Pos_Izq_bp": "N/A", "Pos_Izq_cM": "N/A",
                "Marcador_Der": "N/A", "Pos_Der_bp": "N/A", "Pos_Der_cM": "N/A",
                "Distancia_Marcador_Cercano_bp": "N/A", "Estado": "Sin coordenadas físicas"
            })
            continue

        gene_mid = (start + end) // 2
        chr_norm  = normalize_chr(chr_name)
        chr_markers = markers_by_chr.get(chr_norm, [])

        if not chr_markers:
            results.append({
                "Gene": gene_id, "Funcion": func, "Categoria_KaKs": category,
                "Chr_Phys": chr_name, "Start_Phys": start, "End_Phys": end,
                "Gene_Mid_Phys": gene_mid, "LG_Asignado": "N/A", "Pos_cM": "N/A",
                "Marcador_Izq": "N/A", "Pos_Izq_bp": "N/A", "Pos_Izq_cM": "N/A",
                "Marcador_Der": "N/A", "Pos_Der_bp": "N/A", "Pos_Der_cM": "N/A",
                "Distancia_Marcador_Cercano_bp": "N/A",
                "Estado": f"Sin marcadores en {chr_name} (chr_norm={chr_norm})"
            })
            continue

        cm_interp, lg, left, right, dist = interpolate_cm(gene_mid, chr_markers)

        results.append({
            "Gene":                          gene_id,
            "Funcion":                       func,
            "Categoria_KaKs":               category,
            "Chr_Phys":                     chr_name,
            "Start_Phys":                   start,
            "End_Phys":                     end,
            "Gene_Mid_Phys":               gene_mid,
            "LG_Asignado":                  lg if lg else "N/A",
            "Pos_cM":                       cm_interp if cm_interp is not None else "N/A",
            "Marcador_Izq":                 left[3] if left else "N/A",
            "Pos_Izq_bp":                   left[0] if left else "N/A",
            "Pos_Izq_cM":                   left[1] if left else "N/A",
            "Marcador_Der":                 right[3] if right else "N/A",
            "Pos_Der_bp":                   right[0] if right else "N/A",
            "Pos_Der_cM":                   right[1] if right else "N/A",
            "Distancia_Marcador_Cercano_bp": dist if dist is not None else "N/A",
            "Estado":                       "Interpolado" if (left and right) else
                                            ("Extrapolado_izq" if left else "Extrapolado_der")
        })

    # Print summary
    print("\n" + "="*90)
    print(f"{'Gen':<20} {'Chr':<7} {'Mid(bp)':<12} {'LG':>6} {'cM':>8} {'Dist(bp)':>10}  Función")
    print("="*90)
    for r in results:
        print(f"{r['Gene']:<20} {str(r['Chr_Phys']):<7} {str(r['Gene_Mid_Phys']):<12} "
              f"{str(r['LG_Asignado']):>6} {str(r['Pos_cM']):>8} "
              f"{str(r['Distancia_Marcador_Cercano_bp']):>10}  {r['Funcion'][:40]}")
    print("="*90)

    # Write TSV
    fieldnames = ["Gene","Funcion","Categoria_KaKs","Chr_Phys","Start_Phys","End_Phys",
                  "Gene_Mid_Phys","LG_Asignado","Pos_cM","Marcador_Izq","Pos_Izq_bp",
                  "Pos_Izq_cM","Marcador_Der","Pos_Der_bp","Pos_Der_cM",
                  "Distancia_Marcador_Cercano_bp","Estado"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"\n[Phase 2] Results written to: {output_path}")

    # Stats
    positioned = sum(1 for r in results if r["LG_Asignado"] != "N/A")
    no_markers = sum(1 for r in results if "Sin marcadores" in str(r["Estado"]))
    print(f"[Phase 2] {positioned}/{len(results)} genes positioned on the map.")
    if no_markers:
        print(f"[Phase 2] {no_markers} genes on chromosomes with no map markers.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Position candidate genes on the genetic linkage map (Phase 2)."
    )
    parser.add_argument("--map",    required=True, help="Genetic map .map file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--gff",    default=None,  help="GFF3 to override gene coordinates (optional)")
    args = parser.parse_args()

    run(args.map, args.output, args.gff)
