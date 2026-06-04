#!/usr/bin/env python3
"""
Exporta el VCF de la población biparental al formato JoinMap 4 CP (Cross-Pollinator).
Solo exporta los marcadores que están en el mapa genético (.map).

Formato JoinMap CP:
  - name = <population_name>
  - popt = CP
  - nloc = número de loci
  - nind = número de individuos
  - Para cada marcador: nombre, código de segregación, genotipos por individuo
  - Códigos: <lmxll> = padre 1 simplex, <nnxnp> = padre 2 simplex,
             <hkxhk> = ambos simplex, <efxeg> = otros
  - En una población biparental donde padre1 es simplex (10,0→padre2 ausente),
    los marcadores segregan como <lmxll> (si el padre informativo es P1)
    o <nnxnp> (si el informativo es P2).

Para marcadores single-dose desde ACN:
  - Individuo con dosage=1 → 'lm' (heterocigoto para padre simplex)
  - Individuo con dosage=0 → 'll' (homocigoto nulo)
  - Dato faltante → '--'
"""

import argparse
import csv
import math
from collections import defaultdict


def load_map_markers(map_path):
    """Retorna set de IDs de marcadores en el mapa."""
    markers = {}
    with open(map_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            markers[row['Marker']] = row['LinkageGroup']
    return markers


def run(map_path, vcf_path, output_path, pop_name="Biparental_CC1940"):
    print(f"[JoinMap] Cargando marcadores del mapa: {map_path}")
    map_markers = load_map_markers(map_path)
    print(f"  {len(map_markers)} marcadores a exportar")

    samples = []
    loci = []      # list of (marker_id, lg, genotypes_list)

    print(f"[JoinMap] Leyendo VCF: {vcf_path}")
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                cols = line.strip().split('\t')
                samples = cols[9:]
                continue
            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue

            chr_ = cols[0]
            pos  = cols[1]
            mid  = f"{chr_}_{pos}"
            if mid not in map_markers:
                continue

            fmt = cols[8].split(':')
            acn_idx = next((i for i,f in enumerate(fmt) if f == 'ACN'), -1)

            genotypes = []
            for g in cols[9:]:
                if g.startswith('.'):
                    genotypes.append('--')
                    continue
                fields = g.split(':')
                if acn_idx != -1 and len(fields) > acn_idx:
                    parts = fields[acn_idx].split(',')
                    try:
                        alt = int(float(parts[-1]))
                        # Single-dose: 1 → lm (heterozygote), 0 → ll (nulliplex)
                        genotypes.append('lm' if alt >= 1 else 'll')
                    except (ValueError, IndexError):
                        genotypes.append('--')
                else:
                    genotypes.append('--')

            # Check informativeness: need both lm and ll present
            n_lm = genotypes.count('lm')
            n_ll = genotypes.count('ll')
            if n_lm == 0 or n_ll == 0:
                continue  # monomorphic or all missing

            loci.append((mid, map_markers[mid], genotypes))

    print(f"[JoinMap] {len(loci)} loci exportables, {len(samples)} individuos")

    with open(output_path, 'w') as f:
        # Header
        f.write(f"name = {pop_name}\n")
        f.write(f"popt = CP\n")
        f.write(f"nloc = {len(loci)}\n")
        f.write(f"nind = {len(samples)}\n\n")

        # Loci data
        for mid, lg, genos in loci:
            # Segregation type: <lmxll> for simplex parent 1
            seg = "<lmxll>"
            line = f"{mid}\t{seg}\t" + "\t".join(genos) + "\n"
            f.write(line)

        # Individual names
        f.write("\nindividual names:\n")
        for s in samples:
            f.write(f"{s}\n")

        # LG assignment as comments
        f.write("\n; Linkage Group assignments from BioJava:\n")
        for mid, lg, _ in loci:
            f.write(f"; {mid} -> {lg}\n")

    print(f"[JoinMap] Archivo exportado: {output_path}")
    print(f"[JoinMap] Abre este archivo en JoinMap 4/5 con tipo de población CP.")
    print(f"[JoinMap] En JoinMap: File > Open > selecciona el .loc")
    print(f"[JoinMap] Usa LOD=3.0 y REC=0.35 para comparar grupos con BioJava.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Exporta mapa biparental a formato JoinMap CP.")
    p.add_argument("--map",    required=True, help="Archivo .map del mapa genético de BioJava")
    p.add_argument("--vcf",    required=True, help="VCF de la población biparental")
    p.add_argument("--output", required=True, help="Archivo de salida .loc para JoinMap")
    p.add_argument("--pop-name", default="Biparental_CC1940", help="Nombre de la población")
    args = p.parse_args()
    run(args.map, args.vcf, args.output, args.pop_name)
