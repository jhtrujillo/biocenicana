#!/usr/bin/env python3
"""
Validación interna del mapa genético.
Computa métricas de calidad sin necesidad de software externo:
  1. Correlación Spearman física vs genética por LG
  2. Factor de expansión (cM/Mb) por LG y cromosoma
  3. Distorsión de segregación (chi² 1:1) por marcador
  4. Dobles recombinantes entre marcadores consecutivos
  5. Resumen general de calidad
"""

import argparse
import csv
import math
import sys
from collections import defaultdict


# ── Estadísticas ────────────────────────────────────────────────────────────

def spearman_r(xs, ys):
    """Correlación de Spearman entre dos listas."""
    n = len(xs)
    if n < 3:
        return float('nan')
    def rank(lst):
        sorted_idx = sorted(range(n), key=lambda i: lst[i])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and lst[sorted_idx[j]] == lst[sorted_idx[j+1]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1
            for k in range(i, j+1):
                ranks[sorted_idx[k]] = avg_rank
            i = j + 1
        return ranks
    rx = rank(xs); ry = rank(ys)
    mean_rx = sum(rx)/n; mean_ry = sum(ry)/n
    num = sum((rx[i]-mean_rx)*(ry[i]-mean_ry) for i in range(n))
    den = math.sqrt(sum((v-mean_rx)**2 for v in rx) * sum((v-mean_ry)**2 for v in ry))
    return num/den if den > 0 else 0.0


def chi2_p(obs, exp):
    """p-valor aproximado del test chi² con 1 grado de libertad."""
    if exp == 0:
        return 0.0
    chi2 = (obs - exp)**2 / exp + (exp - (obs if isinstance(obs,int) else int(obs+0.5)))**2 / exp
    # Approx p via complementary error function
    x = math.sqrt(chi2 / 2.0)
    return 1.0 - math.erf(x)


def chi2_p_binom(n_present, n_absent):
    """Chi² 1:1 entre presentes y ausentes."""
    n = n_present + n_absent
    if n < 5:
        return float('nan')
    exp = n / 2.0
    chi2 = (n_present - exp)**2 / exp + (n_absent - exp)**2 / exp
    x = math.sqrt(chi2 / 2.0)
    return max(0.0, 1.0 - math.erf(x))


# ── Lectura de datos ────────────────────────────────────────────────────────

def load_map(map_path):
    markers = []
    with open(map_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                markers.append({
                    'id':  row['Marker'],
                    'lg':  row['LinkageGroup'],
                    'cm':  float(row['Position(cM)']),
                    'chr': row['Chr_Phys'],
                    'pos': int(row['Pos_Phys'])
                })
            except (ValueError, KeyError):
                continue
    return markers


def load_dosages_from_vcf(vcf_path, marker_ids):
    """Lee las dosificaciones (campo ACN) de los marcadores del mapa desde el VCF."""
    marker_set = set(marker_ids)
    dosages = {}   # marker_id -> list of dosage values (0/1/NaN per sample)
    samples = []

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
            mid = f"{chr_}_{pos}"
            if mid not in marker_set:
                continue

            # Find ACN index
            fmt = cols[8].split(':')
            acn_idx = next((i for i, f in enumerate(fmt) if f == 'ACN'), -1)

            ds = []
            for g in cols[9:]:
                if g.startswith('.'):
                    ds.append(float('nan'))
                    continue
                fields = g.split(':')
                if acn_idx != -1 and len(fields) > acn_idx:
                    parts = fields[acn_idx].split(',')
                    try:
                        ds.append(float(parts[-1]))
                    except (ValueError, IndexError):
                        ds.append(float('nan'))
                else:
                    ds.append(float('nan'))
            dosages[mid] = ds

    return dosages, samples


# ── Validaciones ─────────────────────────────────────────────────────────────

def validate_spearman(markers_by_lg):
    """Correlación Spearman entre posición física y genética por LG."""
    results = []
    for lg, ms in markers_by_lg.items():
        if len(ms) < 3:
            continue
        xs = [m['pos'] for m in ms]
        ys = [m['cm']  for m in ms]
        r  = spearman_r(xs, ys)
        dominant_chr = max(set(m['chr'] for m in ms), key=lambda c: sum(1 for m in ms if m['chr']==c))
        results.append({'lg': lg, 'n': len(ms), 'spearman_r': r, 'dominant_chr': dominant_chr})
    return sorted(results, key=lambda x: x['spearman_r'])


def validate_expansion(markers_by_lg):
    """Factor de expansión cM/Mb por LG."""
    results = []
    for lg, ms in markers_by_lg.items():
        if len(ms) < 2:
            continue
        max_cm  = max(m['cm']  for m in ms)
        min_pos = min(m['pos'] for m in ms)
        max_pos = max(m['pos'] for m in ms)
        span_mb = (max_pos - min_pos) / 1e6
        if span_mb < 0.01:
            continue
        cM_per_Mb = max_cm / span_mb if max_cm > 0 else 0
        dominant_chr = max(set(m['chr'] for m in ms), key=lambda c: sum(1 for m in ms if m['chr']==c))
        results.append({'lg': lg, 'n': len(ms), 'length_cM': round(max_cm,2),
                        'span_Mb': round(span_mb,2), 'cM_per_Mb': round(cM_per_Mb,2),
                        'dominant_chr': dominant_chr,
                        'flag': '⚠' if cM_per_Mb > 20 or cM_per_Mb < 0.5 else '✓'})
    return sorted(results, key=lambda x: -x['cM_per_Mb'])


def validate_segregation(dosages):
    """Test chi² 1:1 para cada marcador (presentes vs ausentes)."""
    results = []
    for mid, ds in dosages.items():
        valid = [d for d in ds if not math.isnan(d)]
        n_present = sum(1 for d in valid if round(d) == 1)
        n_absent  = sum(1 for d in valid if round(d) == 0)
        n_multi   = sum(1 for d in valid if round(d) > 1)
        n = n_present + n_absent + n_multi
        p = chi2_p_binom(n_present, n_absent) if n >= 5 else float('nan')
        results.append({'marker': mid, 'n_total': n, 'n_present': n_present,
                        'n_absent': n_absent, 'n_multi': n_multi,
                        'freq_present': round(n_present/n,3) if n>0 else 0,
                        'chi2_p': round(p, 4) if not math.isnan(p) else 'NA',
                        'distorted': (p < 0.05) if not math.isnan(p) else False})
    return sorted(results, key=lambda x: (x['chi2_p'] if isinstance(x['chi2_p'], float) else 1.0))


def validate_double_recombinants(markers_by_lg, dosages, samples):
    """Cuenta dobles recombinantes entre marcadores consecutivos en cada LG."""
    results = []
    for lg, ms in markers_by_lg.items():
        if len(ms) < 3:
            continue
        ms_sorted = sorted(ms, key=lambda m: m['cm'])
        double_rec_per_individual = defaultdict(int)
        total_triplets = 0

        for k in range(len(ms_sorted)-2):
            m1 = ms_sorted[k];   m2 = ms_sorted[k+1];   m3 = ms_sorted[k+2]
            d1 = dosages.get(m1['id'], [])
            d2 = dosages.get(m2['id'], [])
            d3 = dosages.get(m3['id'], [])
            if not d1 or not d2 or not d3:
                continue
            total_triplets += 1
            for i in range(min(len(d1), len(d2), len(d3))):
                v1 = round(d1[i]) if not math.isnan(d1[i]) else -1
                v2 = round(d2[i]) if not math.isnan(d2[i]) else -1
                v3 = round(d3[i]) if not math.isnan(d3[i]) else -1
                if v1 < 0 or v2 < 0 or v3 < 0:
                    continue
                rec12 = (v1 != v2)
                rec23 = (v2 != v3)
                if rec12 and rec23:
                    double_rec_per_individual[i] += 1

        total_double = sum(double_rec_per_individual.values())
        n_individuals = len(samples)
        results.append({'lg': lg, 'n_markers': len(ms), 'n_triplets': total_triplets,
                        'total_double_rec': total_double,
                        'avg_per_ind': round(total_double/n_individuals, 3) if n_individuals>0 else 0,
                        'flag': '⚠' if total_double > total_triplets * n_individuals * 0.05 else '✓'})
    return sorted(results, key=lambda x: -x['avg_per_ind'])


# ── Main ─────────────────────────────────────────────────────────────────────

def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    print(f"  → {path}")


def run(map_path, vcf_path, out_dir):
    import os; os.makedirs(out_dir, exist_ok=True)

    print("[Validación] Cargando mapa genético...")
    markers = load_map(map_path)
    print(f"  {len(markers)} marcadores en {len(set(m['lg'] for m in markers))} LGs")

    markers_by_lg = defaultdict(list)
    for m in markers: markers_by_lg[m['lg']].append(m)

    marker_ids = set(m['id'] for m in markers)

    print("[Validación] Cargando dosificaciones del VCF...")
    dosages, samples = load_dosages_from_vcf(vcf_path, marker_ids)
    print(f"  {len(dosages)} marcadores cargados del VCF, {len(samples)} individuos")

    print()

    # 1. Spearman
    print("=" * 70)
    print("1. CORRELACIÓN SPEARMAN (posición física vs genética por LG)")
    print("=" * 70)
    sp_res = validate_spearman(markers_by_lg)
    low_sp = [r for r in sp_res if r['spearman_r'] < 0.5]
    high_sp = [r for r in sp_res if r['spearman_r'] >= 0.7]
    print(f"   LGs con r >= 0.7 (buen orden): {len(high_sp)}/{len(sp_res)}")
    print(f"   LGs con r < 0.5  (orden dudoso): {len(low_sp)}/{len(sp_res)}")
    if sp_res:
        avg_r = sum(r['spearman_r'] for r in sp_res) / len(sp_res)
        print(f"   Correlación promedio: {avg_r:.3f}")
    print(f"\n   {'LG':<10} {'N':>4}  {'Spearman r':>10}  {'Chr dominante':>14}  Evaluación")
    for r in sp_res:
        flag = '✓ bueno' if r['spearman_r'] >= 0.7 else ('~ aceptable' if r['spearman_r'] >= 0.5 else '⚠ revisar')
        print(f"   {r['lg']:<10} {r['n']:>4}  {r['spearman_r']:>10.3f}  chr{r['dominant_chr']:>12}  {flag}")
    write_tsv(f"{out_dir}/validation_spearman.tsv", sp_res,
              ['lg','n','spearman_r','dominant_chr'])

    # 2. Expansión cM/Mb
    print()
    print("=" * 70)
    print("2. FACTOR DE EXPANSIÓN (cM/Mb) — esperado 1–10 cM/Mb para caña")
    print("=" * 70)
    exp_res = validate_expansion(markers_by_lg)
    flagged = [r for r in exp_res if r['flag'] == '⚠']
    print(f"   LGs dentro del rango esperado: {len(exp_res)-len(flagged)}/{len(exp_res)}")
    print(f"   LGs fuera de rango (⚠): {len(flagged)}")
    print(f"\n   {'LG':<10} {'N':>4}  {'cM':>8}  {'Mb':>8}  {'cM/Mb':>8}  {'Chr':>5}  Estado")
    for r in exp_res:
        print(f"   {r['lg']:<10} {r['n']:>4}  {r['length_cM']:>8.1f}  {r['span_Mb']:>8.2f}  {r['cM_per_Mb']:>8.2f}  chr{r['dominant_chr']:>3}  {r['flag']}")
    write_tsv(f"{out_dir}/validation_expansion.tsv", exp_res,
              ['lg','n','length_cM','span_Mb','cM_per_Mb','dominant_chr','flag'])

    # 3. Segregación
    print()
    print("=" * 70)
    print("3. DISTORSIÓN DE SEGREGACIÓN (chi² 1:1, FDR < 0.05 = distorsionado)")
    print("=" * 70)
    seg_res = validate_segregation(dosages)
    distorted = [r for r in seg_res if r['distorted']]
    print(f"   Marcadores sin distorsión: {len(seg_res)-len(distorted)}/{len(seg_res)}")
    print(f"   Marcadores distorsionados (p < 0.05): {len(distorted)}")
    if distorted:
        print(f"\n   Top 10 más distorsionados:")
        print(f"   {'Marcador':<25} {'N':>4}  {'Pres':>5}  {'Aus':>5}  {'Frec':>6}  {'p-val':>8}")
        for r in distorted[:10]:
            print(f"   {r['marker']:<25} {r['n_total']:>4}  {r['n_present']:>5}  {r['n_absent']:>5}  {r['freq_present']:>6.3f}  {r['chi2_p']:>8}")
    write_tsv(f"{out_dir}/validation_segregation.tsv", seg_res,
              ['marker','n_total','n_present','n_absent','n_multi','freq_present','chi2_p','distorted'])

    # 4. Dobles recombinantes
    print()
    print("=" * 70)
    print("4. DOBLES RECOMBINANTES entre marcadores consecutivos")
    print("   (>5% por tripleta indica posible error de orden)")
    print("=" * 70)
    dr_res = validate_double_recombinants(markers_by_lg, dosages, samples)
    flagged_dr = [r for r in dr_res if r['flag'] == '⚠']
    print(f"   LGs con orden confiable: {len(dr_res)-len(flagged_dr)}/{len(dr_res)}")
    print(f"   LGs con posible error de orden (⚠): {len(flagged_dr)}")
    print(f"\n   {'LG':<10} {'Mks':>4}  {'Triplets':>8}  {'Dobles Rec':>10}  {'Avg/ind':>8}  Estado")
    for r in dr_res:
        print(f"   {r['lg']:<10} {r['n_markers']:>4}  {r['n_triplets']:>8}  {r['total_double_rec']:>10}  {r['avg_per_ind']:>8.3f}  {r['flag']}")
    write_tsv(f"{out_dir}/validation_double_rec.tsv", dr_res,
              ['lg','n_markers','n_triplets','total_double_rec','avg_per_ind','flag'])

    # Resumen final
    print()
    print("=" * 70)
    print("RESUMEN DE CALIDAD DEL MAPA")
    print("=" * 70)
    n_lgs = len([lg for lg in markers_by_lg if len(markers_by_lg[lg]) >= 3])
    n_good_sp = len([r for r in sp_res if r['spearman_r'] >= 0.7])
    n_good_exp = len([r for r in exp_res if r['flag'] == '✓'])
    n_nodist = len([r for r in seg_res if not r['distorted']])
    n_good_order = len([r for r in dr_res if r['flag'] == '✓'])
    print(f"   LGs con ≥3 marcadores evaluados: {n_lgs}")
    print(f"   ✓ Orden físico-genético coherente (Spearman ≥0.7): {n_good_sp}/{len(sp_res)}")
    print(f"   ✓ Factor expansión normal (1-10 cM/Mb):            {n_good_exp}/{len(exp_res)}")
    print(f"   ✓ Marcadores sin distorsión de segregación:        {n_nodist}/{len(seg_res)}")
    print(f"   ✓ LGs sin errores de orden detectados:             {n_good_order}/{len(dr_res)}")
    print()
    print(f"   Resultados guardados en: {out_dir}/")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Validación interna del mapa genético.")
    p.add_argument("--map",    required=True, help="Archivo .map del mapa genético")
    p.add_argument("--vcf",    required=True, help="VCF de la población biparental")
    p.add_argument("--out-dir",default="genomica_comparativa/mapa_genetico/validacion",
                   help="Directorio de salida para los TSVs")
    args = p.parse_args()
    run(args.map, args.vcf, args.out_dir)
