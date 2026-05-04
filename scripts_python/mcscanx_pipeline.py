import os
import subprocess
import sys

def format_gff_for_mcscanx(gff_in, gff_out):
    """Converts standard GFF3 to MCScanX format: chr gene start end"""
    print("  - Formateando GFF: {0}".format(gff_in))
    with open(gff_in, 'r') as fin, open(gff_out, 'w') as fout:
        for line in fin:
            if line.startswith('#') or '\tgene\t' not in line:
                continue
            parts = line.split('\t')
            chrom = parts[0]
            start = parts[3]
            end = parts[4]
            # Extract ID
            attr = parts[8]
            gene_id = ""
            if 'ID=' in attr:
                gene_id = attr.split('ID=')[1].split(';')[0]
            
            if gene_id:
                fout.write("{0}\t{1}\t{2}\t{3}\n".format(chrom, gene_id, start, end))

def run_pipeline(genome1_name, gff1, cds1, genome2_name, gff2, cds2, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    base_name = "{0}_vs_{1}".format(genome1_name, genome2_name)
    
    # 1. Format GFFs
    merged_gff = os.path.join(output_dir, "{0}.gff".format(base_name))
    format_gff_for_mcscanx(gff1, merged_gff)
    
    # Append second GFF
    temp_gff2 = "temp_gff2.txt"
    format_gff_for_mcscanx(gff2, temp_gff2)
    with open(merged_gff, 'a') as f:
        with open(temp_gff2, 'r') as f2:
            f.write(f2.read())
    if os.path.exists(temp_gff2):
        os.remove(temp_gff2)

    print("\n!!! REQUIERE BLAST/DIAMOND !!!")
    print("Para completar este paso, ejecuta:")
    print("  makeblastdb -in {0} -dbtype nucl".format(cds1))
    print("  blastn -query {0} -db {1} -out {2}.blast -evalue 1e-10 -outfmt 6 -num_threads 8".format(cds2, cds1, os.path.join(output_dir, base_name)))
    
    print("\n[INFO] Una vez tengas el archivo .blast, ejecuta MCScanX:")
    print("  MCScanX {0}".format(os.path.join(output_dir, base_name)))

if __name__ == "__main__":
    run_pipeline(
        "sorghum", "benchmarks/genomas/sorghum/sorghum.gff3", "benchmarks/genomas/sorghum/sorghum.cds.fa",
        "1940", "benchmarks/genomas/1940/CC-01-1940.gff3", "benchmarks/genomas/1940/CC-01-1940.cds.fna",
        "benchmarks/genomica_comparativa/sorghum_vs_1940"
    )
