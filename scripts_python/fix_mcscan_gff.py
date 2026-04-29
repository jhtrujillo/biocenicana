
import re

def simplify_gff(input_file, output_file, genome_prefix=""):
    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            # For Sorghum, we need the transcript/protein ID (EER...)
            if "sorghum" in input_file.lower():
                if parts[2] == 'mRNA' or parts[2] == 'transcript':
                    attributes = parts[8]
                    # Try to find transcript_id or ID
                    m = re.search(r'transcript_id=([^;]+)', attributes)
                    if not m:
                        m = re.search(r'ID=transcript:([^;]+)', attributes)
                    
                    if m:
                        gene_id = m.group(1)
                        out.write(f"{parts[0]}\t{gene_id}\t{parts[3]}\t{parts[4]}\n")
            
            # For 1940, we need the mRNA ID (CC01t...)
            else:
                if parts[2] == 'mRNA':
                    attributes = parts[8]
                    m = re.search(r'ID=([^;]+)', attributes)
                    if m:
                        gene_id = m.group(1)
                        out.write(f"{parts[0]}\t{gene_id}\t{parts[3]}\t{parts[4]}\n")

# Main execution
simplify_gff('benchmarks/genomas/sorghum/sorghum.gff3', 'temp_sorghum.gff')
simplify_gff('benchmarks/genomas/1940/CC-01-1940.gff3', 'temp_1940.gff')

# Combine them into the final file for MCScanX
with open('benchmarks/genomica_comparativa/sorghum_vs_1940/sorghum_vs_1940.gff', 'w') as out:
    with open('temp_sorghum.gff') as f: out.write(f.read())
    with open('temp_1940.gff') as f: out.write(f.read())

print("Archivo GFF para MCScanX regenerado con éxito.")
