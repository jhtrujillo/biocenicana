import random
import os

FUNCTIONS = [
    ("kinase domain protein", "GO:0004672", "PF00069"),
    ("disease resistance protein (NBS-LRR)", "GO:0043207", "PF00931"),
    ("cytochrome P450", "GO:0004497", "PF00067"),
    ("transcription factor MYB", "GO:0003677", "PF00249"),
    ("heat shock protein 70", "GO:0009408", "PF00012"),
    ("UDP-glucosyltransferase", "GO:0008194", "PF00201"),
    ("ATP binding cassette transporter", "GO:0042626", "PF00005"),
    ("photosystem II reaction center", "GO:0009523", "PF00124"),
    ("hypothetical protein", "", ""),
    ("unknown protein", "", "")
]

def generate_mock_gff_and_annot(gff_file, annot_file, prefix, num_genes=50, chrom="Chr01"):
    with open(gff_file, 'w') as gf, open(annot_file, 'w') as af:
        gf.write("##gff-version 3\n")
        af.write("GeneID\tDescription\tGO_Term\tPFAM\n")
        for i in range(1, num_genes + 1):
            start = i * 1000
            end = start + 500
            gene_id = f"{prefix}_{chrom}_g{i:04d}"
            
            func, go, pfam = random.choice(FUNCTIONS)
            attr = f"ID={gene_id};Name={gene_id};description={func}"
            if go: attr += f";Ontology_term={go}"
            if pfam: attr += f";Dbxref=PFAM:{pfam}"
            
            gf.write(f"{chrom}\tMockSource\tgene\t{start}\t{end}\t.\t+\t.\t{attr}\n")
            gf.write(f"{chrom}\tMockSource\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}.1;Parent={gene_id}\n")
            af.write(f"{gene_id}\t{func}\t{go}\t{pfam}\n")

def generate_mock_fasta(filename, prefix, num_genes=50, chrom="Chr01"):
    with open(filename, 'w') as f:
        for i in range(1, num_genes + 1):
            gene_id = f"{prefix}_{chrom}_g{i:04d}.1"
            seq = "".join(random.choices("ACGT", k=300))
            f.write(f">{gene_id}\n{seq}\n")

def generate_mock_collinearity(filename, pref1, pref2, num_genes=50, chrom1="Chr01", chrom2="Chr01"):
    with open(filename, 'w') as f:
        f.write("## Alignment 0: score=1000.0 e_value=0.0 N=50\n")
        for i in range(1, num_genes + 1):
            g1 = f"{pref1}_{chrom1}_g{i:04d}"
            # Simulate some conservation with slight variation
            g2 = f"{pref2}_{chrom2}_g{i:04d}"
            f.write(f"{i-1}: {g1} {g2} 0.0\n")

def generate_mock_kaks(filename, pref1, pref2, num_genes=50, chrom1="Chr01", chrom2="Chr01"):
    with open(filename, 'w') as f:
        f.write("#Gene1\tGene2\tKa\tKs\tKa/Ks\n")
        for i in range(1, num_genes + 1):
            g1 = f"{pref1}_{chrom1}_g{i:04d}"
            g2 = f"{pref2}_{chrom2}_g{i:04d}"
            ks = random.uniform(0.01, 0.5)
            # Most genes are purifying (Ka/Ks < 1)
            # Some are neutral (Ka/Ks ~ 1)
            # A few are positive (Ka/Ks > 1)
            r_type = random.random()
            if r_type < 0.8: kaks = random.uniform(0.05, 0.5)  # Purifying
            elif r_type < 0.95: kaks = random.uniform(0.8, 1.2) # Neutral
            else: kaks = random.uniform(1.5, 3.5)              # Positive
            ka = ks * kaks
            f.write(f"{g1}\t{g2}\t{ka:.4f}\t{ks:.4f}\t{kaks:.4f}\n")

# Directory for simulation
sim_dir = "simulation_data"
if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)

# Genomes: R570, CC01, Spontaneum
genomes = ["R570", "CC01", "Spont"]

print("[Sim] Generating GFFs, Annotations and FASTA...")
for g in genomes:
    generate_mock_gff_and_annot(f"{sim_dir}/{g}.gff", f"{sim_dir}/{g}_annot.tsv", g)
    generate_mock_fasta(f"{sim_dir}/{g}.cds.fa", g)

print("[Sim] Generating Collinearity and Ka/Ks...")
# R570 vs CC01 (highly syntenic)
generate_mock_collinearity(f"{sim_dir}/R570_vs_CC01.collinearity", "R570", "CC01")
generate_mock_kaks(f"{sim_dir}/R570_vs_CC01.kaks.tsv", "R570", "CC01")

# R570 vs Spontaneum (mostly syntenic)
generate_mock_collinearity(f"{sim_dir}/R570_vs_Spont.collinearity", "R570", "Spont")
generate_mock_kaks(f"{sim_dir}/R570_vs_Spont.kaks.tsv", "R570", "Spont")

print(f"[Sim] Done! Data available in '{sim_dir}'")
