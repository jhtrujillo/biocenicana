import random
import os

def generate_vcf(filename, num_samples, num_snps):
    chromosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
    bases = ['A', 'C', 'G', 'T']

    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, 'w') as f:
        # Write VCF Headers
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=SyntheticGenerator\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")

        # Write Column Headers
        header_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        sample_names = [f"Sample_{i+1}" for i in range(num_samples)]
        f.write('\t'.join(header_cols + sample_names) + '\n')

        # Generate SNPs
        for i in range(num_snps):
            chrom = random.choice(chromosomes)
            pos = i * 1000 + random.randint(1, 100)
            snp_id = f"snp_{i+1}"
            
            ref = random.choice(bases)
            alt_candidates = [b for b in bases if b != ref]
            alt = random.choice(alt_candidates)
            
            # Make ~10% of variants Indels to test Ts/Tv vs Indels chart
            if random.random() < 0.1:
                alt = ref + random.choice(bases)
            
            qual = "."
            filter_val = "PASS"
            info = "."
            fmt = "GT:DP"

            row = [chrom, str(pos), snp_id, ref, alt, qual, filter_val, info, fmt]

            # Generate sample data
            for _ in range(num_samples):
                # 5% missingness rate globally
                if random.random() < 0.05:
                    row.append("./.:.")
                else:
                    # Random genotype 0/0, 0/1, 1/1
                    gt = random.choice(["0/0", "0/1", "1/1"])
                    
                    # Poisson-like depth centered around 30, with some variance
                    dp = max(1, int(random.gauss(30, 10)))
                    row.append(f"{gt}:{dp}")

            f.write('\t'.join(row) + '\n')

if __name__ == "__main__":
    vcf_path = "example_vcfs/large_synthetic.vcf"
    print(f"Generating VCF: {vcf_path} with 200 samples and 5000 SNPs...")
    generate_vcf(vcf_path, 200, 5000)
    print("Done!")
