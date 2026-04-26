import random
chroms = ["Chr01"]
with open("simulation_data/CC01_sim.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n")
    for _ in range(200):
        pos = random.randint(1000, 50000)
        f.write(f"Chr01\t{pos}\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\n")
print("VCF simulation done.")
