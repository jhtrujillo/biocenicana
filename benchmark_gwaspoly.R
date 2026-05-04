library(GWASpoly)

# Load data
pheno <- read.csv("potato_pheno.csv")
geno <- read.csv("potato_geno.csv")

# GWASpoly requires a specific format
# The first column is marker name, then chrom, then bp, then genotypes
# My potato_geno.csv is already in this format (marker, chrom, bp, ...)

# Prepare GWASpoly data
data <- set.GWASpoly(geno = geno, 
                     pheno = pheno, 
                     ploidy = 4, 
                     format = "numeric", 
                     n.clusters = 1)

# Set kinship (VanRaden is default in GWASpoly if not specified, 
# but we can compute it explicitly to match BioJava)
data <- set.K(data)

# Run GWAS with EMMAX (P3D)
# BioJava uses P3D by default in EMMAX.
# We test Additive and Simplex Dominance (SimplexDom)
data <- GWASpoly(data = data, 
                  models = c("additive", "1-dom"), 
                  traits = "vine.maturity", 
                  params = NULL)

# Extract results
results <- get.GWASpoly(data)
write.csv(results, "gwaspoly_results.csv", row.names = FALSE)

# Print top results for comparison
print(head(results[order(results$P),]), 20)
