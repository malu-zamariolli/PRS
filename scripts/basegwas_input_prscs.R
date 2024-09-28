# Script to generate PRS-CS input (with p-value)
# SNP	A1	A2	BETA	P (tab delimited)
library(data.table)

# Load input files
prsice_gwas <- fread(snakemake@input[["prsice_gwas"]], header = TRUE)

# Select columns
out <- prsice_gwas[, c(3, 4, 5, 8, 6)]
out$A1 <- toupper(out$A1)
out$A2 <- toupper(out$A2)

# Save
fwrite(out, snakemake@output[["out"]],
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)