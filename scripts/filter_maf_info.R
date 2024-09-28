gwas_base = snakemake@input[["gwas_base"]]

# Filter GWAS for MAF and INFO (imputation score)
library(data.table)
# Read in file
dat <- fread(gwas_base)
# Filter out SNPs
result <- dat[INFO > 0.8 & MAF > 0.01]
# Output the gz file
fwrite(result, snakemake@output[[1]], sep="\t")