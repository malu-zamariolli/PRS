# Format BaseGWAS for PRS
library(data.table)
library(dplyr)

# Output needs to look like this: 
# CHR, BP, SNP, A1, A2, P, INFO (case sensitive) and OR / BETA (case insensitive), MAF

# Load input
gwas <- fread(snakemake@input[["gwas"]], header = T)

# Convert logOR to OR
gwas$OR <- exp(gwas$LogOR)

# Create missing columns
gwas$Dummy <- c("DummyColumn")
gwas$Dummy1 <- c("DummayColumn1")
gwas$INFO <- c(1) # So nothing gets filtered out

# Select columns
gwas <- gwas %>% select("Dummy", "Dummy1", "MarkerName", "A1", "A2", "P", "INFO", "OR", "Freq")

# Rename columns
colnames(gwas) <- c("Dummy", "Dummy1", "SNP", "A1", "A2", "P", "INFO", "OR", "MAF")

print(paste0("Number of SNPs in DepressionHoward: ", nrow(gwas)))

# Save output
output <- snakemake@output[["out"]]
data.table::fwrite(gwas, file = output, quote = FALSE, sep = " ")