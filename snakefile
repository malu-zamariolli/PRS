#--------------------------------------- Polygenic Risk Scores ---------------------------------------------------------------------------#
# PRSice2 and PRS-CS 


#---------------------------------------  Key-value pairs --------------------------------------------------------------------------------#
# A. 
    # Define number of individuals per GWAS (for rule prs_cs)
import pandas as pd
df = pd.read_csv("data/info/BaseGWAS_TargetTrait_info.csv")
n_gwas = df.set_index("BaseGWAS")["N_GWAS"].to_dict()
df_combination = df.set_index("BaseTarget")["BaseGWAS"].to_dict()
    # Define arguments --stat and --binary target for rule prsice
df_stat = df.set_index("BaseTarget")["STAT"].to_dict()
df_bin_target = df.set_index("BaseTarget")["BinaryTargetOption"].to_dict()

#---------------------------------------------  Definition of WILDCARDS ------------------------------------------------------------------#
# B. Define wildcard
DATASETS,EXT = glob_wildcards("../QC/Binary_final/{dataset}_QCed_final.{ext}")
BASEGWAS = df["BaseGWAS"].tolist()
BASETARGET = df["BaseTarget"].tolist()
# df was defined above (key-value pairs)

#---------------------------------------------  Placeholder input files ------------------------------------------------------------------#
# Create placeholder file (empty file to guide trait and covariate rule)
empty_file_names = expand("empty/{BaseTarget}.txt", zip, BaseTarget = BASETARGET)

for file_name in empty_file_names:
    # Create any required directories
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    # Create the empty file
    with open(file_name, 'w') as f:
        pass  # An empty file

#---------------------------------------------  List expected output ------------------------------------------------------------------#
# This needs to be done, otherwise snakemake will expect output files with all wildcards combinations
from itertools import product
combinations_A_B = list(zip(BASEGWAS, BASETARGET))
combinations_A_B_C = [(BaseGWAS, BaseTarget, dataset) for (BaseGWAS, BaseTarget), dataset in product(combinations_A_B, DATASETS)]

# Expand the filenames with combined A, B, and C
prsice_output = [f"Results/PRSice2/{BaseGWAS}_{BaseTarget}_{dataset}.best" for BaseGWAS, BaseTarget, dataset in combinations_A_B_C]

coef_output = [f"Results/PRSCS/{dataset}_{BaseGWAS}_{BaseTarget}_regression_coef.txt" for BaseGWAS, BaseTarget, dataset in combinations_A_B_C]
r2_output = [f"Results/PRSCS/{dataset}_{BaseGWAS}_{BaseTarget}_regression_r_squared.txt" for BaseGWAS, BaseTarget, dataset in combinations_A_B_C]

     
#---------------------------------------  Input Function -------------------------------------------------------------------------------#
# C. Define input function for rule input_gwas
import glob
from pathlib import Path

def wildcard_input(wildcards):
    # Define the directory where the files are located
    data_dir = "../BaseDataGWAS"
    # Construct the wildcard pattern to search for files
    wildcard_pattern = f"{data_dir}/gwas_{wildcards.BaseGWAS}*"
    # Use glob to find files that match the pattern
    input_files = list(Path(".").glob(wildcard_pattern))
    if not input_files:
        raise Exception(f"No matching input files found for wildcard: {wildcards.BaseGWAS}")
    return input_files

#---------------------------------------  Rule all -------------------------------------------------------------------------------------#
# D. Rule all
rule all:
    input:
        expand("Results/PRSCS/{dataset}_{BaseGWAS}_PRSCS_Score.profile", dataset=DATASETS,BaseGWAS=BASEGWAS),
        prsice_output,
        coef_output,
        r2_output 
   # Rule all - listing files that are output and not input to another rule

# ------------------------------------- Base Data GWAS ---------------------------------------------------------------------------------#
rule input_gwas:
    input: 
        gwas = wildcard_input
    output:
        out = "data/BaseGWAS_input/gwas_{BaseGWAS}.gz"
    conda: "./envs/R.yaml"
    script:
        "scripts/BaseDataGWAS/{wildcards.BaseGWAS}_input.R"
 
rule filter_maf_info:
    input: 
        gwas_base = "data/BaseGWAS_input/gwas_{BaseGWAS}.gz"
    output: 
      temp("data/BaseGWAS_input/baseGWAS_files/gwas_filter_{BaseGWAS}.gz")
    conda: "./envs/R.yaml"
    script: 
        "scripts/filter_maf_info.R"
    # Filter BaseGWAS to keep SNPs with INFO > 0.8 (imputation score) & MAF > 0.01

rule duplicate:
    input: 
        gwas_base_filtered = "data/BaseGWAS_input/baseGWAS_files/gwas_filter_{BaseGWAS}.gz"
    output: 
        temp("data/BaseGWAS_input/baseGWAS_files/gwas_nodup_{BaseGWAS}.gz")
    shell:
        "sh scripts/duplicate_SNPs.sh {input.gwas_base_filtered} {output}"
    #Remove duplicated SNPs

rule ambiguous:
    input:
        gwas_base_nodup = "data/BaseGWAS_input/baseGWAS_files/gwas_nodup_{BaseGWAS}.gz"
    output:
        "data/BaseGWAS_input/gwas_prsice_{BaseGWAS}.gz"
    shell:
        """
        gunzip -c {input.gwas_base_nodup} |\
        awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {{print}}' |\
        gzip > {output}
        """
    # Remove ambiguous SNPs: if strand is unknown (+/-), SNPs with complementary alleles, either C/G or A/T can't be resolved)

# -------------------------------------------- Target Data ------------------------------------------------------------------------------#

# QC was previously done
# Duplicated SNPs are not expected (if it happens, add a rule to resolve them)

# --------------------------------------------- PRSice-2 --------------------------------------------------------------------------------#
rule covariate:
    input:
        pcs = "../QC/PCA/{dataset}_pca.eigenvec",
        ids = "../QC/Binary_final/{dataset}_QCed_final.fam",
        info = empty_file_names
    output:
        cov_prs = "data/covariates/{dataset}_{BaseTarget}.covariate"
    conda: "./envs/R.yaml"
    script:
        "scripts/Traits_covariates/{wildcards.dataset}_{wildcards.BaseTarget}_covariates.R"
# Create file with covariates for PRSice

rule target_trait:
    input:
        info = "data/info/BaseGWAS_TargetTrait_info.csv",
        ids = "../QC/Binary_final/{dataset}_QCed_final.fam",
        comb = empty_file_names
    output:
        trait_prs = "data/TargetTrait/{dataset}_prsinput.{BaseTarget}",
        plot = "data/TargetTrait/{dataset}_prsinput_{BaseTarget}.pdf"
    conda: "./envs/R.yaml"
    script: 
        "scripts/Traits_covariates/{wildcards.dataset}_{wildcards.BaseTarget}_target.R"
# Create file with Target trait for PRSice

rule prsice:
    input:
        multiext("../QC/Binary_final/{dataset}_QCed_final", ".bed", ".bim", ".fam"),
        base_gwas = "data/BaseGWAS_input/gwas_prsice_{BaseGWAS}.gz",
        pheno = "data/TargetTrait/{dataset}_prsinput.{BaseTarget}",
        cov = "data/covariates/{dataset}_{BaseTarget}.covariate"
    output:
        best_fit = "Results/PRSice2/{BaseGWAS}_{BaseTarget}_{dataset}.best"
    params:
        t = "../QC/Binary_final/{dataset}_QCed_final",
        o = "Results/PRSice2/{BaseGWAS}_{BaseTarget}_{dataset}",
        stat = lambda wildcards: df_stat[wildcards.BaseTarget],
        bin_target = lambda wildcards: df_bin_target[wildcards.BaseTarget]
    shell:
        """
        Rscript /home/t-adminsouzamz@lan.afip.com.br/PRSice/PRSice.R \
        --prsice /home/t-adminsouzamz@lan.afip.com.br/PRSice/PRSice_linux \
        --base {input.base_gwas} \
        --target {params.t} \
        --pheno {input.pheno} \
        --cov {input.cov} \
        --stat {params.stat} \
        --binary-target {params.bin_target} \
        --base-maf MAF:0.01 \
        --base-info INFO:0.8 \
        --print-snp \
        --out {params.o}
        """

# ---------------------------------------------- PRS-CS ----------------------------------------------------------------------------------#

rule input_prscs:
    input:
        prsice_gwas = "data/BaseGWAS_input/gwas_prsice_{BaseGWAS}.gz"
    output:
        out = temp("data/BaseGWAS_input/gwas_prscs_{BaseGWAS}.txt")
    conda: "./envs/R.yaml"
    script: 
        "scripts/basegwas_input_prscs.R"
#select and order columns according to PRS-CS 

rule prs_cs:
    input:
        ref = "data/LD_ref/ldblk_1kg_eur",
        gwas = "data/BaseGWAS_input/gwas_prscs_{BaseGWAS}.txt",
        target = "../QC/Binary_final/{dataset}_QCed_final.bim"
    output:
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr1.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr2.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr3.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr4.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr5.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr6.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr7.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr8.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr9.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr10.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr11.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr12.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr13.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr14.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr15.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr16.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr17.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr18.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr19.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr20.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr21.txt"),
        temp("Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr22.txt")
    threads: 1 # Spicy dependency by default uses all available threads
    log:
        "Results/PRSCS/logs/{dataset}_{BaseGWAS}.log"
    params:
        t = "../QC/Binary_final/{dataset}_QCed_final",
        n_gwas_p = lambda wildcards: n_gwas[wildcards.BaseGWAS],
        out = "Results/PRSCS/{dataset}_{BaseGWAS}"
    conda: "./envs/prscs.yaml"
    shell:
        """
        PRScs.py  \
            --ref_dir={input.ref} \
            --bim_prefix={params.t} \
            --sst_file={input.gwas} \
            --n_gwas={params.n_gwas_p} \
            --seed=42 \
            --out_dir={params.out} > {log} 2>&1
        """
# PRS-CS adjusts GWAS effect sizes (LD patterns)
# path to the directory that contains information on the LD reference panel (the snpinfo file and hdf5 files)
# bim_prefix  Full path and the prefix of the bim file for the target dataset. This file is used to provide a list of SNPs that are available in the target dataset.
# sst_file  Full path and the file name of the GWAS summary statistics

rule combine_chr:
    input:
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr1.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr2.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr3.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr4.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr5.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr6.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr7.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr8.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr9.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr10.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr11.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr12.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr13.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr14.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr15.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr16.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr17.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr18.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr19.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr20.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr21.txt",
        "Results/PRSCS/{dataset}_{BaseGWAS}_pst_eff_a1_b0.5_phiauto_chr22.txt"
    output:
        "Results/PRSCS/{dataset}_{BaseGWAS}_allchr.txt" 
    shell:
        """
        cat {input} > {output}
        """
# Combine all chr in one file
rule score:
    input:
        multiext("../QC/Binary_final/{dataset}_QCed_final", ".bed", ".bim", ".fam"),
        score = "Results/PRSCS/{dataset}_{BaseGWAS}_allchr.txt", 
    output:
        "Results/PRSCS/{dataset}_{BaseGWAS}_PRSCS_Score.profile"
    params:
        t = "../QC/Binary_final/{dataset}_QCed_final",
        o = "Results/PRSCS/{dataset}_{BaseGWAS}_PRSCS_Score"
    conda: "./envs/plink.yaml"
    shell:
        """
        plink \
            --bfile {params.t} \
            --score {input.score} 2 4 6 sum center \
            --out {params.o}
        """
# generate polygenic score in plink using adjusted effect values from PRS-CS
# 2 4 6: SNP ID,  effective allele information; adjusted effect size estimate

rule regression:
    input:
        trait = "data/TargetTrait/{dataset}_prsinput.{BaseTarget}",
        covariates = "data/covariates/{dataset}_{BaseTarget}.covariate",
        info = "data/info/BaseGWAS_TargetTrait_info.csv",
        prs_score = "Results/PRSCS/{dataset}_{BaseGWAS}_PRSCS_Score.profile"
    output: 
        coefficients = "Results/PRSCS/{dataset}_{BaseGWAS}_{BaseTarget}_regression_coef.txt",
        r2 = "Results/PRSCS/{dataset}_{BaseGWAS}_{BaseTarget}_regression_r_squared.txt"
    conda: "./envs/R.yaml"
    script:
        "scripts/prscs_regression.R"

