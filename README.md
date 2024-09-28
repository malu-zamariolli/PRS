# Polygenic Risk Scores (PRS)

# Usage

This Github contains the workflow pipeline that was used to run polygenic risk scores with two approaches: [PRS-CS](https://github.com/getian107/PRScs) and [PRSice2](https://choishingwan.github.io/PRSice/). The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

To run the pipeline with snakemake:
```console
cd  PRS/
mamba activate snakemake
snakemake --use-conda --cores 1 -s snakefile
# number of cores can be modified accordingly
```

# Input data
For this study, we need:
- GWAS summary statistics (Base GWAS data) - PRS calculation is based on this data
- Plink binary files for the samples in which PRS will be applied (Target data)
- Information file describing combinations of Base GWAS data and Target data to be applied, as well as other parameters
- File with traits/phenotypes of interest in the target sample
- File with covariates of the target sample to be considered

## Base GWAS data
GWAS summary statistics comes from various studies. Each GWAS is processed by a individual script (rule input:gwas in snakemake) to obtain the following columns:

**CHR, BP, SNP, A1, A2, P, INFO,  OR / BETA, MAF**

- columns are white space delimited
- When MAF and INFO columns are not present, they are created with dummy data to not interfere in the following steps
- When CHR and BP columns are not present, dummy data is created to not interfere in the following steps. In this case, columns should not be named CHR and BP. (If CHR and BP are present in the column header, PRSice2 will identify them and fail) 


Files should be named in order to be recognized by defined wildcards in snakemake:

*gwas_{BaseGWAS}.{ext} -> gwas_DepressionHoward.gz*

## Target data
Plink binary files (.bed, .bim, .fam) after quality control steps.

Files should be named in order to be recognized by defined wildcards in snakemake:

*{dataset}.bed -> EUR.bed*

*{dataset}.bim -> EUR.bim*

*{dataset}.fam -> EUR.fam*


## Information file (.csv)
A .csv file with information regarding the study to be performed. Some columns are used as parameters in the pipeline

|BaseTarget|BaseGWAS|N_GWAS|TargetTrait|TraitColumn|Covariates|BinaryTargetOption|STAT|
| -------- | --- | --- | --- | --- | --- | --- | --- |
|DepHeight|DepressionHoward|807553|Height|Altura|Age-age^2-PC1-10-sex|F|OR|

- BaseTarget: Unique identifier to each study to be performed. Ensures the right combination of BaseGWAS data and Target Trait
- N_GWAS: Number of individuals in the base GWAS. It is a parameter for PRS-CS tool
- TraitColumn: Name of the column with the phenotype of interest
- Covariates: covariates to be considered
- BinaryTargetOption: If T, target trait is binary. If F, target trait is continuous. It's used as a parameter in the PRSice2 tool
- STAT: if OR, base GWAS reports OR (binary). If BETA, base GWAS reports beta (continuous). It's used as a parameter in the PRSice2 tool

## Traits/Phenotypes
A file with the target trait or phenotype to be investigated in the target sample is generated in the analysis. Since the source of the trait may vary, a script per study needs to be written and will be called in the snakemake pipeline with wildcards

*{dataset}_{BaseTarget}_target.R -> EUR_DepHeight_target.R*

The resulting file should have the following columns **(without header)** and delimiter should be **white space**

**FID, IID, TargetTrait**

## Covariates
A file with the covariates from the target sample to be considered is generated in the analysis. Since the source of the covariates may vary, a script per study needs to be written and will be called in the snakemake pipeline with wildcards

*{dataset}_{BaseTarget}_covariates.R -> EUR_DepHeight_covariates.R*

The resulting file should start with the columns **FID	IID** followed by columns with covariates **(with header)** and delimiter should be **tab**

# Tools

- [PRS-CS](https://github.com/getian107/PRScs), [Plinkv1.9](https://www.cog-genomics.org/plink/) and R packages can be used from environments created in the snakemake pipeline
    ../envs/plink.yaml
    ../envs/prscs.yaml
    ../envs/R.yaml
- [PRSice2](https://choishingwan.github.io/PRSice/) needs to be [installed](https://github.com/choishingwan/PRSice/) or binary needs to be downloaded. 

# References
Some steps in this pipeline were adapted from [PRS-tutorial](https://choishingwan.github.io/PRS-Tutorial/prsice/).
