# Mapping genetic determinants of TF activity (aQTLs)

## Input data:
1) Phenotype files: Inferred TF activity for each tissue and each TF from the first step.
2) Genotype files: [GTEx Protected Access Data](https://gtexportal.org/home/protectedDataAccess)
3) Covariates: Main technical covariates are extracted from [GTEx subject phenotypes ](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt)

## Scripts:
1) Genome-wide associate studies: `GWAS_mapping.sh`   
2) Finemapping: `finemapping.sh`   

## Dependencies:
- plink v2.0   
- gcta v1.93.1   
- ldstore v1.1   
- finemap v1.3.1   