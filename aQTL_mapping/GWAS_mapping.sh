#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --array=1-49


cd ${dir}/fastgwa

tissue=$(sed -n $SLURM_ARRAY_TASK_ID\p tissues.txt)
echo -e "\n${tissue}\n"


echo -e "\n###################### Conver genotype vcf files to PLINK format (each tissue) #######################\n"

module load plink/2.0

vcf_file=$(echo "vcf.files/GTEx_v8_WholeGenomeSeq_${tissue}.vcf.gz")

mkdir -p plink/${tissue}

for i in $(seq 1 22);
do 
plink2 --vcf ${vcf_file} \
  --make-bed \
  --vcf-half-call 'missing' \
  --chr ${i} \
  --out plink/${tissue}/${tissue}.chr${i}
done

for i in $(seq 1 22); do echo "plink/${tissue}/${tissue}.chr${i}"; done > plink/${tissue}.geno_list.txt 



echo -e "\n###################### Estimate GRM (genetic relationship matrix) among individuals #######################\n"

gcta64 --mbfile plink/${tissue}.geno_list.txt \
  --make-grm \
  --sparse-cutoff 0.05 \
  --thread-num 8 \
  --out sparse_grm/${tissue}.auto



echo -e "\n###################### Run fastGWA #######################\n"

module load gcta

for TF in $(cat TFs.txt);
do
  gcta64  --mbfile plink/${tissue}.geno_list.txt \
    --maf 0.05 \
    --fastGWA-mlm \
    --grm-sparse sparse_grm/${tissue}.auto \
    --pheno phenotypes/${TF}.${tissue}.txt \
    --covar covariates/covar.min.${tissue}.txt \
    --qcovar covariates/qcovar.min.${tissue}.txt \
    --thread-num 8 \
    --out output/${TF}.${tissue}.MAF005;
done


