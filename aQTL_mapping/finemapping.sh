#!/bin/bash
#
#SBATCH --mem=50G
#SBATCH --array=1-55


TF=$(sed -n $SLURM_ARRAY_TASK_ID\p TFs.txt)
echo -e "TF:${TF} \n"
tissue=$1
echo -e "tissue:${tissue}"


echo -e "\n###################### Identify conditionally independent lead variants  #######################\n"

cd ${dir}/COJO

module load gcta

FILE=${dir}/fastGWA/output/${TF}.${tissue}.MAF005.fastGWA
if [ -f "$FILE" ]; then  
  awk '{ print $2, $4, $5, $7, $8, $9, $10, $6}' ${FILE}  > gwas_summary/${TF}.${tissue}.ma;
  gcta64  --mbfile ${dir}/fastGWA/plink/${tissue}.geno_list.txt  \
    --maf 0.05 \
    --cojo-p 5e-6 \
    --cojo-file gwas_summary/${TF}.${tissue}.ma \
    --cojo-slct \
    --out output/${TF}.${tissue}
  fi
done


echo -e "\n###################### Finemapping  #######################\n"

cd ${dir}/FINEMAP
module load plink/2.0
module load finemap/1.3.1
module load ldstore/1.1


if [ ! -d finemap/${tissue} ]; then mkdir finemap/${tissue}; fi
if [ ! -d bgen/${tissue} ]; then mkdir bgen/${tissue}; fi
if [ ! -d bcor/${tissue} ]; then mkdir bcor/${tissue}; fi


if [ ! -d finemap/${tissue}/${TF} ]; then mkdir finemap/${tissue}/${TF}; fi
if [ ! -d bgen/${tissue}/${TF} ]; then mkdir bgen/${tissue}/${TF}; fi
if [ ! -d bcor/${tissue}/${TF} ]; then mkdir bcor/${tissue}/${TF}; fi


loci=$(sed '1d' ${dir}/COJO/output/${TF}.${tissue}.jma.cojo | awk '{print $2}')


for lead in ${loci};
do

echo -e "\n########## locus: $lead ###########\n"

chr=$(echo $lead | cut -d'_' -f1)
bp=$(echo $lead | cut -d'_' -f2)
prefix=${TF}.${tissue}.${chr}.${bp}


echo -e "\n### Make Z file:\n"

awk -v i=$(echo $chr | sed 's/chr//g') -v hit=$bp \
'NR==1 {print "rsid chromosome position allele1 allele2 maf beta se"; next} { if(($1==i) && ($3>=hit-1e6) && ($3<=hit+1e6)) { print $2,$1,$3,$4,$5,$7<0.5?$7:1-$7,$8,$9 } }' \
${dir}/fastGWA/output/${TF}.${tissue}.MAF005.fastGWA > finemap/${tissue}/${TF}/${prefix}.z


echo -e "\n### Extract genotypes to the bgen file:\n"

awk '{print $1}' finemap/${tissue}/${TF}/${prefix}.z | sed "1d" > bgen/${tissue}/${TF}/${prefix}.snpIDs

plink2 --bfile ${dir}/fastGWA/plink/${tissue}/${tissue}.${chr} \
  --extract bgen/${tissue}/${TF}/${prefix}.snpIDs \
  --export bgen-1.1 \
  --out bgen/${tissue}/${TF}/${prefix}


echo -e "\n### Compute the LD matrix:\n"

ldstore --bgen bgen/${tissue}/${TF}/${prefix}.bgen \
  --bcor bcor/${tissue}/${TF}/${prefix}.bcor  \
  --n-variants-chunk 800\
  --n-threads 2

ldstore --bcor bcor/${tissue}/${TF}/${prefix}.bcor --merge 2

ldstore --bcor bcor/${tissue}/${TF}/${prefix}.bcor --matrix bcor/${tissue}/${TF}/${prefix}.matrix

ldstore --bcor bcor/${tissue}/${TF}/${prefix}.bcor --meta bcor/${tissue}/${TF}/${prefix}.meta

perl -p -e "s/^\s+//g;" bcor/${tissue}/${TF}/${prefix}.matrix | perl -p -e "s/[ ]+/ /g;" > finemap/${tissue}/${TF}/${prefix}.ld


echo -e "\n### Fine-map:\n"

echo "z;ld;snp;config;cred;log;n_samples" > finemap/${tissue}/${TF}/${prefix}
echo -e "finemap/${tissue}/${TF}/${prefix}.z;finemap/${tissue}/${TF}/${prefix}.ld;finemap/${tissue}/${TF}/${prefix}.snp;finemap/${tissue}/${TF}/${prefix}.config;finemap/${tissue}/${TF}/${prefix}.cred;finemap/${tissue}/${TF}/${prefix}.log;$(cat bgen/${tissue}/${TF}/${prefix}.sample | grep "GTEX" | wc -l)" >> finemap/${tissue}/${TF}/${prefix}

finemap --sss --in-files finemap/${tissue}/${TF}/${prefix} --n-causal-snps 1 --log


echo -e "\n### Filter: causal variants\n"

#Signal: the top variant to at least be 5% plausibly causal and have a log10BF > 2
#criteria: keep SNPs with a log10 Bayes factor > 2 as a starting point, 
awk -v p0=$(sed -n 2p finemap/${tissue}/${TF}/${prefix}.snp | cut -d' ' -f11) \
-v lbf0=$(sed -n 2p finemap/${tissue}/${TF}/${prefix}.snp | cut -d' ' -f12) \
'{if (p0>=0.05 && lbf0>=2 && $12>=2) print}' finemap/${tissue}/${TF}/${prefix}.snp > finemap/${tissue}/${TF}/${prefix}.snp.filtered

done


