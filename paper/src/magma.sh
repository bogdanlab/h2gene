#!/bin/bash -l
#$ -cwd
#$ -l rh7,h_data=12G,h_rt=20:00:00
#$ -j y
#$ -o ./job_out


# for trait in $(ls out/estimate/); do qsub magma.sh $trait; done
trait=$1

sumstats_dir=/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/sumstats
magma_dir=/u/project/pasaniuc/kangchen/software/magma
# ${magma_dir}/magma --annotate window=10,10 --snp-loc ${magma_dir}/g1000_eur.bim --gene-loc ${magma_dir}/NCBI37.3.gene.loc --out out/magma/step1

python sumstats2pval.py \
    --sumstats ${sumstats_dir}/${trait}.tsv.gz \
    --pfile out/magma/step2/${trait}.pval

${magma_dir}/magma \
    --bfile ${magma_dir}/g1000_eur \
    --pval out/magma/step2/${trait}.pval use='SNP,P' ncol='N' \
    --gene-annot out/magma/step1.genes.annot \
    --out out/magma/step2/${trait}
