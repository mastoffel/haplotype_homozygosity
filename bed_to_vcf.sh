#!/bin/bash

#$ -V
#$ -cwd
#$ -N bed_to_vcf
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID/
mkdir -p $SCRATCH

plink --bfile data/plink/sheep --recode vcf --chr 1-26 --sheep --snps-only just-acgt --out $SCRATCH/sheep

rsync -av $SCRATCH data/

rm -rf /scratch/$USER/$JOB_ID

