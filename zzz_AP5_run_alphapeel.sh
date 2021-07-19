#!/bin/bash

#$ -V
#$ -cwd
#$ -N alphapeel_wi_chr
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID/
mkdir -p $SCRATCH

#sconda haplotype_homozygosity
Rscript phasing.R $SRCATCH
#conda deactivate
rsync -av $SCRATCH output/

rm -rf /scratch/$USER/$JOB_ID

