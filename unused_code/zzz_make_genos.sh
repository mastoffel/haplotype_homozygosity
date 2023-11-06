#!/bin/bash

#$ -V
#$ -cwd
#$ -N tab_gt_per_chr
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID/
mkdir -p $SCRATCH

#sconda haplotype_homozygosity

Rscript create_gts_per_chr.R $SCRATCH

#conda deactivate
rsync -av $SCRATCH data/genos 

rm -rf /scratch/$USER/$JOB_ID

