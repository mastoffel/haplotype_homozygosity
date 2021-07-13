#!/bin/bash

#$ -V
#$ -cwd
#$ -N shapeit
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH

shapeit4 --input data/sheep.vcf.gz --region 1 --thread 8  --output $SCRATCH/sheep_phased.vcf.gz
rsync -av $SCRATCH/* output/

rm -rf /scratch/$USER/$JOB_ID

