#!/bin/bash

#$ -V
#$ -cwd
#$ -N index
#$ -o o_files/
#$ -e e_files/

bgzip -c data/sheep.vcf > data/sheep.vcf.gz
bcftools index data/sheep.vcf.gz

