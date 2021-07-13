#!/bin/bash

#$ -V
#$ -cwd
#$ -N shapeit
#$ -o o_files/
#$ -e e_files/


SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH

for CHR in {1..26}
do
        shapeit4 --input data/sheep.vcf.gz --region $CHR --thread 8 --map data/genetic_map/chr_$CHR.gmap.gz --log $SCRATCH/phased_chr_$CHR.log --output $SCRATCH/sheep_phased_chr_$CHR.vcf.gz
done
        
rsync -av $SCRATCH/* output/

rm -rf /scratch/$USER/$JOB_ID

