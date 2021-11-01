#!/bin/bash

for i in {1..26}
do
   bcftools convert --hapsample --vcf-ids output/phased/sheep_phased_chr_$i.vcf.gz -o output/phased_matrix/sheep_phased_chr_$i
done