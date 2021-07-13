import vcf

vcf_reader = vcf.Reader(open('data/sheep_phased_20.vcf', 'r'))

for record in vcf_reader:
        print record
        
        
