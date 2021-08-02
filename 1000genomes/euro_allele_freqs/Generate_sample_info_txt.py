# This script extracts a sample's genotype information from the merged VCF "EUR_sans_FIN_filtered_keep_w_scores.vcf.gz" that contains information regarding all samples
# Example command to run script: python Generate_sample_info_txt.py --sample NA10865

import vcf
import argparse

# Load arguments from command line
parser = argparse.ArgumentParser()
parser.add_argument("--sample")
    
args = parser.parse_args()

vcf_reader = vcf.Reader(filename='EUR_sans_FIN_filtered_keep_w_scores.vcf.gz')
i = 0
results = open(args.sample + '.txt', "w")
results.write("chrom" + "\t" + "start" + "\t" + "genotype1" + "\t" + "genotype2" + "\n")
for record in vcf_reader:
    i = i + 1

    call = record.genotype(args.sample)
    if call['REPCN'] is None:
        geno1 = 0
        geno2 = 0
    else:
        geno1 = call['REPCN'][0]
        geno2 = call['REPCN'][1]

    results.write(record.CHROM[3:] + '\t' + str(record.POS) + '\t' + str(geno1) + '\t' + str(geno2) + '\n')

results.close()