'''Script to run LRT for per-locus method to get p value testing whether a model 
of selection fits better than a model without selection for each STR 
'''

### Imports ###

import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from LRT_functions import *

### Main function ###
def main():
    # Load parameters
    constant_het = float(sys.argv[1])
    denom_het = int(sys.argv[2])
    constant_common = int(sys.argv[3])
    denom_common = int(sys.argv[4])
    eps_bins = float(sys.argv[5])
    inFile = sys.argv[6]
    use_het = sys.argv[7]
    use_common = sys.argv[8]
    use_bins = sys.argv[9]
    num_bins = int(sys.argv[10])
    model = sys.argv[11]
    LRT_num_sims = int(sys.argv[12])
    
    filename = str(constant_het) + "_" + str(denom_het) + "_" + str(eps_bins) + \
               "_" + use_het + use_common + use_bins + str(num_bins) + "_" + model 
    
    outFile = '/storage/BonnieH/selection_project/per_locus/final_results/' + filename + '.txt'
    
    allele_freqs_file = open(inFile, 'r')
    header = allele_freqs_file.readline().strip()
    results = open(outFile, "w")

    # Write results header
    results.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "period" + "\t" +"optimal_ru" + '\t'+ \
                  'motif' + "\t" + "coding" + "\t" + 'intron' + '\t' + 'UTR5' + '\t' + 'UTR3' + '\t' + \
                  'promoter5kb' + '\t' + 'intergenic'  + "\t" + "gene" + "\t" + "het" + "\t" + "common" + \
                  "\t" + "bins" + "\t" + "ABC_s_median" + "\t" + "ABC_s_95%_CI" + "\t" + "Num_s_accepted" + \
                  "\t" + "Likelihood_0" + "\t" + "Likelihood_s" + "\t" + "LR" + "\t" + "LogLR" + "\t" + "LRT_p_value" + "\n")
   
    for line in allele_freqs_file:
        # Get information from line 
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        per = int(info[3])
        optimal_ru = int(info[4])
        motif = info[5]
        coding = info[6]
        intron = info[7]
        UTR5 = info[8]
        UTR3 = info[9]
        promoter5kb = info[10]
        intergenic = info[11]
        gene = info[12]
        obs_het = float(info[13])
        obs_common = int(info[14])
        obs_bins = info[15]
        s_ABC = float(info[16])
        ABC_conf_int = info[17]
        num_accepted = int(info[18])
         
        results.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(per) + '\t' + str(optimal_ru) + \
                      '\t' + motif + '\t' + coding + '\t'  + intron + '\t' + UTR5 + '\t' + UTR3 + '\t' + \
                      promoter5kb + '\t' + intergenic + '\t' + gene + '\t' + str(round(obs_het, 7)) + '\t' + \
                      str(obs_common) + '\t' + obs_bins + '\t')
        
        # Write N/A if <10 s accepted during ABC
        if s_ABC == -1:
            s_ABC = 'N/A'
            ABC_conf_int = 'N/A'
            num_accepted = '<10'
            likelihood_0 = 'N/A'
            likelihood_s_ABC = 'N/A'
            LR = 'N/A'
            LogLR = 'N/A'
            pval = 'N/A'
        
            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(num_accepted) + '\t')
            results.write(str(likelihood_0) + '\t' + str(likelihood_s_ABC) + '\t' + \
                          str(LR) + '\t' + str(LogLR) + '\t' + str(pval) + '\n')
         
        # Perform LRT if can get posterior estimate of s
        else:
            results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(num_accepted) + '\t')
        
            if per == 3 and optimal_ru > 12:
                optimal_ru = 12
            if per == 3 and optimal_ru < 5:
                optimal_ru = 5

            if per == 4 and optimal_ru > 10:
                optimal_ru = 10
            if per == 4 and optimal_ru < 7:
                optimal_ru = 7

            if per == 2 and optimal_ru > 20:
                optimal_ru = 20
            if per == 2 and optimal_ru < 11:
                optimal_ru = 11

            s_ABC_round = get_LRT_bin(s_ABC)

            # Perform LRT

            lrtFile = '/gymreklab-tscc/bonnieh/lrt/results/' + model + '_prelim' + '/' + str(per) + '_' + str(optimal_ru) + '_freqs.txt'

            ### Get list of s values in LRT simulations file ###
            s_list_available = []
            lrtFile = '/gymreklab-tscc/bonnieh/lrt/results/' + 'euro_prelim' + '/' + str(per) + '_' + str(optimal_ru) + '_freqs.txt'
            lrt_file = open(lrtFile, 'r')
            header = lrt_file.readline().strip()

            for line in lrt_file:
                info = line.strip().split('\t')
                s = float(info[0])
                if s not in s_list_available:
                    s_list_available.append(s)

            # Get nearest s in LRT file to perform LRT
            if s_ABC_round not in s_list_available:
                s_ABC_round = getNearestS(s_ABC_round, s_list_available)

            # Get LRT summary statistic tables for s = 0
            freqs_list_raw_0 = GetLRTListFreq(lrtFile, 0)
            LRT_table_0_het = []
            LRT_table_0_common = []
            LRT_table_0_bins = []

            # Get summary statistics from allele frequencies
            for freq_string in freqs_list_raw_0:

                obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
                LRT_table_0_het.append(obs_het) 
                LRT_table_0_common.append(obs_common) 
                LRT_table_0_bins.append(obs_bins)

            # Get LRT summary statistic tables for s = s_ABC_round
            freqs_list_raw_s = GetLRTListFreq(lrtFile, s_ABC_round)

            LRT_table_s_het = []
            LRT_table_s_common = []
            LRT_table_s_bins = []

            # Get summary statistics from allele frequencies
            for freq_string in freqs_list_raw_s:

                obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
                LRT_table_s_het.append(obs_het) 
                LRT_table_s_common.append(obs_common) 
                LRT_table_s_bins.append(obs_bins)

            # Perform LRT
            likelihood_0, likelihood_s_ABC, LR, LogLR, pval = LikelihoodRatioTest(LRT_table_0_het, LRT_table_s_het, \
                                    LRT_table_0_common, LRT_table_s_common, LRT_table_0_bins, LRT_table_s_bins, LRT_num_sims, \
                                    obs_het, obs_common, obs_bins, constant_het, denom_het, \
                                    constant_common, denom_common, eps_bins, use_het, use_common, use_bins)

            results.write(str(round(likelihood_0, 7)) + '\t' + str(round(likelihood_s_ABC, 7)) + '\t' + \
                              str(round(LR, 7)) + '\t' + str(round(LogLR ,7)) + '\t' + str(round(pval, 7)) + '\n')

    allele_freqs_file.close()
            
    results.close()

if __name__ == '__main__':
    main()