### Imports ###

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scipy.stats import geom
import copy
from matplotlib import pyplot as plt
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
    file_type = ''
    if inFile == '/storage/BonnieH/selection_project/ssc_files/allele_freqs/trinuc_100.txt':
        file_type = '_test'
    if inFile == '/storage/BonnieH/selection_project/ssc_files/allele_freqs/trinuc.txt':
        file_type = '_trinuc'
    if inFile == '/storage/BonnieH/selection_project/ssc_files/allele_freqs/allele_freqs_filt.txt':
        file_type = '_all_per'
    filename = str(constant_het) + "_" + str(denom_het) + "_" + str(eps_bins) + \
               "_" + use_het + use_common + use_bins + str(num_bins) + "_" + model + file_type
    outFile = '/storage/BonnieH/selection_project/per_locus/results/' + filename + '.txt'
    figFile = '/storage/BonnieH/selection_project/per_locus/results/figs/' + filename + '.png'
    statsFile = '/storage/BonnieH/selection_project/per_locus/results/stats/' + filename + '.txt'
    
    allele_freqs_file = open(inFile, 'r')
    stats_file = open(statsFile, 'w')
    results = open(outFile, "w")

    # Write results header
    results.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "period" + '\t' + "optimal_ru" + "\t" + "motif" + '\t' + "coding" + \
                  '\t'  + 'intron' + '\t' + 'UTR5' + '\t' + 'UTR3' + '\t' + 'promoter5kb' + '\t' + 'intergenic'  + "\t" + "gene" + \
                  "\t" + "het" + "\t" + "common" + "\t" + "bins" + "\t" + "ABC_s_median" + "\t" + "ABC_s_95%_CI" + "\t" + \
                  "Num_s_accepted" + "\n")

    total_lines_acc = 0
    total_lines = 0
    lower_0 = 0
    s_acc = []
    
    for line in allele_freqs_file:
        # Get information from line
        total_lines = total_lines + 1
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        per = int(info[4])
        motif = info[5]
        coding = info[6]
        intron = info[7]
        UTR5 = info[8]
        UTR3 = info[9]
        promoter5kb = info[10]
        intergenic = info[11]
        gene = info[12]

        opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
        
        #for i in range(0, num_zeros_to_add):
            #allele_freqs.append(0)
            #allele_freqs.insert(0, 0)
        freq_string = ','.join(str(round(item, 5)) for item in allele_freqs)
        if len(allele_freqs) < num_bins:
            num_zeros_to_add = int((num_bins - len(allele_freqs))/2)
            for i in range(0, num_zeros_to_add):
                freq_string = '0.0,' + freq_string
                freq_string = freq_string + ',0.0'
        obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
        results.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(per) + '\t' + str(opt_allele) + '\t' + motif + \
                      '\t' + coding + '\t'  + intron + '\t' + UTR5 + '\t' + UTR3 + '\t' + promoter5kb + '\t' + intergenic + \
                      '\t' + gene + '\t' + str(round(obs_het, 7)) + '\t' + str(obs_common) + '\t' + \
                      ','.join(str(round(item,4)) for item in obs_bins) + '\t')
        
        if per ==3 and opt_allele > 12: 
            opt_allele = 12
        if per ==3 and opt_allele < 5:
            opt_allele = 5
            
        if per ==4 and opt_allele > 10:
            opt_allele = 10
        if per == 4 and opt_allele < 7:
            opt_allele = 7
        
        if per == 2 and opt_allele > 20:
            opt_allele = 20
        if per == 2 and opt_allele < 11:
            opt_allele = 11
            
        abcFile = '/gymreklab-tscc/bonnieh/abc/results/' + model +'/' + str(per) + '_' + str(opt_allele) + '.txt' 
        
        # Read abcFile line by line and place in lookup table in the form of a list
        abc_list = GetABCList(abcFile, num_bins)
        
        # Perform ABC
        s_ABC, lower_bound, upper_bound, num_accepted, s_accepted = Get_S_ABC(abc_list, 
                                       obs_het, obs_common, obs_bins, constant_het, 
                                       denom_het, constant_common, denom_common, eps_bins, use_het, 
                                       use_common, use_bins)
        if s_ABC != -1:
            if s_ABC < 10**-5:
                s_ABC = 0
            else:
                s_ABC = round(s_ABC, 5)

            if lower_bound < 10**-5:
                lower_bound = 0
            else:
                lower_bound = round(lower_bound, 5)

            if upper_bound < 10**-5:
                upper_bound = 0
            else:
                upper_bound = round(upper_bound, 5)

            if lower_bound == 0:
                lower_0 = lower_0 + 1

        if lower_bound != -1:
            total_lines_acc = total_lines_acc + 1
            
        if num_accepted >= 10:
            s_acc.append(num_accepted)
        
        ABC_conf_int = '(' + str(lower_bound) + ' , ' + str(upper_bound) + ')'

        results.write(str(s_ABC) + '\t' + ABC_conf_int + '\t' + str(num_accepted) + '\n')
           
        s_ABC_round = get_LRT_bin(s_ABC)
        
    allele_freqs_file.close()
            
    results.close()
    
    stats_file.write("Total lines: " + str(total_lines) + "\n")
    percent_acc = total_lines_acc/total_lines
    stats_file.write("Total accepted lines: " + str(total_lines_acc) + " " + str(percent_acc) + "\n")
    percent_0 = lower_0/total_lines
    stats_file.write("Lower bound for s is zero: " + str(lower_0) + " " + str(percent_0) + "\n")
    
    for i in range(0, len(s_acc)):
        s_acc[i] = s_acc[i]/10000 * 100
    
    s_acc_med = np.median(s_acc)
    s_acc_mean = np.mean(s_acc)

    stats_file.write("s_acc median: " + str(s_acc_med) + "\n")
    stats_file.write("s_acc mean: " + str(s_acc_mean) + "\n")
 
    plt.figure(1)
    buckets = list(np.arange(0, 101, 1)) 
    
    plt.hist(s_acc, bins=buckets, weights=np.ones(len(s_acc)) / len(s_acc))
    #print(s_acc)

    plt.title('Percent s accepted during ABC %s model \n const_het %.5f, denom_het = %d, num_bins = %d \n Summ stats used: het = %s common alleles = %s bins = %s'%(model, constant_het, denom_het, num_bins, use_het, use_common, use_bins))
    plt.xlabel("% s accepted")
    plt.ylabel("Frequency")
    plt.savefig(figFile, bbox_inches='tight')
    
if __name__ == '__main__':
    main()