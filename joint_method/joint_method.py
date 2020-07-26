# Script to run ABC for all-locus/joint method 

### Imports ###

import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from Joint_method_functions import *
PLOTDIR = '/storage/BonnieH/selection_project/joint_method/results/'

# Main function
def main():
    inFile = '/storage/BonnieH/selection_project/ssc_files/allele_freqs/allele_freqs_filt.txt'
    allele_freqs_file = open(inFile, 'r')
    header = allele_freqs_file.readline().strip()
    column = int(sys.argv[1]) # Which column of file to jointly estimate s on
    mot = sys.argv[2] # Which motif to jointly estimate s on
    period = int(sys.argv[3])
    opt_thresh = int(sys.argv[4]) # Threshold for optimal allele
    model = sys.argv[5] # First file from which to obtain s values and summary statistics 
    eps_mean = int(sys.argv[6]) # Denominator for epsilon for mean
    eps_var = int(sys.argv[7]) # Denominator for epsilon for variance
    eps_med = int(sys.argv[8]) # Denominator for epsilon for median
    num_sims = int(sys.argv[9]) # Number of simulations
    outFolder = sys.argv[10] # Name of outfolder
    perform_validation = int(sys.argv[11]) # Whether to perform validation
    if perform_validation == 0:
        perform_validation = True
    else:
        perform_validation = False
        
    ### Validation inputs ###
    k_val = float(sys.argv[12])
    theta_val = float(sys.argv[13])
    num_bins = int(sys.argv[14])
    
    # Naming file
    filename = PLOTDIR + 'no_validation/'
    if perform_validation == True:
        filename = PLOTDIR + 'validation/' 
        
    filename = filename + outFolder + '/'
    solution_file = open(filename + 'per' + str(period) + '_' + str(column) + '_' + mot + '_thresh_' + \
                         str(opt_thresh) + '_k_' + str(k_val) + '_theta_' + str(theta_val) + '_sims_' + \
                         str(num_sims) + '_' + str(eps_mean) + str(eps_var) + str(eps_med) + '.txt', 'w')
    
    solution_file.write("Num sims: " + str(num_sims) + '\n')
    
    motif_dic = {} # Dictionary of motifs; Key motif; Value number of motifs
    
    # Get all motifs
    for line in allele_freqs_file:
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        per = int(info[4])
        motif = info[column]
        opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
        if opt_allele >= opt_thresh and per == period:
            if motif not in motif_dic:
                motif_dic[motif] = 1
            else:
                motif_dic[motif] = motif_dic[motif] + 1
              
    allele_freqs_file.close()
    
    # Only look at ~200 loci in the STR class -> faster runtime
    divisor = int(motif_dic[mot]/200)+1
       
    allele_freqs_file = open(inFile, 'r')
    header = allele_freqs_file.readline().strip()
    
    obs_het_distr = []
    obs_common_distr = []
    opt_allele_list = []
    count = 0
   
    # Get STRs (represented by optimal allele list) to estimate s jointly on
    # Get observed het distribution on these STRs
    for line in allele_freqs_file:
        # Get information from line
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        per = int(info[4])
        
        motif = info[column]
            
        if motif == mot and per == period:
            opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
            freq_string = ','.join(str(round(item, 5)) for item in allele_freqs)
            if len(allele_freqs) < num_bins:
                num_zeros_to_add = int((num_bins - len(allele_freqs))/2)
                for i in range(0, num_zeros_to_add):
                    freq_string = '0.0,' + freq_string
                    freq_string = freq_string + ',0.0'
            obs_het, obs_common, obs_bins = GetSummStats(freq_string, num_bins)
            if opt_allele >= opt_thresh:
                count = count + 1
                if count%divisor==0: 
                    obs_het_distr.append(obs_het)
                    obs_common_distr.append(obs_common)
                    opt_allele_list.append((per, opt_allele))
    
    obs_mean_het = np.mean(obs_het_distr)
    obs_var_het = np.var(obs_het_distr)
    obs_med_het = np.median(obs_het_distr)
    obs_het_stats = [obs_mean_het, obs_var_het, obs_med_het]
    
    obs_mean_common = np.mean(obs_common_distr)
    obs_var_common = np.var(obs_common_distr)
    obs_med_common = np.median(obs_common_distr)
    obs_common_stats = [obs_mean_common, obs_var_common, obs_med_common]
    eps_het = 
    eps_common = 
    if perform_validation == True:
        toAdd, het_list, common_list = EstimateParam(opt_allele_list, k_val, theta_val, obs_het_stats, obs_common_stats, \
                                        model, eps_het, eps_common) 
        
        # TODO: Add use_common parameter

        obs_mean = np.mean(het_list)
        obs_var = np.var(het_list)
        obs_med = np.median(het_list)
        
    solution_file.write('Number of loci used: ' + str(len(obs_het_distr)) + ' k used: ' + str(k_val) + ' theta used: ' + str(theta_val) + '\n')
    solution_file.write('Column: ' + str(column) +  ' Value of column: ' + mot + '\n')
    solution_file.write('Mean of observed heterozygosity: ' + str(obs_mean) + '\n')
    solution_file.write('Variance of observed heterozygosity: ' + str(obs_var) + '\n')
    solution_file.write('Median of observed heterozygosity: ' + str(obs_med) + '\n')
       
    accepted_params = []
    
    for i in range(0, num_sims):
            
        k = np.random.uniform()
        prior_type = random.randint(0,1)
        mu, sigma = np.log(0.002), 1.8
        mean = np.random.lognormal(mu, sigma)
        theta = mean/k

        toAdd, het_list = EstimateParam(opt_allele_list, k, theta, obs_mean, obs_var, obs_vec, obs_med, \
                                        model, eps_mean, eps_var, eps_med)
        if toAdd == True:
            accepted_params.append((k, theta))
            
    k_list = []
    theta_list = []

    sort_k = sorted(accepted_params, key=lambda x: x[0])
    sort_theta = sorted(accepted_params, key=lambda x: x[1])
    sort_mean = sorted(accepted_params, key=lambda x: x[0]*x[1])

    num_accepted = len(accepted_params)
    middle_index = int(num_accepted/2)
    solution_file.write('Sorted by k - Median k, theta:\n')
    solution_file.write(str(sort_k[middle_index][0]) + ',' + str(sort_k[middle_index][1]) + '\n')
    solution_file.write('Sorted by theta - Median k, theta:\n')
    solution_file.write(str(sort_theta[middle_index][0]) + ',' + str(sort_theta[middle_index][1])+ '\n')
    solution_file.write('Sorted by mean - Median k, theta:\n')
    solution_file.write(str(sort_mean[middle_index][0]) + ',' + str(sort_mean[middle_index][1])+ '\n')
    solution_file.write('Min k, theta: \n')
    solution_file.write(str(sort_mean[0][0]) + ',' + str(sort_mean[0][1])+ '\n')
    solution_file.write('Max k, theta: \n')
    solution_file.write(str(sort_mean[num_accepted-1][0]) + ',' + str(sort_mean[num_accepted-1][1])+ '\n')
        
    for combo in accepted_params:
        k = combo[0]
        theta = combo[1]
        k_list.append(k)
        theta_list.append(theta)
            
    k_med = np.median(k_list)
    k_lower_bound = np.percentile(k_list, 2.5)
    k_upper_bound = np.percentile(k_list, 97.5)

    theta_med = np.median(theta_list)
    theta_lower_bound = np.percentile(theta_list, 2.5)
    theta_upper_bound = np.percentile(theta_list, 97.5)
      
    k_tup = (k_lower_bound, k_med, k_upper_bound)
        
    theta_tup = (theta_lower_bound, theta_med, theta_upper_bound)
        
    allele_freqs_file.close()
    solution_file.write('Estimate of k: ' + ', '.join(str(item) for item in k_tup) + '\n')
    solution_file.write('Estimate of theta: ' + ', '.join(str(item) for item in theta_tup) + '\n')
    solution_file.write('Number k,theta pairs accepted: ' + str(len(accepted_params)) + '\n')
    solution_file.write('Accepted params: ' + ', '.join(str(item) for item in accepted_params))

    solution_file.close()
    
if __name__ == '__main__':
    main()