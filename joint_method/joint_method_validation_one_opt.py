# Script to run ABC for all-locus/joint method 

### Imports ###

import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from Joint_method_functions import *
PLOTDIR = '/storage/BonnieH/selection_project/joint_method/results/'
import time
import statistics
# Main function
def main():
    
    period = int(sys.argv[1]) # Which period to jointly estimate s on
    opt_allele = int(sys.argv[2]) # Threshold for optimal allele
    model = sys.argv[3] # File from which to obtain s values and summary statistics 
    
    eps_mean_het = float(sys.argv[4]) # Denominator for epsilon for mean of het distr
    eps_var_het = float(sys.argv[5]) # Denominator for epsilon for variance of het distr
    eps_med_het = float(sys.argv[6]) # Denominator for epsilon for median of het distr
    
    eps_mean_common = float(sys.argv[7]) # Denominator for epsilon for mean of common distr
    eps_var_common = float(sys.argv[8]) # Denominator for epsilon for variance of common distr
    eps_med_common = float(sys.argv[9]) # Denominator for epsilon for median of common distr
    
    num_sims = int(sys.argv[10]) # Number of simulations
    outFolder = sys.argv[11] # Name of outfolder
    use_common_alleles = int(sys.argv[12]) # Whether to use common alleles
    if use_common_alleles == 0:
        use_common_alleles = True
    else:
        use_common_alleles = False
        
    k_val = float(sys.argv[13])
    theta_val = float(sys.argv[14])
    validation_mean = float(sys.argv[15])
    
    
    num_bins = 0
    
    # Naming file
    
    filename = PLOTDIR + 'validation_one_sim/' 
        
    filename = filename + outFolder + '/'
    solution_file = open(filename + 'per_%d_opt_%d_k_%.4f_theta_%.4f_sims_%d_het_eps_%d_%d_%d_comm_eps_%d_%d_%d_val_mean_%.5f.txt'%(period, opt_allele, k_val, theta_val, num_sims, eps_mean_het, eps_var_het, eps_med_het, eps_mean_common, eps_var_common, eps_med_common, validation_mean), 'w') 
    
    solution_file.write("Num sims: " + str(num_sims) + '\n')
   
    
    opt_allele_list = [(period,opt_allele)] * 200
    count = 0
   
    eps_het = [eps_mean_het, eps_var_het, eps_med_het]
    eps_common = [eps_mean_common, eps_var_common, eps_med_common]
    
    # Get ABC tables
    ABC_tables = {}
    opt_allele_dic = {}
    opt_allele_dic[3] = [5,6,7,8,9,10,11,12]
    opt_allele_dic[2] = [11,12,13,14,15,16,17,18,19,20]
    opt_allele_dic[4] = [7,8,9,10]
    
    for opt_allele in opt_allele_dic[period]:
        file = '/gymreklab-tscc/bonnieh/abc/results/'+model+'/' + str(period) + '_' + str(opt_allele) + '.txt' 
        table = GetABCList(file, num_bins)
        dic_summ_stats = {}
        for combo in table:
            s_round = get_LRT_bin(combo[0])
            if s_round not in dic_summ_stats:
                dic_summ_stats[s_round] = []
            dic_summ_stats[s_round].append([combo[1], combo[2]])
                                              
        ABC_tables[opt_allele] = dic_summ_stats
        
    # Perform validation
    obs_het_stats = [1,1,1]
    obs_common_stats = [1,1,1]
    het_list, common_list = EstimateParam(ABC_tables, opt_allele_list, k_val, theta_val, obs_het_stats, \
                                              obs_common_stats, model, eps_het, eps_common, use_common_alleles, True) 
        
        # TODO: Add use_common parameter

    obs_mean_het = np.mean(het_list)
    obs_var_het = np.var(het_list)
    obs_med_het = np.median(het_list)
    obs_het_stats = [obs_mean_het, obs_var_het, obs_med_het]

    obs_mean_common = np.mean(common_list)
    obs_var_common = np.var(common_list)
    obs_med_common = np.median(common_list)
    obs_common_stats = [obs_mean_common, obs_var_common, obs_med_common]

        
    solution_file.write('Number of loci used: ' + str(len(het_list)) + ' k used: ' + str(k_val) + ' theta used: ' + str(theta_val) + '\n')
    solution_file.write('Mean of k, theta used: ' + str(k_val*theta_val) + '\n')
    
    solution_file.write('Mean of observed heterozygosity: ' + str(obs_mean_het) + '\n')
    solution_file.write('Variance of observed heterozygosity: ' + str(obs_var_het) + '\n')
    solution_file.write('Median of observed heterozygosity: ' + str(obs_med_het) + '\n')
    solution_file.write('Mean of observed number of common alleles: ' + str(obs_mean_common) + '\n')
    solution_file.write('Variance of observed number of common alleles: ' + str(obs_var_common) + '\n')
    solution_file.write('Median of observed number of common alleles: ' + str(obs_med_common) + '\n')
       
    accepted_params = []
    total_time_1 = 0
    total_time_2 = 0
    all_time = 0
    for i in range(0, num_sims):
            
        k = np.random.uniform() 
        #prior_type = random.randint(0,1)
        #mu, sigma = np.log(0.002), 1.8
        mu, sigma = np.log(0.0003), np.log(30)
        mean = np.random.lognormal(mu, sigma)
        theta = mean/k
        #t1 = time.time()
        toAdd, het_summ_list, accepted = EstimateParam(ABC_tables, opt_allele_list, k, theta, obs_het_stats, \
                              obs_common_stats, model, eps_het, eps_common, use_common_alleles)
        #t2 = time.time()
        #all_time = all_time + t2-t1
        #print(toAdd)
        #total_time_1  = total_time_1 + time1
        #total_time_2 = total_time_2 + time2
        if toAdd == True:
            accepted_params.append((k, theta))
            
        solution_file.write('k %.5f theta %.5f '%(k, theta) + str(het_summ_list) + ' ' + str(accepted) + '\n')
    solution_file.write('Number k,theta pairs accepted: ' + str(len(accepted_params)) + '\n')
    #k_list = []
    #theta_list = []

    #sort_k = sorted(accepted_params, key=lambda x: x[0])
    #sort_theta = sorted(accepted_params, key=lambda x: x[1])
    sort_mean = sorted(accepted_params, key=lambda x: x[0]*x[1])
    list_of_means = []
    for pair in sort_mean:
        list_of_means.append(pair[0]*pair[1])

    num_accepted = len(accepted_params)
    middle_index = int(num_accepted/2)
    #solution_file.write('Sorted by k - Median k, theta:\n')
    #solution_file.write(str(sort_k[middle_index][0]) + ',' + str(sort_k[middle_index][1]) + '\n')
    #solution_file.write('Sorted by theta - Median k, theta:\n')
    #solution_file.write(str(sort_theta[middle_index][0]) + ',' + str(sort_theta[middle_index][1])+ '\n')
    solution_file.write('Sorted by mean - Median k, theta:\n')
    solution_file.write(str(sort_mean[middle_index][0]) + ',' + str(sort_mean[middle_index][1])+ '\n')
    solution_file.write('Median mean: ' + str(np.median(list_of_means)) + '\n')
    mean_lower_bound = np.percentile(list_of_means, 2.5)
    mean_upper_bound = np.percentile(list_of_means, 97.5)
    solution_file.write('2.5 percentile mean: ' + str(mean_lower_bound) + '\n')
    solution_file.write('97.5 percentile mean: ' + str(mean_upper_bound) + '\n')
    #solution_file.write('Min k, theta: \n')
    #solution_file.write(str(sort_mean[0][0]) + ',' + str(sort_mean[0][1])+ '\n')
    #solution_file.write('Max k, theta: \n')
    #solution_file.write(str(sort_mean[num_accepted-1][0]) + ',' + str(sort_mean[num_accepted-1][1])+ '\n')
        
    #for combo in accepted_params:
        #k = combo[0]
        #theta = combo[1]
        #k_list.append(k)
        #theta_list.append(theta)
            
    #k_med = np.median(k_list)
    #k_lower_bound = np.percentile(k_list, 2.5)
    #k_upper_bound = np.percentile(k_list, 97.5)

    #theta_med = np.median(theta_list)
    #theta_lower_bound = np.percentile(theta_list, 2.5)
    #theta_upper_bound = np.percentile(theta_list, 97.5)
      
    #k_tup = (k_lower_bound, k_med, k_upper_bound)
        
    #theta_tup = (theta_lower_bound, theta_med, theta_upper_bound)
        
    
    #solution_file.write('Estimate of k: ' + ', '.join(str(item) for item in k_tup) + '\n')
    #solution_file.write('Estimate of theta: ' + ', '.join(str(item) for item in theta_tup) + '\n')
    
    solution_file.write('Accepted params: ' + ', '.join(str(item) for item in accepted_params))

    solution_file.close()
    
if __name__ == '__main__':
    main()