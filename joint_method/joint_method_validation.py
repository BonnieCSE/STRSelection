# Script to validate ABC for all-locus/joint method 

### Imports ###

import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from Joint_method_functions import *
PLOTDIR = '/storage/BonnieH/selection_project/joint_method/results/'
import time
import statistics

def main():
    
    period = int(sys.argv[1]) # Which period to jointly estimate s on
    model = sys.argv[2] # File from which to obtain s values and summary statistics 
    eps_mean_het = float(sys.argv[3]) # Denominator for epsilon for mean of het distr
    eps_var_het = float(sys.argv[4]) # Denominator for epsilon for variance of het distr
    eps_med_het = float(sys.argv[5]) # Denominator for epsilon for median of het distr
    eps_mean_common = float(sys.argv[6]) # Denominator for epsilon for mean of common distr
    eps_var_common = float(sys.argv[7]) # Denominator for epsilon for variance of common distr
    eps_med_common = float(sys.argv[8]) # Denominator for epsilon for median of common distr
    num_sims = int(sys.argv[9]) # Number of simulations
    num_trials = int(sys.argv[10]) # Number of estimates to obtain
    outFolder = sys.argv[11] # Name of outfolder
    use_common_alleles = int(sys.argv[12]) # Whether to use common alleles
    if use_common_alleles == 0:
        use_common_alleles = True
    else:
        use_common_alleles = False
        
    
    #perform_validation = int(sys.argv[15]) # Whether to perform validation (0 for true)
    #if perform_validation == 0:
        #perform_validation = True
    k_val = float(sys.argv[13])
    theta_val = float(sys.argv[14])
    validation_mean = float(sys.argv[15])
    opt_allele_to_use = int(sys.argv[16])
    num_loci = int(sys.argv[17])
    
    num_bins = 0
    
    # Naming file
    if opt_allele_to_use == 0:
        filename = PLOTDIR + 'validation_multiple/' 
    else:
        filename = PLOTDIR + 'validation_multiple_one_opt/' 
        
    filename = filename + outFolder + '/'
    solution_file = open(filename + 'per_%d_opt_%d_k_%.4f_theta_%.4f_sims_%d_%s_het_eps_%d_%d_%d_comm_eps_%d_%d_%d_val_mean_%.5f.txt'%(period, opt_allele_to_use, k_val, theta_val, num_sims, model, eps_mean_het, eps_var_het, eps_med_het, eps_mean_common, eps_var_common, eps_med_common, validation_mean), 'w') 
    
    solution_file.write("Num sims: " + str(num_sims) + '\n')
    solution_file.write('Number of loci used: ' + str(num_loci) + ' k used: ' + str(k_val) + ' theta used: ' + str(theta_val) + ' opt allele used: ' + str(opt_allele_to_use) + '\n')
    solution_file.write('Mean of k, theta used: ' + str(k_val*theta_val) + '\n')
    
    eps_het = [eps_mean_het, eps_var_het, eps_med_het]
    eps_common = [eps_mean_common, eps_var_common, eps_med_common]
    
    # Get ABC tables
    ABC_tables = {}
    opt_allele_dic = {}
    opt_allele_dic[3] = [5,6,7,8,9,10,11,12,13]
    opt_allele_dic[2] = [11,12,13,14,15,16,17,18,19,20]
    opt_allele_dic[4] = [7,8,9,10]
    
    opt_allele_dic_w_per = {}
    opt_allele_dic_w_per[3] = [(3,5),(3,6),(3,7),(3,8),(3,9),(3,10),(3,11),(3,12),(3,13)]
    opt_allele_dic_w_per[2] = [(2,11),(2,12),(2,13),(2,14),(2,15),(2,16),(2,17),(2,18),(2,19),(2,20)]
    opt_allele_dic_w_per[4] = [(4,7),(4,8),(4,9),(4,10)]
    
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
        
    list_of_medians = []
    num_accepted_pairs = []
    list_med_acc = []
    
    gt_mean = []
    gt_var = []
    gt_med = []
    # Perform validation
    if opt_allele_to_use == 0:
        opt_allele_sub_list = random.choices(opt_allele_dic_w_per[period], k=num_loci)
    else:
        opt_allele_sub_list = [(period, opt_allele_to_use)] * num_loci
    for i in range(0, num_trials):
        #opt_allele_sub_list = random.sample(opt_allele_list, 1000)
        obs_het_stats = [1,1,1]
        obs_common_stats = [1,1,1]
        
        
        het_list, common_list = EstimateParam(ABC_tables, opt_allele_sub_list, k_val, theta_val, obs_het_stats, \
                                              obs_common_stats, model, eps_het, eps_common, use_common_alleles, True) 
        
        # TODO: Add use_common parameter

        obs_mean_het = np.mean(het_list)
        obs_var_het = np.var(het_list)
        obs_med_het = np.median(het_list)
        obs_het_stats = [obs_mean_het, obs_var_het, obs_med_het]
        gt_mean.append(obs_mean_het)
        gt_var.append(obs_var_het)
        gt_med.append(obs_med_het)
        obs_mean_common = np.mean(common_list)
        obs_var_common = np.var(common_list)
        obs_med_common = np.median(common_list)
        obs_common_stats = [obs_mean_common, obs_var_common, obs_med_common]
       
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
            t1 = time.time()
            toAdd, time1, time2 = EstimateParam(ABC_tables, opt_allele_sub_list, k, theta, obs_het_stats, \
                                  obs_common_stats, model, eps_het, eps_common, use_common_alleles)
            t2 = time.time()
            all_time = all_time + t2-t1
            #print(toAdd)
            total_time_1  = total_time_1 + time1
            total_time_2 = total_time_2 + time2
            if toAdd == True:
                accepted_params.append((k, theta))

        
        #solution_file.write('Number k,theta pairs accepted: ' + str(len(accepted_params)) + '\n')
        #k_list = []
        #theta_list = []

        #sort_k = sorted(accepted_params, key=lambda x: x[0])
        #sort_theta = sorted(accepted_params, key=lambda x: x[1])
        sort_mean = sorted(accepted_params, key=lambda x: x[0]*x[1])
        
        num_accepted = len(accepted_params)
        num_accepted_pairs.append(num_accepted)
        middle_index = int(num_accepted/2)
        if num_accepted > 0:
            median_mean = sort_mean[middle_index][0] * sort_mean[middle_index][1]
            list_of_medians.append(median_mean)
            list_med_acc.append((median_mean, sort_mean[middle_index][0], sort_mean[middle_index][1]))
        
    mean_of_medians = np.mean(list_of_medians)
    stdev = np.std(list_of_medians)
    solution_file.write('Mean, median of number of accepted pairs: ' + str(np.mean(num_accepted_pairs)) + ',' + str(np.median(num_accepted_pairs)) + '\n')
    solution_file.write('Mean of medians: ' + str(mean_of_medians) + '\n')
    solution_file.write('Standard deviation: ' + str(stdev) + '\n')
        

    solution_file.write('Accepted medians: ' + ', '.join(str(item) for item in list_med_acc) + '\n') #list_of_medians
    solution_file.write('Ground truth means: ' + ','.join(str(item) for item in gt_mean)+ '\n')
    solution_file.write('Ground truth vars: ' + ','.join(str(item) for item in gt_var)+ '\n')
    solution_file.write('Ground truth meds: ' + ','.join(str(item) for item in gt_med)+ '\n')
    solution_file.close()
    
if __name__ == '__main__':
    main()