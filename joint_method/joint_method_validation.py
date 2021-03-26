# Script to validate joint setting of SISTR 

### Imports ###

import sys
sys.path.append("/projects/ps-gymreklab/bonnieh/helper_functions")
from Joint_method_functions import *
PLOTDIR = '/projects/ps-gymreklab/bonnieh/joint_method/'
import argparse
import time
import statistics
from scipy import stats

def getStats(accepted):
    list_means = []
    for elem in accepted:
        list_means.append(elem[0]*elem[1])
    med_value = np.median(list_means)
    lower = np.percentile(list_means, 2.5)
    upper = np.percentile(list_means, 97.5)
    return med_value, lower, upper

# Given file, get ABC lookup table
def getABCTable(file):
    num_bins = 0
    
    # Return list of lists that contain [s, het, common]
    table = GetABCList(file, num_bins)
    dic_summ_stats = {}
        
    # Fill in dic_summ_stats: Key is s, value is list of het, number of common alleles pairs for given s value
    for combo in table:
        s_round = get_LRT_bin(combo[0])
        if s_round not in dic_summ_stats:
            dic_summ_stats[s_round] = []
        dic_summ_stats[s_round].append([combo[1], combo[2]]) # Append het, number of common alleles for s value
    
    return dic_summ_stats
    
def main():
    
    # Load arguments from command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--period", type=int) # Which period to jointly estimate s on
    parser.add_argument("--sim-model", default="eurodem_prior2") # File from which to obtain s values and summary statistics 
    parser.add_argument("--num-sims", type=int, default=2000) # Number of simulations
    parser.add_argument("--num-trials", type=int, default=20) # Number of estimates to obtain
    parser.add_argument("--out-folder", default="results_validation_top_x") # Name of output folder
    parser.add_argument("--out-subfolder") # Name of output subfolder
    parser.add_argument("--a-val", type=float)
    parser.add_argument("--validation-mean", type=float)
    parser.add_argument("--opt-allele-to-validate", type=int) 
    parser.add_argument("--num-loci", type=int, default=1000)
    parser.add_argument("--motif-format", action="store_true")
    parser.add_argument("--header", action="store_true")
    parser.add_argument("--eps-het", type=float, default=0.011) 
    parser.add_argument("--genotyping-errors", action="store_true")
    parser.add_argument("--normalized", action="store_true")
    parser.add_argument("--normalization-constant", type=float, default=0.1)
    parser.add_argument("--ks-test", action="store_true")
    parser.add_argument("--ks-test-threshold", type=float, default=0.05)
    parser.add_argument("--top-x", action="store_true")
    parser.add_argument("--perc-acc", type=float, default=1)
                     
    args = parser.parse_args()
    
    # Naming output file
    outfolder_name = PLOTDIR + args.out_folder + '/' + args.out_subfolder + '/'
    solution_file = open(outfolder_name + 'per_%d_opt_%d_a_%.9f_num_loci_%d_val_mean_%.8f.txt'%(args.period, args.opt_allele_to_validate, args.a_val, args.num_loci, args.validation_mean), 'w') 
    
    solution_file.write("Num sims: " + str(args.num_sims) + '\n')
    b_val = args.validation_mean/args.a_val
    solution_file.write('Number of loci used: ' + str(args.num_loci) + ' a used: ' + str(args.a_val) + ' b used: ' + str(round(b_val,9)) + ' opt allele used: ' + str(args.opt_allele_to_validate) + '\n')
    solution_file.write('Mean of a, b used: ' + str(args.validation_mean) + '\n')
    
    # Preprocess ABC lookup table
    
    opt_allele_dic = {}
    opt_allele_dic[2] = [11,20] #,12,13,14,15,16,17,18,19,
    opt_allele_dic[3] = [5,13] #,6,7,8,9,10,11,12,
    opt_allele_dic[4] = [7,10] #,8,9
    
    # Get ABC tables without genotyping errors
    ABC_tables = {}
   
    for opt_allele in opt_allele_dic[args.period]:
        file = '/projects/ps-gymreklab/bonnieh/abc/results/' + args.sim_model + '/' + str(args.period) + '_' + str(opt_allele) + '.txt' 
        ABC_tables[opt_allele] = getABCTable(file)
        
    # Get ABC tables with genotyping errors
    ABC_tables_w_errors = {}
    
    for opt_allele in opt_allele_dic[args.period]:
        file = '/projects/ps-gymreklab/bonnieh/abc/results/' + args.sim_model + '_errors_low/' + str(args.period) + '_' + str(opt_allele) + '.txt'  
        ABC_tables_w_errors[opt_allele] = getABCTable(file)
       
    # List of medians of mean of accepted a,b pairs
    list_of_medians = []
    
    # List of number of a,b pairs accepted
    num_accepted_pairs = []
    
    # List of tuples with information about a,b inference (median mean, a, b)
    list_med_acc = []
    
    accepted_params_write1 = []
    accepted_params_write2 = []
    accepted_params_write3 = []
    
    # Perform validation
    
    opt_allele_list = [(args.period, args.opt_allele_to_validate)] * args.num_loci
    
    for i in range(0, args.num_trials):
        
        # Get ground truth heterozygosity and number of common alleles distributions 
        obs_het_list, obs_common_list = GetLists(ABC_tables, opt_allele_list, args.a_val, b_val)
        
        if args.genotyping_errors == True:
            obs_het_list, obs_common_list = GetLists(ABC_tables_w_errors, opt_allele_list, args.a_val, b_val)
        
        accepted_params = []
        
        list_all_params = []
        list_all_dist = []
        # Number of simulations for ABC
        for j in range(0, args.num_sims):

            a = np.random.uniform() 
            mu, sigma = np.log(0.0003), np.log(30)
            mean = np.random.lognormal(mu, sigma)
            b = mean/a
            
            if args.ks_test == True:
                het_list, common_list = GetLists(ABC_tables, opt_allele_list, a, b)
                ks_stat, pval = stats.ks_2samp(obs_het_list, het_list)
                
                if args.top_x == True:
                    list_all_params.append((a,b))
                    list_all_dist.append((ks_stat))
                
                else:
                    if pval > args.ks_test_threshold:
                        accepted_params.append((a,b))
                
            else:
                use_common = False
                eps_common = 0.1
                to_accept, mean_of_differences, diff_vector = EstimateParamBinAgnostic(ABC_tables, opt_allele_list, a, b, obs_het_list, obs_common_list, args.eps_het, eps_common, use_common)

                if args.top_x == True:
                    list_all_params.append((a,b))
                    list_all_dist.append((mean_of_differences))
                else:
                    if to_accept == True:
                        accepted_params.append((a, b))
                    
        # Implement top x% 
        if args.top_x == True:
            top_x_num = int(args.perc_acc * args.num_sims / 100)
            top_idx = np.argsort(list_all_dist)[:top_x_num]
            accepted_params = [list_all_params[i] for i in top_idx]
        if i == 0:
            accepted_params_write1 = accepted_params
        if i == 1:
            accepted_params_write2 = accepted_params
        if i == 2:
            accepted_params_write3 = accepted_params
        
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
    solution_file.write('Mean, median of number of accepted pairs: ' + str(np.mean(num_accepted_pairs)) + ',' + str(np.median(num_accepted_pairs)) + ' ' + ';'.join(str(item) for item in num_accepted_pairs) + ' Num accepted out of num trials: ' + str(len(list_of_medians)) + '\n')
    solution_file.write('Mean of medians: ' + str(mean_of_medians) + '\n')
    solution_file.write('Standard deviation: ' + str(stdev) + '\n')
        
    solution_file.write('Accepted medians: ' + ', '.join(str(item) for item in list_med_acc) + '\n')
    solution_file.write('Median of medians: ' + str(np.median(list_of_medians)) + '\n')
    solution_file.write('Lower bound (95 CI) of medians: ' + str(np.percentile(list_of_medians, 2.5)) + '\n')
    solution_file.write('Upper bound (95 CI) of medians: ' + str(np.percentile(list_of_medians, 97.5)) +'\n')
    solution_file.write('Filler line: ' + '\n')
    solution_file.write('Filler line: ' + '\n')
    solution_file.write('List of medians: ' + ','.join(str(item) for item in list_of_medians)+ '\n')
    solution_file.write('Trial 1 list accepted: ' + ','.join(str(item) for item in accepted_params_write1)+ '\n')
    med_value, lower, upper = getStats(accepted_params_write1)
    solution_file.write('Trial 1 stats: ' + str(med_value) + ',' + str(lower) + ',' + str(upper) + '\n')
    solution_file.write('Trial 2 list accepted: ' + ','.join(str(item) for item in accepted_params_write2)+ '\n')
    med_value, lower, upper = getStats(accepted_params_write2)
    solution_file.write('Trial 2 stats: ' + str(med_value) + ',' + str(lower) + ',' + str(upper) + '\n')
    solution_file.write('Trial 3 list accepted: ' + ','.join(str(item) for item in accepted_params_write3)+ '\n')
    med_value, lower, upper = getStats(accepted_params_write3)
    solution_file.write('Trial 3 stats: ' + str(med_value) + ',' + str(lower) + ',' + str(upper) + '\n')
    solution_file.close()
    
if __name__ == '__main__':
    main()