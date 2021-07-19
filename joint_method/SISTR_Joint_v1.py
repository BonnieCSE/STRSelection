# Script to run SISTR Joint Method Version 1

### Imports ###

import argparse
import math
import random
from scipy import stats
from Joint_method_functions import *

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
    
### Main function ###
def main():
    
    # Load arguments from command line
    # See documentation for explanation of each option
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--in-file") 
    parser.add_argument("--out-folder") 
    parser.add_argument("--lookup-table-folder", default = 'sistr_joint_resources/') 
    parser.add_argument("--period", type=int) 
    parser.add_argument("--sim-model") 
    parser.add_argument("--motif-to-use")
    parser.add_argument("--opt-allele-to-use", type=int) 
    parser.add_argument("--num-sims", type=int, default=50000) 
    parser.add_argument("--num-loci", type=int, default=1000)
    parser.add_argument("--motif-format", action="store_true")
    parser.add_argument("--header", action="store_true")
    parser.add_argument("--column", type=int, default=6)
    
    parser.add_argument("--normalized", action="store_true") # Old setting: Whether to use normalized values when comparing two heterozygosity distributions  
    parser.add_argument("--normalization-constant", type=float, default=0.1) # Old setting: Normalization constant if using normalized values when comparing 2 heterozygosity distributions (default value: 0.1)   
    parser.add_argument("--ks-test", action="store_true") # Old setting: Whether to use Kolmogorov-Smirnov test to compare two heterozygosity distributions (default option is taking mean of differences of the two distributions)
    
    # Parameters for accepting top x % method
    parser.add_argument("--top-x", action="store_true") # Whether to use the top x percent method (default option is using arbitrary thresholds to decide whether two heterozyosity distributions are sufficiently similar)
    parser.add_argument("--perc-acc", type=float, default=1)
                     
    args = parser.parse_args()
    
    inFile = args.in_file
    
    # Naming output file
    opt_allele_string = str(args.opt_allele_to_use)
    
    # Add 0 in front of optimal allele if < 10. This allows the files to be named in order from smallest to largest optimal allele if there are multiple files.
    if args.opt_allele_to_use < 10:
        opt_allele_string = '0' + opt_allele_string
    
    solution_file = open(args.out_folder + 'per_%d_%s_%s_col_%d_opt_%s.txt'%(args.period, args.motif_to_use, args.sim_model, args.column, opt_allele_string), 'w') 
    
    solution_file.write("Num sims: " + str(args.num_sims) + '\n')
    
    # Preprocess ABC lookup table
    # Get ABC tables
    ABC_tables = {}
    ABC_tables_master_dic = {}
    
    mut_setting_folder_name = {}
    mut_setting_folder_name[2] = 'eurodem_prior2_dinuc_'
    mut_setting_folder_name[3] = 'eurodem_prior2_trinuc_'
    mut_setting_folder_name[4] = 'eurodem_prior2_tetranuc_'
    
    mut_settings = {}
    mut_settings[2] = ['a','b','c','d','e','f']
    mut_settings[3] = ['a','b','c','d','e','f','g']
    mut_settings[4] = ['a','b','c','d','e','f','g']
    
    if args.sim_model == 'all':
        
        folder_prefix = mut_setting_folder_name[args.period]
        
        for mut_setting in mut_settings[args.period]:
            file = args.lookup_table_folder + folder_prefix + mut_setting + '_1kg_euro/' + str(args.period) + '_' + str(args.opt_allele_to_use) + '.txt'

            ABC_dic = {}
            ABC_dic[args.opt_allele_to_use] = getABCTable(file)
            ABC_tables_master_dic[mut_setting] = ABC_dic
            
    else:
        
        file = args.lookup_table_folder + args.sim_model + '/' + str(args.period) + '_' + str(args.opt_allele_to_use) + '.txt'

        ABC_tables[args.opt_allele_to_use] = getABCTable(file)
        
    # Open input file
    allele_freqs_file = open(inFile, 'r')
    
    obs_het_distr_all = []
    obs_common_distr_all = []
    opt_allele_list = []
    count = 0 # Count number of STRs in the class of loci specified
   
    # Get STRs (represented by optimal allele list) to estimate s jointly on
    # Get observed het and common alleles distribution on these STRs
    for line in allele_freqs_file:
        
        # Get information from line
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        per = int(info[5])
        motif = info[args.column]
        
        if args.column == 6:
            canonical_motif = GetCanonicalRU(motif)
        else:
            canonical_motif = motif
            
        if args.motif_format == True:
            opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start, True, True)
        else:
            opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start, True)
        
        if (canonical_motif == args.motif_to_use or args.motif_to_use == 'all') and per == args.period and args.opt_allele_to_use == opt_allele:
            count = count + 1
            obs_het = 1-sum([item**2 for item in allele_freqs])
            obs_common = len([i for i in allele_freqs if i >= 0.05])
            obs_het_distr_all.append(obs_het)
            obs_common_distr_all.append(obs_common)
            opt_allele_list.append((per, opt_allele))
                   
    allele_freqs_file.close()
   
    # Get a subset of loci if necessary (number of loci > num-loci option)
    opt_allele_sub_list = []
    obs_het_distr = []
    obs_common_distr = []
    total_number_loci = len(opt_allele_list)
    
    # Use all loci since number of loci <= num-loci option
    if len(opt_allele_list) <= args.num_loci:
        opt_allele_sub_list = opt_allele_list
        obs_het_distr = obs_het_distr_all
        obs_common_distr = obs_common_distr_all
      
    # Get heterozygosity and number of common alleles distributions for random subset of num-loci number of loci
    else:
        
        list_random_numbers = random.sample(range(len(opt_allele_list)), args.num_loci)
        for index in list_random_numbers:
            opt_allele_sub_list.append(opt_allele_list[index])
            obs_het_distr.append(obs_het_distr_all[index])
            obs_common_distr.append(obs_common_distr_all[index])
   
    solution_file.write('Number of loci used: ' + str(len(obs_het_distr)) + ' Total number loci: ' + str(total_number_loci) + '\n')
    solution_file.write('Column: ' + str(args.column) +  ' "Motif" of column: ' + args.motif_to_use + ' Optimal allele: ' + str(args.opt_allele_to_use) + '\n')
    solution_file.write('Observed heterozygosity distribution: ' + str(obs_het_distr)+ '\n')
       
    accepted_params = []
    info = []
    
    list_all_params = []
    list_all_dist = []
    list_dist_vec = []
    sim_het_lists = []
    
    eps_het = 0.1
    eps_common = 0.1
    use_common = False
    
    for i in range(0, args.num_sims):
        
        a = np.random.uniform() 
        # Mu - Mean value of the underlying normal distribution
        # Sigma - Standard deviation of the underlying normal distribution
        mu, sigma = np.log(0.0003), np.log(30)
        mean = np.random.lognormal(mu, sigma)
        b = mean/a
        
        # Testing drawing mutation setting from prior
        if args.sim_model == 'all' and args.top_x == True:
            # Always s = 0
            a = 2 
            prior_mut = mut_settings[args.period]
            mut_setting = random.choice(prior_mut)
            ABC_tables_mut = ABC_tables_master_dic[mut_setting]
           
            to_accept, mean_of_differences, diff_vector = EstimateParamBinAgnostic(ABC_tables_mut, opt_allele_sub_list, a, b, obs_het_distr, obs_common_distr, eps_het, eps_common, use_common)

            list_all_dist.append(mean_of_differences)
                
            list_all_params.append(mut_setting)
             
        else:
   
            to_accept, mean_of_differences, diff_vector = EstimateParamBinAgnostic(ABC_tables, opt_allele_sub_list, a, b, obs_het_distr, obs_common_distr, eps_het, eps_common, use_common)

            list_all_params.append((a,b))
            list_all_dist.append((mean_of_differences))
                
            info.append((to_accept, a, b))
        
    # Implement top x% 
    if args.top_x == True:
        top_x_num = int(args.perc_acc * args.num_sims / 100)
        top_idx = np.argsort(list_all_dist)[:top_x_num]
        accepted_params = [list_all_params[i] for i in top_idx]
        
    solution_file.write('Number a,b pairs accepted: ' + str(len(accepted_params)) + '\n')
    
    if args.sim_model == 'all':
        solution_file.write('Accepted params: ' + ', '.join(str(item) for item in accepted_params) + '\n')
        
    else:
        if len(accepted_params) > 0:
            sort_mean = sorted(accepted_params, key=lambda x: x[0]*x[1])
            list_of_means = []
            for pair in sort_mean:
                list_of_means.append(pair[0]*pair[1])

            num_accepted = len(accepted_params)
            middle_index = int(num_accepted/2)

            solution_file.write('Sorted by mean - Median a, b:\n')
            solution_file.write(str(sort_mean[middle_index][0]) + ',' + str(sort_mean[middle_index][1])+ '\n')
            solution_file.write('Median mean: ' + str(np.median(list_of_means)) + '\n')
            mean_lower_bound = np.percentile(list_of_means, 2.5)
            mean_upper_bound = np.percentile(list_of_means, 97.5)
            solution_file.write('2.5 percentile mean: ' + str(mean_lower_bound) + '\n')
            solution_file.write('97.5 percentile mean: ' + str(mean_upper_bound) + '\n')

            solution_file.write('Accepted params: ' + ', '.join(str(item) for item in accepted_params) + '\n')
        
        solution_file.close()
    
if __name__ == '__main__':
    main()