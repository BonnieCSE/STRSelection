# Script to run ABC for all-locus/joint method 

### Imports ###

import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from Joint_method_functions import *
PLOTDIR = '/storage/BonnieH/selection_project/joint_method/results/'
import time
import statistics

def main():
    inFile = '/storage/BonnieH/selection_project/ssc_files/0810/allele_freqs_filt_by_total.txt'
    allele_freqs_file = open(inFile, 'r')
    
    period = int(sys.argv[1]) # Which period to jointly estimate s on
    column = int(sys.argv[2]) # Which column of file to jointly estimate s on
    motif_to_use = sys.argv[3] # Which motif to jointly estimate s on
    
    model = sys.argv[4] # File from which to obtain s values and summary statistics 
    
    eps_mean_het = float(sys.argv[5]) # Denominator for epsilon for mean of het distr
    eps_var_het = float(sys.argv[6]) # Denominator for epsilon for variance of het distr
    eps_med_het = float(sys.argv[7]) # Denominator for epsilon for median of het distr
    
    eps_mean_common = float(sys.argv[8]) # Denominator for epsilon for mean of common distr
    eps_var_common = float(sys.argv[9]) # Denominator for epsilon for variance of common distr
    eps_med_common = float(sys.argv[10]) # Denominator for epsilon for median of common distr
    
    num_sims = int(sys.argv[11]) # Number of simulations
    outFolder = sys.argv[12] # Name of outfolder
    use_common_alleles = int(sys.argv[13]) # Whether to use common alleles
    if use_common_alleles == 0:
        use_common_alleles = True
    else:
        use_common_alleles = False
        
    num_loci = int(sys.argv[14])
    opt_allele_to_validate = int(sys.argv[15])
    
    num_bins = 0
    
    # Naming file
    filename = PLOTDIR + 'results/'
    
    filename = filename + outFolder + '/'
    solution_file = open(filename + 'per_%d_%d_%s_sims_%d_het_eps_%d_%d_%d_comm_eps_%d_%d_%d_num_loci_%d_opt_%d.txt'%(period, column, motif_to_use, num_sims, eps_mean_het, eps_var_het, eps_med_het, eps_mean_common, eps_var_common, eps_med_common, num_loci, opt_allele_to_validate), 'w') 
    
    solution_file.write("Num sims: " + str(num_sims) + '\n')
    
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
    
    motif_dic = {} # Dictionary of motifs; Key motif; Value number of motifs
    
    # Get all motifs
    for line in allele_freqs_file:
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        freqs = info[3]
        per = int(info[5])
        motif = info[column]
        
        if per == period:
            opt_allele = Process_Freqs(freqs, per, end, start, False)
            if opt_allele_to_validate == 0 or opt_allele_to_validate == opt_allele:
                if motif not in motif_dic:
                    motif_dic[motif] = 1
                else:
                    motif_dic[motif] = motif_dic[motif] + 1
              
    allele_freqs_file.close()
    print(motif_dic)
    
    allele_freqs_file = open(inFile, 'r')
    
    obs_het_distr_all = []
    obs_common_distr_all = []
    opt_allele_list = []
    count = 0
   
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
        
        motif = info[column]
        opt_allele, allele_freqs = Process_Freqs(freqs, per, end, start)
        
        if motif == motif_to_use and per == period and (opt_allele_to_validate == 0 or opt_allele_to_validate == opt_allele):
            count = count + 1
            
            
            obs_het = 1-sum([item**2 for item in allele_freqs])
            obs_common = len([i for i in allele_freqs if i >= 0.05])
            obs_het_distr_all.append(obs_het)
            obs_common_distr_all.append(obs_common)
            opt_allele_list.append((per, opt_allele))
                   
    allele_freqs_file.close()
   
    opt_allele_sub_list = []
    obs_het_distr = []
    obs_common_distr = []
    if len(opt_allele_list) <= num_loci:
        opt_allele_sub_list = opt_allele_list
        obs_het_distr = obs_het_distr_all
        obs_common_distr = obs_common_distr_all
        
    else:
        
        divisor = int(len(opt_allele_list)/num_loci) + 2
        num_loci_added = 0
        while num_loci_added < num_loci:
            divisor = divisor * num_loci_added + 5
            if divisor > len(opt_allele_list):
                divisor = divisor%len(opt_allele_list)
            
            toAdd = opt_allele_list.pop(divisor)
            opt_allele_sub_list.append(toAdd)
            toAddHet = obs_het_distr_all.pop(divisor)
            obs_het_distr.append(toAddHet)
            toAddComm = obs_common_distr_all.pop(divisor)
            obs_common_distr.append(toAddComm)
            num_loci_added = num_loci_added + 1
        print('Finished getting random set of opt alleles')
        
    obs_mean_het = np.mean(obs_het_distr)
    obs_var_het = np.var(obs_het_distr)
    obs_med_het = np.median(obs_het_distr)
    obs_het_stats = [obs_mean_het, obs_var_het, obs_med_het]
    
    obs_mean_common = np.mean(obs_common_distr)
    obs_var_common = np.var(obs_common_distr)
    obs_med_common = np.median(obs_common_distr)
    obs_common_stats = [obs_mean_common, obs_var_common, obs_med_common]
    eps_het = [eps_mean_het, eps_var_het, eps_med_het]
    eps_common = [eps_mean_common, eps_var_common, eps_med_common]
        
    solution_file.write('Number of loci used: ' + str(len(obs_het_distr)) + '\n')
    solution_file.write('Column: ' + str(column) +  ' "Motif" of column: ' + motif_to_use + ' Optimal allele: ' + str(opt_allele_to_validate) + '\n')
    solution_file.write('Mean of observed heterozygosity: ' + str(obs_mean_het) + '\n')
    solution_file.write('Variance of observed heterozygosity: ' + str(obs_var_het) + '\n')
    solution_file.write('Median of observed heterozygosity: ' + str(obs_med_het) + '\n')
    solution_file.write('Mean of observed number of common alleles: ' + str(obs_mean_common) + '\n')
    solution_file.write('Variance of observed number of common alleles: ' + str(obs_var_common) + '\n')
    solution_file.write('Median of observed number of common alleles: ' + str(obs_med_common) + '\n')
       
    accepted_params = []
    info = []
    total_time_1 = 0
    total_time_2 = 0
    all_time = 0
    for i in range(0, num_sims):
            
        k = np.random.uniform() 
        
        mu, sigma = np.log(0.0003), np.log(30)
        mean = np.random.lognormal(mu, sigma)
        theta = mean/k
        t1 = time.time()
        toAdd, mean_fit, var_fit, med_fit, sim_mean_het, sim_var_het, sim_med_het = EstimateParam(ABC_tables, \
                                       opt_allele_sub_list, k, theta, obs_het_stats, obs_common_stats, model, \
                                       eps_het, eps_common, use_common_alleles, False, True)
        t2 = time.time()
        all_time = all_time + t2-t1
        
        info.append((mean_fit, var_fit, med_fit, sim_mean_het, sim_var_het, sim_med_het))
        if toAdd == True:
            accepted_params.append((k, theta))
            
    solution_file.write('Number k,theta pairs accepted: ' + str(len(accepted_params)) + '\n')
    
    if len(accepted_params) > 0:
        sort_mean = sorted(accepted_params, key=lambda x: x[0]*x[1])
        list_of_means = []
        for pair in sort_mean:
            list_of_means.append(pair[0]*pair[1])

        num_accepted = len(accepted_params)
        middle_index = int(num_accepted/2)

        solution_file.write('Sorted by mean - Median k, theta:\n')
        solution_file.write(str(sort_mean[middle_index][0]) + ',' + str(sort_mean[middle_index][1])+ '\n')
        solution_file.write('Median mean: ' + str(np.median(list_of_means)) + '\n')
        mean_lower_bound = np.percentile(list_of_means, 2.5)
        mean_upper_bound = np.percentile(list_of_means, 97.5)
        solution_file.write('2.5 percentile mean: ' + str(mean_lower_bound) + '\n')
        solution_file.write('97.5 percentile mean: ' + str(mean_upper_bound) + '\n')

        solution_file.write('Accepted params: ' + ', '.join(str(item) for item in accepted_params) + '\n')
    solution_file.write('Eps mean: '  + str((obs_het_stats[0] + 0.05)/eps_het[0]) + 'Eps var: ' + str((obs_het_stats[1] + 0.05)/eps_het[1]) + 'Eps med: ' + str((obs_het_stats[2] + 0.05)/eps_het[2]) + '\n') # 0.005
    
    solution_file.write('Optimal allele list: ' + ','.join(str(item) for item in opt_allele_sub_list) + '\n')
    solution_file.write('Info: ' + '\n'.join(str(item) for item in info) + '\n')
    solution_file.close()
    
if __name__ == '__main__':
    main()