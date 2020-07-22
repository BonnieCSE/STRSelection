import random
from ABC_functions import *

########## Joint Method Helper Functions ##########

def GetL1Norm(allele_freqs_prev, allele_freqs):
    
    diff_vector = np.absolute(allele_freqs_prev - allele_freqs)
    return np.sum(diff_vector)

def count_range_in_list(list_num, lower, upper):
    count = 0
    for elem in list_num:
        if elem >= lower and elem < upper: 
            count = count + 1
    return count

def getVector(het_distr):
    vec = [] 
    vec.append(count_range_in_list(het_distr, 0, 0.001)/len(het_distr))
    vec.append(count_range_in_list(het_distr, 0.001, 0.01)/len(het_distr))
    vec.append(count_range_in_list(het_distr, 0.01, 0.1)/len(het_distr))
    vec.append(count_range_in_list(het_distr, 0.1, 0.3)/len(het_distr))
    vec.append(count_range_in_list(het_distr, 0.3, 0.6)/len(het_distr))
    vec.append(count_range_in_list(het_distr, 0.6, 1)/len(het_distr))
    return vec

def get_bin_size(s):
    if s <= 0.000001:
        return 0.000001
    return s/10

def GetHet(table1, table2, s, bin_size): 
    het_list = []
    het_list_alt = []
    table_list = [table1, table2]
    for table in table_list:
        for combo in table:
            s_val = combo[0]
            het = combo[1]
            common = combo[2]
            if abs(s-s_val) < bin_size:
                if use_common == False:
                    het_list.append(het)
                if use_common == True and common == 0:    
                    het_list.append(het)
                if use_common == True and abs(common) <= 1:
                    het_list_alt.append(het)
                
    num_accepted = len(het_list)
    if num_accepted > 0:
        index = random.randrange(0, num_accepted, 1)
        return het_list[index]
    else:
        num_accepted = len(het_list_alt)
        index = random.randrange(0, num_accepted, 1)
        return het_list_alt[index]
        
# Get distribution of heterozygosity simulated
# Look up optimal allele, period, s
def EstimateParam(opt_allele_list, shape, scale, obs_mean, obs_var, obs_vec, obs_med, model1, model2, eps_mean, eps_var, eps_med, use_common):
    het_list = []
    for combo in opt_allele_list:
        
        s = np.random.gamma(shape, scale)
        if s > 1:
            s = 1

        # Perform simulation and get het; add to sim het distribution
        # Get from lookup table
        per = combo[0]
        optimal_ru = combo[1]
        if per == 3 and optimal_ru > 12:
            optimal_ru = 12
        if per == 3 and optimal_ru < 5:
            optimal_ru = 5
        if per == 4 and optimal_ru > 10:
            optimal_ru = 10
        if per == 4 and optimal_ru < 7:
            optimal_ru = 7
        if per == 2 and optimal_ru < 11:
            optimal_ru = 11
        if per == 2 and optimal_ru > 20:
            optimal_ru = 20
        bin_size = get_bin_size(s)
        
        file1 = '/gymreklab-tscc/bonnieh/abc/results/'+model1+'/' + str(per) + '_' + str(optimal_ru) + '.txt' 
        file2 = '/gymreklab-tscc/bonnieh/abc/results/'+model2+'/' + str(per) + '_' + str(optimal_ru) + '.txt'
        
        use_bin_type = 'nar'
       
        table1 = GetABCList(file1, use_bin_type)
        table2 = GetABCList(file2, use_bin_type)
        
        het = GetHet(table1, table2,  s, bin_size, use_common) # table3
        
        het_list.append(het)
        
    sim_mean = np.mean(het_list)
    sim_var = np.var(het_list)
    sim_vec = getVector(het_list)
    sim_med = np.median(het_list)
    array_obs = np.array(obs_vec)
    array_sim = np.array(sim_vec)
    
    if abs(obs_mean - sim_mean) < obs_mean/eps_mean and abs(obs_var - sim_var) < obs_var/eps_var and abs(obs_med - sim_med) < obs_med/eps_med:
        return True, het_list
    else:
        return False, het_list