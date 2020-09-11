import random
from LRT_functions import *
import time
import statistics
########## Joint Method Helper Functions ##########

def GetGammaBins(k, theta, num_sims, return_median=False):
    #num_sims = 100
    s = [] # List of s values drawn from gamma distribution
    for i in range(0, num_sims):
        s_val = np.random.gamma(k, theta)
        if s_val > 1:
            s_val = 1
        s.append(s_val)

    # List of binned s values
    # Bins: 0<=s<10^-4, 10^-4<=s<10^-3, 10^-3<=s<10^-2, s>=10^-2
    svals = [0, 0, 0, 0] 
    for elem in s:
        if elem < 10**-4:
            svals[0] += 1
        elif elem < 10**-3:
            svals[1] += 1
        elif elem < 10**-2:
            svals[2] += 1
        elif elem >= 10**-2:
            svals[3] += 1
    
    total = sum(svals)
    for i in range(0, len(svals)):
        svals[i] = svals[i]/total
    if return_median == False:
        return svals
    else:
        return s, np.median(s), np.mean(s)

    

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
    #if s <= 10**-6: 
        #return 10**-6
    if s <= 0.0001:
        return 0.00001
    if s <= 0.001:
        return 0.0001
    if s <= 0.01:
        return 0.001
    if s <= 0.2:
        return 0.01
    return 0.1
    #return s/10

def GetHetandCommon(table, s, bin_size): 
    het_list = []
    common_list = []
    for combo in table:
        s_val = combo[0]
        if abs(s-s_val) < bin_size:
            het = combo[1]
            common = combo[2]
            het_list.append(het)
            common_list.append(common)
               
    num_accepted = len(het_list)
    index = random.randrange(0, num_accepted, 1)
    return het_list[index], common_list[index]

def GetHetCommon(table):
    total = len(table)
    index = random.randrange(0, total, 1)
    return table[index][0], table[index][1]
   
# Get distribution of simulated heterozygosity and common allele
# Look up optimal allele, period, s
def EstimateParam(ABC_tables, opt_allele_list, shape, scale, obs_het_stats, \
                  obs_common_stats, model, eps_het, eps_common, use_common_alleles, return_lists=False, return_info=False):
    het_list = []
    common_list = []
    time1 = 0
    time2 = 0
    for combo in opt_allele_list:
        
        s = np.random.gamma(shape, scale)
        if s > 1:
            s = 1

        # Perform simulation and get het and common allele; add summary statistics to distribution
        # Get from lookup table
        per = combo[0]
        optimal_ru = combo[1]
        if per == 3 and optimal_ru > 13:
            optimal_ru = 13
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
            
        #bin_size = get_bin_size(s)
        
        #file = '/gymreklab-tscc/bonnieh/abc/results/'+model+'/' + str(per) + '_' + str(optimal_ru) + '.txt' 
        
        num_bins = 0
        
        first = time.time()
        #table = ABC_tables[optimal_ru] #GetABCList(file, num_bins)
        dic_summ_stats = ABC_tables[optimal_ru] #GetABCList(file, num_bins)
        second = time.time()
        t1 = second-first
        time1 = time1 + t1
        third = time.time()
        #het, common = GetHetandCommon(table, s, bin_size) 
        s_round = get_LRT_bin(s)
        ## NEW
        if s_round not in dic_summ_stats:
            s_list_available = []
            for elem in dic_summ_stats:
                s_list_available.append(elem)
            s_round = getNearestS(s_round, s_list_available)
        ##
        table = dic_summ_stats[s_round]
        pair = random.choice(table)
        het, common = pair[0], pair[1]#GetHetCommon(table)
        fourth = time.time()
        t2 = fourth-third
        time2 = time2 + t2
        het_list.append(het)
        common_list.append(common)
    
    if return_lists == True:
        return het_list, common_list
        
    else:
        sim_mean_het = np.mean(het_list)
        sim_var_het = np.var(het_list)
        sim_med_het = np.median(het_list)
        
        sim_mean_common = np.mean(common_list)
        sim_var_common = np.var(common_list)
        sim_med_common = np.median(common_list)
        
        #print((shape, scale, sim_mean_het, sim_var_het, sim_med_het))
        '''
        params_accept = [1,1,1,1,1,1]
        if abs(obs_het_stats[0] - sim_mean_het) < obs_het_stats[0]/eps_het[0]:
            params_accept[0] = 0
        
        if abs(obs_het_stats[1] - sim_var_het) < obs_het_stats[1]/eps_het[1]:
            params_accept[1] = 0
            
        if abs(obs_het_stats[2] - sim_med_het) < obs_het_stats[2]/eps_het[2]:
            params_accept[2] = 0
            
        if abs(obs_common_stats[0] - sim_mean_common) < obs_common_stats[0]/eps_common[0]:
            params_accept[3] = 0
            
        if abs(obs_common_stats[1] - sim_var_common) < obs_common_stats[1]/eps_common[1]:
            params_accept[4] = 0
            
        if abs(obs_common_stats[2] - sim_med_common) < obs_het_stats[2]/eps_common[2]:
            params_accept[5] = 0
            
        if use_common_alleles == False:
            if 1 in params_accept[0:3]:
                return False, [sim_mean_het, sim_var_het, sim_med_het], params_accept
            else:
                return True, [sim_mean_het, sim_var_het, sim_med_het], params_accept
                
            
        '''
        
        '''
        if abs(obs_het_stats[0] - sim_mean_het) < obs_het_stats[0]/eps_het[0] and abs(obs_het_stats[1] - sim_var_het) < obs_het_stats[1]/eps_het[1] and abs(obs_het_stats[2] - sim_med_het) < obs_het_stats[2]/eps_het[2]:
            if use_common_alleles == False:
                return True, time1, time2
            if abs(obs_common_stats[0] - sim_mean_common) < obs_common_stats[0]/eps_common[0] and abs(obs_common_stats[1] - sim_var_common) < obs_common_stats[1]/eps_common[1] and abs(obs_common_stats[2] - sim_med_common) < obs_het_stats[2]/eps_common[2]:
                return True, time1, time2
            else:
                return False, time1, time2
        else:
            return False, time1, time2
        '''
        if return_info == False:
            if abs(obs_het_stats[0] - sim_mean_het) < (obs_het_stats[0] + 0.005)/eps_het[0] and abs(obs_het_stats[1] - sim_var_het) < (obs_het_stats[1])/eps_het[1] and abs(obs_het_stats[2] - sim_med_het) < (obs_het_stats[2]+0.005)/eps_het[2]:
                if use_common_alleles == False:
                    return True, time1, time2
                if abs(obs_common_stats[0] - sim_mean_common) < (obs_common_stats[0])/eps_common[0] and abs(obs_common_stats[1] - sim_var_common) < (obs_common_stats[1])/eps_common[1] and abs(obs_common_stats[2] - sim_med_common) < (obs_het_stats[2])/eps_common[2]:
                    return True, time1, time2
                else:
                    return False, time1, time2
            else:
                return False, time1, time2
        
        else:
            params_accept = [1,1,1,1,1,1]
            count_0 = 0
            if abs(obs_het_stats[0] - sim_mean_het) < (obs_het_stats[0] + 0.05)/eps_het[0]:
                params_accept[0] = 0
                count_0 = count_0 + 1
            if abs(obs_het_stats[1] - sim_var_het) < (obs_het_stats[1] + 0.05)/eps_het[1]:
                params_accept[1] = 0
                count_0 = count_0 + 1
            if abs(obs_het_stats[2] - sim_med_het) < (obs_het_stats[2] + 0.05)/eps_het[2]:
                params_accept[2] = 0
                count_0 = count_0 + 1
            if count_0 == 3:
                return True, sim_mean_het, sim_var_het, sim_med_het, params_accept[0], params_accept[1], params_accept[2]
            else:
                return False, sim_mean_het, sim_var_het, sim_med_het, params_accept[0], params_accept[1], params_accept[2]
        '''
        params_accept = [1,1,1,1,1,1]
        if abs(obs_het_stats[0] - sim_mean_het) < (obs_het_stats[0]+0.005)/eps_het[0]:
            params_accept[0] = 0
        
        if abs(obs_het_stats[1] - sim_var_het) < obs_het_stats[1]/eps_het[1]:
            params_accept[1] = 0
            
        if abs(obs_het_stats[2] - sim_med_het) < (obs_het_stats[2]+0.005)/eps_het[2]:
            params_accept[2] = 0
            
        if abs(obs_common_stats[0] - sim_mean_common) < obs_common_stats[0]/eps_common[0]:
            params_accept[3] = 0
            
        if abs(obs_common_stats[1] - sim_var_common) < obs_common_stats[1]/eps_common[1]:
            params_accept[4] = 0
            
        if abs(obs_common_stats[2] - sim_med_common) < obs_het_stats[2]/eps_common[2]:
            params_accept[5] = 0
            
        if use_common_alleles == False:
            if 1 in params_accept[0:3]:
                return False, [sim_mean_het, sim_var_het, sim_med_het], params_accept
            else:
                return True, [sim_mean_het, sim_var_het, sim_med_het], params_accept
        '''
                