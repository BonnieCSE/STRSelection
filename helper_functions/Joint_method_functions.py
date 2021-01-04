# This file contains helper functions for performing the joint setting for SISTR
# to infer the distribution of s across a set of STRs.

########## Imports ##########
import random
from LRT_functions import *
import time
import statistics

########## Joint Method Helper Functions ##########

# Rotate a string by n characterss
def rotate(text, n):
    return text[n:] + text[:n]

# Get reverse, complement, and reverse complement of string
def ReverseComplement(text):
    nuc_map = {}
    nuc_map['A'] = 'T'
    nuc_map['T'] = 'A'
    nuc_map['G'] = 'C'
    nuc_map['C'] = 'G'
    
    reverse = ""
    for nuc in text:
        reverse = nuc + reverse
        
    complement = ""
    for nuc in text:
        comp_map = nuc_map[nuc]
        complement = complement + comp_map

    reverse_complement = "" 
    for nuc in complement: 
         reverse_complement = nuc + reverse_complement
    return reverse,complement, reverse_complement

# Return canonical repeat unit
# The canonical repeat unit is defined as the lexicographically first repeat unit when considering all rotations and strand orientations of the repeat sequence. 
# For example, the canonical repeat unit for the repeat sequence CAGCAGCAGCAG would be AGC.
def GetCanonicalRU(motif, return_all=False):
    all_versions = [motif]
    
    # Get all rotations and strand orientations of the motif
    length = len(motif)
    
    # Rotate length times
    for i in range(0, 1):
        rotated_motif = rotate(motif, i)
        all_versions.append(rotated_motif)
        reverse, complement, reverse_complement = ReverseComplement(rotated_motif)
        all_versions.append(reverse)
        all_versions.append(complement)
        all_versions.append(reverse_complement)
    
    all_versions.sort()
    if return_all == False:
        return(all_versions[0])
    else:
        return all_versions
    
# Function to get bins of s values drawn from gamma distribution
def GetGammaBins(a, b, num_sims, return_median=False):
    
    s = [] # List of s values drawn from gamma distribution
    
    # Draw s values from gamma distribution with parameters a, b
    for i in range(0, num_sims):
        s_val = np.random.gamma(a, b)
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

# Get L1 Norm of two vectors
def GetL1Norm(allele_freqs_prev, allele_freqs):
    diff_vector = np.absolute(allele_freqs_prev - allele_freqs)
    return np.sum(diff_vector)

# Get number of elements in list that are greater than or equal than lower and less than upper
def count_range_in_list(list_num, lower, upper):
    count = 0
    for elem in list_num:
        if elem >= lower and elem < upper: 
            count = count + 1
    return count

# Get vector of heterozygosities given distribution
def getHetVector(distr, fine_grained=False):
    vec = [] 
    
    if fine_grained == False:
        vec.append(count_range_in_list(distr, 0, 0.1)/len(distr))
        vec.append(count_range_in_list(distr, 0.1, 0.4)/len(distr))
        vec.append(count_range_in_list(distr, 0.4, 1)/len(distr))
    else:
        
        vec.append(count_range_in_list(distr, 0, 0.002)/len(distr))
        vec.append(count_range_in_list(distr, 0.002, 0.1)/len(distr))
        vec.append(count_range_in_list(distr, 0.1, 0.4)/len(distr))
        vec.append(count_range_in_list(distr, 0.4, 0.8)/len(distr))
        vec.append(count_range_in_list(distr, 0.8, 1)/len(distr))
        '''
        
        vec.append(count_range_in_list(distr, 0, 0.001)/len(distr))
        vec.append(count_range_in_list(distr, 0.001, 0.02)/len(distr))
        vec.append(count_range_in_list(distr, 0.02, 0.1)/len(distr))
        vec.append(count_range_in_list(distr, 0.1, 0.3)/len(distr))
        vec.append(count_range_in_list(distr, 0.3, 0.6)/len(distr))
        vec.append(count_range_in_list(distr, 0.6, 0.7)/len(distr))
        vec.append(count_range_in_list(distr, 0.7, 0.8)/len(distr))
        vec.append(count_range_in_list(distr, 0.8, 1)/len(distr))
        '''
        '''
        vec.append(count_range_in_list(distr, 0, 0.1)/len(distr))
        vec.append(count_range_in_list(distr, 0.1, 0.3)/len(distr))
        vec.append(count_range_in_list(distr, 0.3, 0.6)/len(distr))
        vec.append(count_range_in_list(distr, 0.6, 1)/len(distr))
        '''
        # Old fine grained vecotr
        '''
        vec.append(count_range_in_list(distr, 0, 0.1)/len(distr))
        vec.append(count_range_in_list(distr, 0.1, 0.3)/len(distr))
        vec.append(count_range_in_list(distr, 0.3, 0.5)/len(distr))
        vec.append(count_range_in_list(distr, 0.5, 0.7)/len(distr))
        vec.append(count_range_in_list(distr, 0.7, 1)/len(distr))
        '''
        # Old vector 
        #vec.append(count_range_in_list(distr, 0, 0.001)/len(distr))
        #vec.append(count_range_in_list(distr, 0.001, 0.01)/len(distr))
        #vec.append(count_range_in_list(distr, 0.01, 0.1)/len(distr))
        
    return vec

# Get vector of number of common alleles given distribution
def getCommonVector(distr):
    vec = [] 
    vec.append(count_range_in_list(distr, 0, 2)/len(distr))
    vec.append(count_range_in_list(distr, 2, 4)/len(distr))
    vec.append(count_range_in_list(distr, 4, 100)/len(distr))
    return vec

# Get bin size given s value
def get_bin_size(s):
    if s <= 0.0001:
        return 0.00001
    if s <= 0.001:
        return 0.0001
    if s <= 0.01:
        return 0.001
    if s <= 0.2:
        return 0.01
    return 0.1
    
    #if s <= 10**-6: 
        #return 10**-6
    #return s/10

'''
# Get heterozygosity and number of common alleles from table given s value and bin size
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
'''

'''
def GetHetCommon(table):
    total = len(table)
    index = random.randrange(0, total, 1)
    return table[index][0], table[index][1]
'''
   
# return_lists = True -> Function to get distribution of simulated heterozygosity and number of common alleles
# Return whether to accept shape, scale pair
def EstimateParam(ABC_tables, opt_allele_list, shape, scale, obs_het_stats, obs_common_stats, obs_het_vec, obs_common_vec, model, \
                  eps_het, eps_common, use_common_alleles, use_bins, eps_bins_het, eps_bins_common, return_lists=False, \
                  return_info=False, het_fine_grained=False):
    
    het_list = []
    common_list = []
    
    for combo in opt_allele_list:
        
        s = np.random.gamma(shape, scale)
        if s > 1:
            s = 1

        # Get heterozygosity and number of common alleles from lookup table; add summary statistics to het_list and common_list
       
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
        
        num_bins = 0
        
        
        dic_summ_stats = ABC_tables[optimal_ru] 
        s_round = get_LRT_bin(s)
        
        if s_round not in dic_summ_stats:
            s_list_available = []
            for elem in dic_summ_stats:
                s_list_available.append(elem)
            s_round = getNearestS(s_round, s_list_available)
        
        table = dic_summ_stats[s_round]
        pair = random.choice(table)
        het, common = pair[0], pair[1]
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
        
        sim_het_vec = getHetVector(het_list, het_fine_grained)
        sim_common_vec = getCommonVector(common_list)
        
        if use_bins == False:
            params_accept = [1,1,1,1,1,1]
            count_0 = 0
            if abs(obs_het_stats[0] - sim_mean_het) < (obs_het_stats[0] + 0.005)/eps_het[0]:
                params_accept[0] = 0
                count_0 = count_0 + 1
            if abs(obs_het_stats[1] - sim_var_het) < (obs_het_stats[1] + 0.005)/eps_het[1]:
                params_accept[1] = 0
                count_0 = count_0 + 1
            if abs(obs_het_stats[2] - sim_med_het) < (obs_het_stats[2] + 0.005)/eps_het[2]:
                params_accept[2] = 0
                count_0 = count_0 + 1
            '''
            if use_common_alleles == False:
                if count_0 == 3:
                    return True, sim_mean_het, sim_var_het, sim_med_het, params_accept[0], params_accept[1], params_accept[2]
                else:
                    return False, sim_mean_het, sim_var_het, sim_med_het, params_accept[0], params_accept[1], params_accept[2]
            '''
            if use_common_alleles == True:
                if abs(obs_common_stats[0] - sim_mean_common) < (obs_common_stats[0] + 0.005)/eps_common[0]:
                    params_accept[3] = 0
                    count_0 = count_0 + 1
                if abs(obs_common_stats[1] - sim_var_common) < (obs_common_stats[1] + 0.005)/eps_common[1]:
                    params_accept[4] = 0
                    count_0 = count_0 + 1
                if abs(obs_common_stats[2] - sim_med_common) < (obs_common_stats[2] + 0.005)/eps_common[2]:
                    params_accept[5] = 0
                    count_0 = count_0 + 1
                if count_0 == 6:
                    return True, sim_mean_het, sim_med_het, sim_mean_common, sim_med_common, params_accept[0], params_accept[1], params_accept[2], params_accept[3], params_accept[4], params_accept[5]
                else:
                    return False, sim_mean_het, sim_med_het, sim_mean_common, sim_med_common, params_accept[0], params_accept[1], params_accept[2], params_accept[3], params_accept[4], params_accept[5]
        
        else:
            info_accept = [1,1,1,1]
            count_0 = 0
            dist_het = GetVectorDistance(sim_het_vec, obs_het_vec)
            if dist_het < eps_bins_het:
                count_0 = count_0 + 1
                info_accept[0] = 0
            dist_common = GetVectorDistance(sim_common_vec, obs_common_vec)
            if dist_common < eps_bins_common:
                count_0 = count_0 + 1
                info_accept[1] = 0
            info_accept[2] = dist_het
            info_accept[3] = dist_common
            if count_0 == 2:
                return True, sim_het_vec[0], sim_het_vec[1], sim_het_vec[2], sim_common_vec[0], sim_common_vec[1], sim_common_vec[2], info_accept[0], info_accept[1], info_accept[2], info_accept[3]
            else:
                return False, sim_het_vec[0], sim_het_vec[1], sim_het_vec[2], sim_common_vec[0], sim_common_vec[1], sim_common_vec[2], info_accept[0], info_accept[1], info_accept[2], info_accept[3]       