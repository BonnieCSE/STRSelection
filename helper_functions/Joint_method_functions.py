# This file contains helper functions for performing the joint setting for SISTR
# to infer the distribution of s across a set of STRs.

########## Imports ##########
from LRT_functions import *
import random
import time
import statistics

########## Joint Method Helper Functions ##########
    
# Get difference vector
def GetDiffVector(vector1, vector2, normalized = False, const = 0):
    diff_vector = []
    for i in range(0, len(vector1)):
        distance = abs(vector1[i] - vector2[i])
        if normalized == False:
            diff_vector.append(distance)
        else:
            if distance == 0:
                diff_vector.append(distance)
            else:
                mean = (vector1[i] + vector2[i])/2 + const

                diff_vector.append(distance/mean)
    return diff_vector

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
                
# Function to get distribution of simulated heterozygosity and number of common alleles
def GetLists(ABC_tables, opt_allele_list, a, b):
    
    het_list = []
    common_list = []
    
    for combo in opt_allele_list:
        
        if a == 2:
            s = 0
        else:
            s = np.random.gamma(a, b)
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
    
    return het_list, common_list
        
# Return whether to accept shape, scale pair
def EstimateParamBinAgnostic(ABC_tables, opt_allele_list, a, b, obs_het_list, obs_common_list, eps_het, eps_common, use_common, normalized=False, const=0):
    
    het_list, common_list = GetLists(ABC_tables, opt_allele_list, a, b)
    
    het_list.sort()
    obs_het_list.sort()
    
    if normalized == False:
        diff_vector = GetDiffVector(het_list, obs_het_list)
        
    else:
        diff_vector = GetDiffVector(het_list, obs_het_list, True, const)
    
    mean_of_differences = np.mean(diff_vector)
    
    to_accept = False
    
    if use_common == False:
        if mean_of_differences < eps_het:
            to_accept = True

        return to_accept, mean_of_differences, diff_vector
    
    else:
        
        if mean_of_differences < eps_het and dist_common < eps_common:
            to_accept = True
            
        return to_accept, mean_of_differences, diff_vector