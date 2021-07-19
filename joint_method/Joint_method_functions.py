# This file contains helper functions for performing SISTR Joint
# to infer the distribution of s across a set of STRs.

########## Imports ##########

import random
import time
import statistics
import numpy as np

########## Joint Method Helper Functions ##########
    
### Process allele frequencies ###
# Return optimal allele repeat units, allele frequencies
def Process_Freqs(freqs, per, end, start, return_freqs=True, motif=False):
    # Process freqs
    freqs_list = freqs.split(',')

    # Create dictionary of freqs: key is allele and value is freq
    freqs_dic = {}
    
    for pair in freqs_list:
        # Processing allele frequencies in non-motif (alternate) format
        if motif == False:
            allele = int(pair.split(':')[0])
            
        # Processing allele frequencies in motif format
        
        else:
            allele = pair.split(':')[0]
            allele = len(allele)
        freq = int(pair.split(':')[1])
        freqs_dic[int(allele/per)] = freq
        
    # Note: Actual population size is actually half the value of the variable pop_size 
    # because pop_size is calculated by adding up diploid allele freqs
    pop_size = 0
    for elem in freqs_dic:
        pop_size = pop_size + freqs_dic[elem]
   
    # Get optimal allele (how many ru away from reference allele) - allele with highest frequency
    opt_allele_rel = max(freqs_dic, key=freqs_dic.get)

    # Get info about reference allele (length in base pairs, repeat units)
    # Get optimal allele repeat unit length
    if motif == False:
        ref_length_bp = end - start + 1
        ref_length_ru = int(ref_length_bp/per)
        opt_allele = ref_length_ru + opt_allele_rel
    
    else:
        opt_allele = opt_allele_rel
    
    if return_freqs == False:
        return opt_allele
    
    else:
        
        # Get allele frequencies
        freqs_dic_final = {}
        allele_list = []
        
        for allele in freqs_dic:
            new_allele = allele - opt_allele_rel
            freqs_dic_final[new_allele] = freqs_dic[allele]
            allele_list.append(new_allele)

        # Get highest absolute value in list
        max_allele = abs(max(allele_list, key=abs))

        # Put allele freqs in list
        allele_freqs = [0] * (2*max_allele+1) 
        for key in freqs_dic_final:
            allele_freq = freqs_dic_final[key]/pop_size
            allele_freqs[max_allele+key] = allele_freq

        return opt_allele, allele_freqs

### Get list of s with corresponding summary statistics ###
def GetABCList(abcFile, num_bins):
    abc_file = open(abcFile, 'r')
    header = abc_file.readline().strip().split('\t')
    
    abc_list = []
    
    if num_bins == 0:
        for line in abc_file:
            info = line.strip().split('\t')
            s = float(info[0])
            het = float(info[1])
            common = int(info[2])
            stats_list = [s, het, common]
            abc_list.append(stats_list)
     
    # Get bins summary statistic
    else:
        # Get column number of freqs column in file
        freqs_column = 0
        for i in range(0, len(header)):
            if header[i] == 'freqs':
                freqs_column = i

        for line in abc_file:
            info = line.strip().split('\t')
            s = float(info[0])
            freq_string = info[freqs_column]
            allele_freqs = [float(freq) for freq in freq_string.split(',')]

            abc_het = 1-sum([item**2 for item in allele_freqs])
            abc_common = len([i for i in allele_freqs if i >= 0.05]) 
            abc_bins = GetBins(allele_freqs, num_bins)

            stats_list = [s, abc_het, abc_common, abc_bins]
            abc_list.append(stats_list)
        
    abc_file.close()
    return abc_list
    
# Out of all s values in list, get s value closest to given s
def getNearestS(s_ABC_round, s_list_available):
    min_dist = 100000000
    nearest_s = -2
    for elem in s_list_available:
        dist = abs(s_ABC_round - elem)
        if dist < min_dist:
            min_dist = dist
            nearest_s = elem
    return nearest_s
    
# Get LRT bin for given s
def get_LRT_bin(s):
    if s < 0.00001:
        return 0
    if s <= 0.0001:
        return round(s, 5)
    if s <= 0.001:
        return round(s, 4)
    if s <= 0.01:
        return round(s, 3)
    if s <= 0.1:
        return round(s,2)
    return round(s,1)
    
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
    return reverse, complement, reverse_complement

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