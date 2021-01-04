# Generate new abc lookup tables from eurodem_prior2 with simulated genotyping errors
# New name: eurodem_prior2_errors

# Imports 
import sys
sys.path.append("/storage/BonnieH/selection_project/helper_functions")
from ABC_functions import *
import random

def SimErrors(allele_freqs):
    
    end_samp_n = 6500
    
    allele_freqs_nums = []
    
    for elem in allele_freqs:
        count = int(elem*end_samp_n)
        allele_freqs_nums.append(count)
        
    allele_freqs_errors = [0] * len(allele_freqs)
    
    for i in range(0, len(allele_freqs_nums)):
        elem = allele_freqs_nums[i]
        for j in range(0, elem):
            random_num = random.randint(1,101)
            if random_num == 2 and i != 0 and i != 24:
                rand_num_2 = random.randint(1,3)
                if rand_num_2 == 1:
                    new_index = i + 1
                else:
                    new_index = i - 1
                allele_freqs_errors[new_index] = allele_freqs_errors[new_index] + 1
            else:
                allele_freqs_errors[i] = allele_freqs_errors[i] + 1
                
    allele_freqs_errors_final = []
    for elem in allele_freqs_errors:
        frac = elem/end_samp_n
        allele_freqs_errors_final.append(frac)
            
    return allele_freqs_errors_final

def WriteErrorsFile():
    pers = [2,3,4]
    opt_allele_dic = {}
    opt_allele_dic[2] = np.arange(11,21,1)
    opt_allele_dic[3] = np.arange(5,14,1)
    opt_allele_dic[4] = np.arange(7,11,1)
    
    folder = 'eurodem_prior2'
    
    for per in pers:
        for opt_allele in opt_allele_dic[per]:
            
            outFile = '/gymreklab-tscc/bonnieh/abc/results/' + folder + '_errors' + '/' + str(per) + '_' + str(opt_allele) + '.txt'
            results = open(outFile, "w")

            # Write results header
            results.write("s" + "\t" + "het" + "\t" + "common" +"\t" + "bins" + "\t" + "freqs" + "\n")
            
            abcFile  = '/gymreklab-tscc/bonnieh/abc/results/' + folder + '/' + str(per) + '_' + str(opt_allele) + '.txt'
            abc_file = open(abcFile, 'r')
            header = abc_file.readline().strip().split('\t')
            
            freqs_column = 0
            for i in range(0, len(header)):
                if header[i] == 'freqs':
                    freqs_column = i
            
            count = 0
            for line in abc_file:
                count = count + 1
                #if count%1000 == 0:
                    #print(count)
                info = line.strip().split('\t')
                
                s = float(info[0])
                
                freq_string = info[freqs_column]
                
                allele_freqs = [float(freq) for freq in freq_string.split(',')]
                allele_freqs_errors = SimErrors(allele_freqs)
                    
                het = 1-sum([item**2 for item in allele_freqs_errors])
                common = len([i for i in allele_freqs_errors if i >= 0.05]) 
                bins = GetBins(allele_freqs_errors, 5)
                results.write(str(s) + "\t" + str(het) + "\t" + str(common) + "\t" + ','.join(str(item) for item in bins) + "\t" + ','.join(str(item) for item in allele_freqs_errors) + "\n")
                

            abc_file.close()   
            results.close()
            
def main():
    WriteErrorsFile()

if __name__ == '__main__':
    main()    