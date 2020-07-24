import sys 
sys.path.append("/projects/ps-gymreklab/bonnieh/helper_functions")
from Simulation_functions import *
from ABC_functions import *

def main():    
    # Load parameters
    per = int(sys.argv[1])
    opt_allele = int(sys.argv[2])
    s_vals = sys.argv[3]
    LRT_num_sims = int(sys.argv[4])
    s_list = [float(s) for s in s_vals.split(',')]
    outFolder1 = 'const_prelim/'
    outFolder2 = 'euro_prelim/'
    outFile1 = '/projects/ps-gymreklab/bonnieh/lrt/results/' + outFolder1 + str(per) + '_' + str(opt_allele) 
    outFile2 = '/projects/ps-gymreklab/bonnieh/lrt/results/' + outFolder2 + str(per) + '_' + str(opt_allele) 
   
    outFile1a = outFile1 + '_het.txt'
    outFile1b = outFile1 + '_common.txt'
    outFile1c = outFile1 + '_bins_narrow.txt'
    outFile1d = outFile1 + '_freqs.txt'
   
    outFile2a = outFile2 + '_het.txt'
    outFile2b = outFile2 + '_common.txt'
    outFile2c = outFile2 + '_bins_narrow.txt'
    outFile2d = outFile2 + '_freqs.txt'
    
    results1a = open(outFile1a, "w")
    results1b = open(outFile1b, "w")
    results1c = open(outFile1c, "w")
    results1d = open(outFile1d, "w")
    
    results2a = open(outFile2a, "w")
    results2b = open(outFile2b, "w")
    results2c = open(outFile2c, "w")
    results2d = open(outFile2d, "w")
    
    results_list = []
    results_list.append(results1a)
    results_list.append(results1b)
    results_list.append(results1c)
    results_list.append(results1d)
    
    results_list.append(results2a)
    results_list.append(results2b)
    results_list.append(results2c)
    results_list.append(results2d)
    
    # Write results header
    results1a.write("s" + "\t" + "het" +"\n")
    results1b.write("s" + "\t" + "common" +"\n")
    results1c.write("s" + "\t" + "bins_narrow" +"\n")
    results1d.write("s" + "\t" + "freqs" +"\n")
    
    results2a.write("s" + "\t" + "het" +"\n")
    results2b.write("s" + "\t" + "common" +"\n")
    results2c.write("s" + "\t" + "bins_narrow" +"\n")
    results2d.write("s" + "\t" + "freqs" +"\n")
    
    # Period info
    period_info = {}
    
    L2_log = 0.15 
    L3_log = 0.65 
    L4_log = 0.45 

    # mu, beta, p, l
    period_info[2] = [10**-5, 0.3, 0.6, L2_log, 6]
    period_info[3] = [10**-7, 0.3, 0.9, L3_log, 5]
    period_info[4] = [10**-6, 0.3, 0.9, L4_log, 3]

    num_alleles = 25
    N_e = 7300
    max_iter = 26000
    end_samp_n = 8000
    set_start_equal = False
    
    log_mu_prime = np.log10(period_info[per][0])+period_info[per][3]*(opt_allele - period_info[per][4])
    mu_prime = 10**log_mu_prime
    if mu_prime < 10**-8: mu_prime = 10**-8
    if mu_prime > 10**-3: mu_prime = 10**-3
            
    mu = mu_prime
    beta = period_info[per][1]
    p = period_info[per][2]
    L = period_info[per][3]
    
    for s in s_list:
        for result_file in results_list:
            result_file.write(str(s) + "\t")
        
        for i in range(0, LRT_num_sims):
            
            # Simulate allele frequencies
            allele_freqs_const, allele_freqs_euro = Simulate(num_alleles, N_e, mu, beta, p, L, s, max_iter, end_samp_n, set_start_equal)

            # Compute summary statistic of simulated allele frequencies
            het_const = 1-sum([item**2 for item in allele_freqs_const]) 
            het_euro = 1-sum([item**2 for item in allele_freqs_euro]) 

            common_const = (allele_freqs_const>=0.05).sum()
            common_euro = (allele_freqs_euro>=0.05).sum()
            
            bins_const_narrow = GetBinsNarrow(allele_freqs_const)
            bins_euro_narrow = GetBinsNarrow(allele_freqs_euro)
            
            results1a.write(str(het_const))
            results1b.write(str(common_const))
            results1c.write(','.join(str(item) for item in bins_const_narrow))
            results1d.write(','.join(str(item) for item in allele_freqs_const))
            
            results2a.write(str(het_euro))
            results2b.write(str(common_euro))
            results2c.write(','.join(str(item) for item in bins_euro_narrow))
            results2d.write(','.join(str(item) for item in allele_freqs_euro))
            
            if i == LRT_num_sims - 1:
                for result_file in results_list:
                    result_file.write('\n')
            else:
                for result_file in results_list:
                    result_file.write(';')
                
    for result_file in results_list:    
        result_file.close()
    
if __name__ == '__main__':
    main()