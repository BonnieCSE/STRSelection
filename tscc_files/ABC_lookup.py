import sys
sys.path.append("/projects/ps-gymreklab/bonnieh/helper_functions")
from Simulation_functions import *
from ABC_functions import *

def main():
    # Load parameters
    per = int(sys.argv[1]) 
    opt_allele = int(sys.argv[2])
    ABC_num_sims = int(sys.argv[3])
    k = float(sys.argv[4])
    theta = float(sys.argv[5])
    filenum = sys.argv[6]
    outFolder = sys.argv[7]
    outFile1 = '/projects/ps-gymreklab/bonnieh/abc/results/20k_' + outFolder + '/' + str(per) + '_' + str(opt_allele) + '_' + str(filenum) + '.txt'
    outFile2 = '/projects/ps-gymreklab/bonnieh/abc/results/50k_' + outFolder + '/' + str(per) + '_' + str(opt_allele) + '_' + str(filenum) + '.txt'
    outFile3 = '/projects/ps-gymreklab/bonnieh/abc/results/eurodem_' + outFolder + '/' + str(per) + '_' + str(opt_allele) + '_' + str(filenum) + '.txt'
    results1 = open(outFile1, "w")
    results2 = open(outFile2, "w")
    results3 = open(outFile3, "w")
    
    # Write results header
    results1.write("s" + "\t" + "het" + "\t" + "common" +"\t" + "bins" + "\t" + "freqs" + "\n")
    results2.write("s" + "\t" + "het" + "\t" + "common" +"\t" + "bins" + "\t" + "freqs" + "\n")
    results3.write("s" + "\t" + "het" + "\t" + "common" +"\t" + "bins" + "\t" + "freqs" + "\n")

    # Period info
    period_info = {}

    L1_log = 0.04
    L2_log = 0.15 
    L3_log = 0.65 
    L4_log = 0.45 

    # mu, beta, p, l
    period_info[1] = [10**-4.2, 0.5, 1, L1_log, 13]
    period_info[2] = [10**-5, 0.3, 0.6, L2_log, 6]
    period_info[3] = [10**-7, 0.3, 0.9, L3_log, 5]
    period_info[4] = [10**-6, 0.3, 0.9, L4_log, 3]

    number_alleles = 25
    n_effec = 7300
    max_iter = 55920
    end_samp_n = 8000
    
    log_mu_prime = np.log10(period_info[per][0])+period_info[per][3]*(opt_allele - period_info[per][4])
    mu_prime = 10**log_mu_prime
    if mu_prime < 10**-8: mu_prime = 10**-8 # 10**-8
    if mu_prime > 10**-3: mu_prime = 10**-3

    mu = mu_prime
    beta = period_info[per][1]
    p_param = period_info[per][2]
    L = period_info[per][3]

    for i in range(0, ABC_num_sims):

        # Draw s from prior
        # Assume only one possible prior value for other parameters (mu, beta, p, l)
        s = 0
        if k != -1:
            s = np.random.gamma(k, theta)
        else:
            s = np.random.uniform()
        if s > 1:
            s = 1

        # Simulate allele frequencies
        allele_freqs_20k, allele_freqs_50k, allele_freqs_euro = Simulate(number_alleles, n_effec, mu, beta, p_param, L, s, max_iter, end_samp_n)

        # Compute summary statistic of simulated allele frequencies
        het_20k = 1-sum([item**2 for item in allele_freqs_20k])
        het_50k = 1-sum([item**2 for item in allele_freqs_50k])
        het_euro = 1-sum([item**2 for item in allele_freqs_euro]) 

        common_20k = (allele_freqs_20k>=0.05).sum()
        common_50k = (allele_freqs_50k>=0.05).sum()
        common_euro = (allele_freqs_euro>=0.05).sum()
            
        bins_20k = GetBins(allele_freqs_20k, 5)
        bins_50k = GetBins(allele_freqs_50k, 5)
        bins_euro = GetBins(allele_freqs_euro, 5)
            
        # Write summary statistics and allele frequencies to file
        results1.write(str(s) + "\t" + str(het_20k) + "\t" + str(common_20k) + "\t" + ','.join(str(item) for item in bins_20k) + "\t" + ','.join(str(item) for item in allele_freqs_20k) + "\n")

        results2.write(str(s) + "\t" + str(het_50k) + "\t" + str(common_50k) + "\t" + ','.join(str(item) for item in bins_50k) + "\t" + ','.join(str(item) for item in allele_freqs_50k) + "\n")
        
        results3.write(str(s) + "\t" + str(het_euro) + "\t" + str(common_euro) + "\t" + ','.join(str(item) for item in bins_euro) + "\t" + ','.join(str(item) for item in allele_freqs_euro) + "\n")
        
    results1.close()
    results2.close()
    results3.close()

if __name__ == '__main__':
    main()