# Running SISTR Joint  

## SISTR Joint Script: SISTR_joint_method.py

## SISTR Joint Usage  
Example command:
```
python SISTR_joint_method.py \
     --period 2 \
     --sim-model eurodem_prior2_dinuc_a \
     --out-subfolder motif_analysis \
     --opt-allele-to-use 11 \
     --motif-to-use AT \
     --file-to-use allele_freqs_test.txt \
```

Required parameters:  
* __`--period <int>`__ which period to jointly estimate s on (e.g. either `2`, `3`, or `4`) 
* __`--sim-model <string>`__ mutation model from which to obtain s values and summary statistics (e.g. `eurodem_prior2_dinuc_a`) OR use `all` to run all mutation models for a given period  
* __`--out-subfolder <string>`__ name of output subfolder  
* __`--opt-allele-to-use <int>`__ optimal allele to jointly estimate s on  
* __`--motif-to-use <string>`__ "motif" to jointly estimate s on (does not necessarily need to be a repeat motif (e.g. AT), this parameter indicates which value of the column specified should be used to determine whether to include a particular locus in the joint analysis)  
* __`--file-to-use <string>`__ name of input file containing allele frequency data (see below for format specifications)    

Optional parameters:  
* __`--num-sims <int>`__ number of simulations to perform (default value: 50000) 
* __`--num-loci <int>`__ maximum number of loci for which to jointly estimate s on; if the number of loci is less than this value, all loci will be used  (default value: 1000)
* __`--motif-format`__ whether allele frequencies of the input data is in motif format  
* __`--header`__ whether input file has a header  
* __`--column <int>`__ column of input file to use to determine which loci should be included in the joint analysis (default value: 6)    
* __`--normalized`__ whether to use normalized values when comparing two heterozygosity distributions  
* __`--normalized-constant <float>`__ normalization constant if using normalized values when comparing 2 heterozygosity distributions (default value: 0.1)   
* __`--ks-test`__ whether to use Kolmogorov-Smirnov test to compare two heterozygosity distributions (default option is taking mean of differences of the two distributions)
* __`--top-x`__ whether to use the top x percent method (default option is using arbitrary thresholds to decide whether two heterozyosity distributions are sufficiently similar)
* __`--perc-acc <float>`__ what percent of top similar heterozygosity distributions to use (default value: 1)

## Input file format:
The input file should be a tab delimited file with 7 columns: chrom, start, end, allele_freqs, total, period, motif    
Note: The 7th column does not have to be a motif. It can be any annotation you want to stratify the loci by.  
Each row should be a separate locus.    
See example input file `allele_freqs_test.txt`.  

### Explanation of each column:
* chrom: chromosome number  
* start: starting position of locus in a reference genome  
* end: ending position of locus in a reference genome  
* allele_freqs: allele frequency data (see note below for required format)
* total: sample size (number of alleles) of empirical allele frequencies  
* period: number of base pairs in the repeat motif (e.g. period 2 refers to dinucleotides)  
* motif: repeat motif (or any other annotation)  

Note: Required format for allele frequency data  
Data for each allele present in the population is represented by two numbers separated by a colon. The first number is the allele size in base pairs relative to a reference allele size (which is calculated as `end - start + 1`). The second number is the number of copies of the allele in the population. Each allele is separated by a comma.  
   
Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|-2:100,0:6000,4:400|6500|2|AT| 
   
At this locus, the reference allele length is 26 base pairs or 13 repeat units. In the population, there are 3 alleles present: 100 copies of 12 repeat units, 6000 copies of 13 repeat units, and 400 copies of 15 repeat units.

## Output file format:
Num sims: XXX  
Number of loci used: XXX Total number loci: XXX  
Column: XXX "Motif" of column: XXX Optimal allele: XXX  
Observed heterozygosity distribution: [XXX, XXX, ... , XXX]  
Ignore this line  
Ignore this line  
Ignore this line  
Ignore this line  
Ignore this line  
Number a,b pairs accepted: XXX  
Sorted by mean - Median a, b:  
XXX,XXX  
Median mean: XXX  
2.5 percentile mean: XXX  
97.5 percentile mean: XXX  
Accepted params: (XXX, XXX), (XXX, XXX), ... , (XXX, XXX)  
Ignore the rest of the file  
See example output file `test_results.txt`.  