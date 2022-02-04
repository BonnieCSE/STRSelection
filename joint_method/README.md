# SISTR Joint   

SISTR Joint is an extension of SISTR (https://github.com/BonnieCSE/SISTR) that enables joint estimation of the distribution of selection coefficients across a set of short tandem repeats.  

It takes allele frequency data (per-locus frequencies of each allele length) for a particular class of loci as input. It assumes s for each STR is drawn from a gamma distribution with parameters (a, b) and outputs a posterior estimate of (a, b).  

# Installation
SISTR Joint uses Python3 and the following libraries in addition to the Python Standard Library: SciPy, NumPy. 

You can obtain SISTR Joint by cloning the Github repository:

```
git clone https://github.com/BonnieCSE/STRSelection
cd joint_method
```

## Note: Before running SISTR Joint, lookup tables must first be generated. 
**We recommend using the available lookup tables found [here](https://drive.google.com/drive/folders/1fjEWxmvQ38p397-9NUCWYAUoqpto7N_n?usp=sharing) by downloading the entire folder `sistr_joint_resources/` and saving it at the same directory where you are running SISTR Joint.**  

## SISTR Joint Usage  
Example command:
```
python SISTR_Joint_v1.py \
     --in-file allele_freqs_test.txt \
     --out-folder test_results/ \
     --period 2 \
     --sim-model eurodem_prior2_dinuc_e_1kg_euro \
     --motif-to-use AG \
     --opt-allele-to-use 20 \
     --motif-format \
     --top-x 
     
python SISTR_Joint_v1.py \
     --in-file allele_freqs_test.txt \
     --out-folder test_results/ \
     --period 2 \
     --sim-model all \
     --motif-to-use AG \
     --opt-allele-to-use 20 \
     --motif-format \
     --top-x
```

Required parameters:  
* __`--in-file <string>`__ name of input file containing allele frequency data (see below for format specifications)  
* __`--out-folder <string>`__ name of output folder  
* __`--period <int>`__ which period to jointly estimate s on (e.g. either `2`, `3`, or `4`)  
* __`--sim-model <string>`__ mutation model from which to obtain s values and summary statistics (e.g. `eurodem_prior2_dinuc_a_1kg_euro` OR use `all` to run all mutation models for a given period)  
* __`--motif-to-use <string>`__ "motif" to jointly estimate s on (does not necessarily need to be a repeat motif (e.g. AT), this parameter indicates which value of the column specified should be used to determine whether to include a particular locus in the joint analysis)  
* __`--opt-allele-to-use <int>`__ optimal allele to jointly estimate s on  
* __`--top-x`__ use the top x percent method to decide whether two heterozyosity distributions are sufficiently similar  

Optional parameters:  
* __`--lookup-table-folder <string>`__ name of folder with lookup tables of simulations (default value: sistr_joint_resources/)
* __`--num-sims <int>`__ number of simulations to perform (default value: 50000) 
* __`--num-loci <int>`__ maximum number of loci for which to jointly estimate s on; if the number of loci is less than this value, all loci will be used  (default value: 1000)
* __`--motif-format`__ whether allele frequencies of the input data is in motif format  
* __`--header`__ whether input file has a header  
* __`--column <int>`__ column of input file to use to determine which loci should be included in the joint analysis. Columns are 0-indexed, with the first column of the file as column 0 (default value: 6)  
* __`--perc-acc <float>`__ what percent of top similar heterozygosity distributions to use (default value: 1)

## Input file format:
The input file should be a tab delimited file with 7 columns: chrom, start, end, allele_freqs, total, period, motif    
Each row should be a separate locus.  
See example input file `allele_freqs_test.txt`. 

Notes:  
* The 7th column does not have to be a motif. It can be any annotation you want to stratify the loci by (e.g. 1.0 for coding and 0.0 for noncoding).  
* The file can contain loci of different classes. SISTR Joint will extract the relavant loci since as input, it requires the period, motif, and optimal allele of the class of loci to infer the distribution of selection coefficients on.

### Explanation of each column:
* chrom: chromosome number  
* start: starting position of locus in a reference genome  
* end: ending position of locus in a reference genome  
* allele_freqs: allele frequency data (see note below for required format)
* total: sample size (number of alleles) of empirical allele frequencies  
* period: number of base pairs in the repeat motif (e.g. period 2 refers to dinucleotides)  
* motif: repeat motif (or any other annotation)  

#### Note: SISTR Joint supports two different formats for the allele frequency data: motif format and alternate format. 

**Motif format:**   
Data for each allele present in the population is represented by the STR allele and the number of copies of the allele in the population separated by a colon. Each allele is separated by a comma.  

Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|ATATATATATATATATATATATAT:100,ATATATATATATATATATATATATAT:6000,ATATATATATATATATATATATATATAT:400|6500|2|AT| 

**Alternate format:**  

Data for each allele present in the population is represented by two numbers separated by a colon. The first number is the allele size in base pairs relative to a reference allele size (which is calculated as `end - start + 1`). The second number is the number of copies of the allele in the population. Each allele is separated by a comma.  
   
Example: 
| chrom | start | end | allele_freqs | total | period | motif |
|---|---|---|---|---|---|---|
|1|1000|1025|-2:100,0:6000,4:400|6500|2|AT| 
   
At this locus, the reference allele length is 26 base pairs or 13 repeat units. In the population, there are 3 alleles present: 100 copies of 12 repeat units, 6000 copies of 13 repeat units, and 400 copies of 15 repeat units.  

## Output file format:  
The output file is a .txt file with 11 lines.  
File name: per_\[period\]\_\[motif-to-use\]\_\[sim-model\]\_col_\[column\]\_opt_\[opt-allele\].txt
See below for an explanation of each line.  
See example output files generated by the example commands in the `test_results` folder.  

### Explanation of each line
Num sims: XXX  
Number of loci used: XXX Total number loci: XXX  
Column: XXX "Motif" of column: XXX Optimal allele: XXX  
Observed heterozygosity distribution: [XXX, XXX, ... , XXX]  
Number a,b pairs accepted: XXX  
Sorted by mean - Median a, b:  
XXX,XXX  
Median mean: XXX  
2.5 percentile mean: XXX  
97.5 percentile mean: XXX  
Accepted params: (XXX, XXX), (XXX, XXX), ... , (XXX, XXX)  

**Note: If `all` is used for the `--sim-model` option, the resulting output file will only have 6 lines. Lines 6-10 above will not be included.**