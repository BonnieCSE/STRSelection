# Running SISTR with Motif-specific Mutation Rates

## Note: Before running SISTR, lookup tables must first be generated. 
**When running SISTR with motif-specific mutation rates, we recommend using the available lookup tables found [here](https://drive.google.com/drive/folders/1mjZeKIYVa6FrtDS5LU4Hu-8lzUZJh3C2?usp=sharing) by downloading the entire folder `sistr_resources_motif_specific/` and saving it at the same directory where you are running SISTR.**  

## SISTR Motif-specific Usage
Example command:
```
# Using motif format for input allele frequencies
python SISTR_motif_specific_mut_rates.py \
     --in-file allele_freqs_test_motif_format.txt \
     --out-file test_results_motif_format.txt \
     --motif-format \
     --lrt-num-sims 200
```  

## See Main SISTR per-locus repo for more information about parameters, input file format, and output file format: https://github.com/BonnieCSE/SISTR/tree/master/sistr  

## 1000 Genomes Per-population Preliminary Analysis  

See `SISTR_results_v2` folder  