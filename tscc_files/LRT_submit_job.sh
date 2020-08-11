#!/bin/bash
file="/projects/ps-gymreklab/bonnieh/lrt/LRT_combos_per3_opt13.txt"
while IFS=$'\t' read -r -a line; do
    PERIOD=${line[0]}
    OPTIMAL_ALLELE=${line[1]}
    FILENUM=${line[2]}
    S_VALS=${line[3]}
    JOBNAME="LRT_${PERIOD}_${OPTIMAL_ALLELE}_${FILENUM}"
    NUM_SIMS=2000
    
    echo Running LRT for $JOBNAME 
    cat LRT_job.sh | sed "s/jobname/${JOBNAME}/g" | sed "s#per#${PERIOD}#g" | sed "s#opt#${OPTIMAL_ALLELE}#g" | sed "s#filenum#${FILENUM}#g" |sed "s#s_vals#${S_VALS}#g" | sed "s#num_sims#${NUM_SIMS}#g" | qsub
done < "$file"
